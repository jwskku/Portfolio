library(dplyr)
library(lubridate)
library(purrr)
library(zoo)
library(tidyr)
library(ggplot2)
library(dtw)
library(proxy)
library(tidyr)
library(dtwclust)   # DTW 클러스터링 패키지 추가
library(foreach)
library(doParallel)
setwd("C:/weather")
raw_data <- readxl::read_excel('train_heat.xlsx')
raw_test <- read.csv('test_heat.csv')
new_name <- colnames(raw_test)

# 변수명 정리
data <- raw_data %>% select(-...1) %>%setNames(new_name)

# 결측치 na로 설정
data <- data %>%
  mutate(across(everything(), ~ ifelse(.x == -99, NA_real_, .x)))

# 결측 수 확인
print(colSums(is.na(data)))

# 날짜로 변경
data <- data %>%
  mutate(TM = ymd_h(as.character(TM), tz = "Asia/Seoul"))

# branch_ID별 그룹화 후 시간대로 정리.
data <- data %>%
  group_by(branch_ID) %>%
  arrange(TM)

# 연속na길이 구하기 함수
calculate_total_na_streak_fixed <- function(x) {
  rle_x <- rle(is.na(x))
  
  # rep(벡터, 반복횟수) 형태로 사용. times 인자를 명시하지 않는 것이 핵심입니다.
  total_lengths <- rep(rle_x$lengths, rle_x$lengths)
  
  is_na_block <- rep(rle_x$values, rle_x$lengths)
  total_lengths[!is_na_block] <- 0
  
  return(total_lengths)
}
# na있는 변수 지정
na_cols <- c("TA", "WD", "WS", "RN_DAY", "RN_HR1", "HM", "SI", "ta_chi", "heat_demand")

data <- data %>%
  select(-ends_with("_na_streak_length"), -ends_with("_na_consecutive"))

list_by_branch <- data %>%
  group_split(branch_ID)

#지사별로 나눠서 na연속성 구하기
data_processed <- map_dfr(list_by_branch, function(df) {
  df %>%
    mutate(across(
      all_of(na_cols),
      calculate_total_na_streak_fixed,
      .names = "{.col}_na_streak_length"
    ))
})

data <- data_processed



## SI 조정하기

print("--- SI 변수의 NA 연속 길이 분포 확인 ---")
print(
  data %>%
    filter(SI_na_streak_length > 0) %>%
    count(SI_na_streak_length, sort = TRUE),n=28
)


# 18시 '이상(>=)' 이거나 8시 '미만(<)'인 NA를 0으로 대치
data <- data %>%
  mutate(SI = if_else(is.na(SI) & (hour(TM) >= 18 | hour(TM) < 8), 0, SI))







print("--- 보간 전 각 변수별 NA 개수 ---")
print(colSums(is.na(data[, na_cols])))


target_cols_no_wd <- c("TA","WS", "RN_HR1", "HM", "SI", "ta_chi", "heat_demand")

data <- data %>%
  group_by(branch_ID)

# 각 대상 변수에 대해 반복 작업 수행
for (col_name in target_cols_no_wd) {
  streak_col_name <- paste0(col_name, "_na_streak_length")
  
  # 임시로 전체 보간된 열을 생성
  interp_col_name <- paste0(col_name, "_interp")
  data <- data %>%
    mutate(!!interp_col_name := na.approx(!!sym(col_name), na.rm = FALSE))
  
  # streak 길이가 1 또는 2인 NA만 보간된 값으로 대치
  data <- data %>%
    mutate(
      !!col_name := if_else(
        is.na(!!sym(col_name)) & !!sym(streak_col_name) %in% c(1, 2),
        !!sym(interp_col_name),
        !!sym(col_name)
      )
    ) %>%
    # 임시 보간열 제거
    select(-all_of(interp_col_name))
}

# 작업 완료 후 그룹 해제
data <- data %>%
  ungroup()

# --- 보간 후 NA 개수 확인 ---
print("--- 보간 후 각 변수별 남은 NA 개수 (WD 제외) ---")
print(colSums(is.na(data[, na_cols])))
















print("--- DTW 기반 클러스터링 시작 (NA 시점 제거 로직) ---")

# --- 1. 클러스터링을 위한 데이터 준비 ---
# 각 branch_ID를 열(column)로, 시간을 행(row)으로 갖는 wide 포맷으로 데이터를 변환
heat_demand_wide <- data %>%
  select(TM, branch_ID, heat_demand) %>%
  pivot_wider(names_from = branch_ID, values_from = heat_demand) %>%
  arrange(TM) %>%
  select(-TM)

# --- 2. NA가 포함된 행(시점) 제거 ---
# 특정 시점에 한 지점이라도 NA가 있으면, 그 시점 전체를 분석에서 제외
heat_demand_wide_complete <- na.omit(heat_demand_wide)
print(paste("원본 시점 개수:", nrow(heat_demand_wide), 
            ", NA 제외 후 분석에 사용될 시점 개수:", nrow(heat_demand_wide_complete)))



# 행과 열을 바꿔 매트릭스 생성
heat_demand_matrix <- t(as.matrix(heat_demand_wide_complete))


# --- 3. DTW 거리 행렬 계산 및 계층적 클러스터링 (이전과 동일) ---
dtw_dist_matrix <- proxy::dist(heat_demand_matrix, method = "DTW")
hclust_result <- hclust(dtw_dist_matrix, method = "ward.D2")
plot(hclust_result) # 덴드로그램 확인

# --- 4. 클러스터 그룹 할당 (이전과 동일) ---
clusters_dtw <- cutree(hclust_result, k = 3)
branch_summary_dtw <- data.frame(
  branch_ID = names(clusters_dtw),
  cluster = clusters_dtw
)


# 새로운 DTW 클러스터 정보를 원본 데이터에 결합
data <- data %>%
  left_join(branch_summary_dtw, by = "branch_ID")

# 새로운 클러스터 할당 결과 확인
print("--- DTW 클러스터링 할당 결과 ---")
print(count(branch_summary_dtw, cluster))







# 보간 대상: WD를 제외한 숫자형 변수들
imputation_targets <- setdiff(target_cols_no_wd, "WD")

# --- 보간 전 NA 개수 확인 ---
print("--- 보간 전 각 변수별 NA 개수 ---")
print(colSums(is.na(data[, imputation_targets])))


# a. 각 시점(TM) & 클러스터(cluster)별 평균값을 미리 계산하여 임시열로 추가
data <- data %>%
  group_by(cluster, TM) %>%
  mutate(across(
    all_of(imputation_targets),
    ~mean(.x, na.rm = TRUE),
    .names = "{.col}_cluster_mean"
  )) %>%
  ungroup()




print(colSums(is.na(data[, imputation_targets])))


# --- 2. '첫 행' 결측치 우선 처리 (RN_HR1 추가) ---
data <- data %>%
  group_by(branch_ID) %>%
  mutate(
    row_num_in_group = row_number(),
    
    TA = if_else(row_num_in_group == 1 & is.na(TA) & TA_na_streak_length >= 3, TA_cluster_mean, TA),
    WS = if_else(row_num_in_group == 1 & is.na(WS) & WS_na_streak_length >= 3, WS_cluster_mean, WS),
    RN_HR1 = if_else(row_num_in_group == 1 & is.na(RN_HR1) & RN_HR1_na_streak_length >= 3, RN_HR1_cluster_mean, RN_HR1), # <-- 로직 추가된 부분
    HM = if_else(row_num_in_group == 1 & is.na(HM) & HM_na_streak_length >= 3, HM_cluster_mean, HM),
    SI = if_else(row_num_in_group == 1 & is.na(SI) & SI_na_streak_length >= 3, SI_cluster_mean, SI),
    ta_chi = if_else(row_num_in_group == 1 & is.na(ta_chi) & ta_chi_na_streak_length >= 3, ta_chi_cluster_mean, ta_chi),
    heat_demand = if_else(row_num_in_group == 1 & is.na(heat_demand) & heat_demand_na_streak_length >= 3, heat_demand_cluster_mean, heat_demand)
    
  ) %>%
  select(-row_num_in_group)







# --- 2. 병렬 처리 클러스터 설정 ---
# 사용 가능한 코어 수 확인 (시스템 안정을 위해 하나는 제외)
num_cores <- detectCores() - 2
# 클러스터 생성
cl <- makeCluster(num_cores)
# 병렬 처리 백엔드 등록
registerDoParallel(cl)

print(paste(num_cores, "개의 코어를 사용하여 병렬 처리를 시작합니다..."))


# --- 3. foreach를 이용한 병렬 반복문 실행 ---

list_by_branch <- data %>%
  group_split(branch_ID)

# b. foreach 반복문을 수정된 옵션으로 실행합니다.
data_processed <- foreach(
  branch_df = list_by_branch, 
  .combine = 'rbind',                   # <-- 'list'에서 'rbind'로 변경
  .export = 'imputation_targets',      # <-- 사용할 변수 명시적 전달
  .packages = 'dplyr'                  # <-- 사용할 패키지 명시적 전달
) %dopar% {
  
  # --- 각 코어에서 실행되는 로직 (이전과 동일) ---
  branch_indices <- 1:nrow(branch_df)
  
  if (length(branch_indices) > 1) {
    for (i in branch_indices[-1]) {
      for (col in imputation_targets) {
        streak_col <- paste0(col, "_na_streak_length")
        if (!is.na(branch_df[i, streak_col]) && branch_df[i, streak_col] >= 3 && is.na(branch_df[i, col])) {
          lag_val <- branch_df[[i - 1, col]]
          mean_col <- paste0(col, "_cluster_mean")
          cluster_mean_val <- branch_df[[i, mean_col]]
          imputed_value <- NA_real_
          if (!is.na(lag_val)) {
            if (!is.nan(cluster_mean_val)) {
              imputed_value <- (lag_val * 0.5) + (cluster_mean_val * 0.5)
            } else {
              imputed_value <- lag_val
            }
          } else {
            if (!is.nan(cluster_mean_val)) {
              imputed_value <- cluster_mean_val
            }
          }
          if (!is.na(imputed_value)) {
            branch_df[[i, col]] <- imputed_value
          }
        }
      }
    }
  }
  
  # 처리된 지점의 데이터프레임을 반환
  return(branch_df)
}

# --- 4. 클러스터 종료 및 결과 할당 ---
stopCluster(cl)

# foreach의 결과가 바로 최종 데이터프레임이 됩니다.
data <- data_processed


# 4. 임시 열들 제거
data <- data %>%
 select(-ends_with("_na_streak_length"),-ends_with("_cluster_mean"),-cluster)


# --- 보간 후 NA 개수 확인 ---
print("--- 병렬 처리 보간 후 각 변수별 남은 NA 개수 ---")
print(colSums(is.na(data[, imputation_targets])))











#--- wd 제외 채워진 데이터
write.csv(data,'train_temporary_dtw.csv',index=FALSE)

data <- read.csv('train_temporary_dtw.csv')







raw_data <- readxl::read_excel('train_heat.xlsx')
raw_test <- read.csv('test_heat.csv')

data <- read.csv('test_heat.csv')


# 결측치 na로 설정
data <- data %>%
  mutate(across(everything(), ~ ifelse(.x == -99, NA_real_, .x)))


# 결측 수 확인
print(colSums(is.na(data)))

# 날짜로 변경
data <- data %>%
  mutate(TM = ymd_h(as.character(TM), tz = "Asia/Seoul"))

# branch_ID별 그룹화 후 시간대로 정리.
data <- data %>%
  group_by(branch_ID) %>%
  arrange(TM)

# 연속na길이 구하기 함수
calculate_total_na_streak_fixed <- function(x) {
  rle_x <- rle(is.na(x))
  
  # rep(벡터, 반복횟수) 형태로 사용. times 인자를 명시하지 않는 것이 핵심입니다.
  total_lengths <- rep(rle_x$lengths, rle_x$lengths)
  
  is_na_block <- rep(rle_x$values, rle_x$lengths)
  total_lengths[!is_na_block] <- 0
  
  return(total_lengths)
}
# na있는 변수 지정
na_cols <- c("TA", "WD", "WS", "RN_DAY", "RN_HR1", "HM", "SI", "ta_chi")

data <- data %>%
  select(-ends_with("_na_streak_length"), -ends_with("_na_consecutive"))

list_by_branch <- data %>%
  group_split(branch_ID)

#지사별로 나눠서 na연속성 구하기
data_processed <- map_dfr(list_by_branch, function(df) {
  df %>%
    mutate(across(
      all_of(na_cols),
      calculate_total_na_streak_fixed,
      .names = "{.col}_na_streak_length"
    ))
})

data <- data_processed



## SI 조정하기

print("--- SI 변수의 NA 연속 길이 분포 확인 ---")
print(
  data %>%
    filter(SI_na_streak_length > 0) %>%
    count(SI_na_streak_length, sort = TRUE),n=28
)


# 18시 '이상(>=)' 이거나 8시 '미만(<)'인 NA를 0으로 대치
data <- data %>%
  mutate(SI = if_else(is.na(SI) & (hour(TM) >= 18 | hour(TM) < 8), 0, SI))





print("--- 보간 전 각 변수별 NA 개수 ---")
print(colSums(is.na(data[, na_cols])))


target_cols_no_wd <- c("TA","WD" ,"WS", "RN_HR1", "HM", "SI", "ta_chi")

data <- data %>%
  group_by(branch_ID)

# 각 대상 변수에 대해 반복 작업 수행
for (col_name in target_cols_no_wd) {
  streak_col_name <- paste0(col_name, "_na_streak_length")
  
  # 임시로 전체 보간된 열을 생성
  interp_col_name <- paste0(col_name, "_interp")
  data <- data %>%
    mutate(!!interp_col_name := na.approx(!!sym(col_name), na.rm = FALSE))
  
  # streak 길이가 1 또는 2인 NA만 보간된 값으로 대치
  data <- data %>%
    mutate(
      !!col_name := if_else(
        is.na(!!sym(col_name)) & !!sym(streak_col_name) %in% c(1, 2),
        !!sym(interp_col_name),
        !!sym(col_name)
      )
    ) %>%
    # 임시 보간열 제거
    select(-all_of(interp_col_name))
}



# 작업 완료 후 그룹 해제
data <- data %>%
  ungroup()

# --- 보간 후 NA 개수 확인 ---
print("--- 보간 후 각 변수별 남은 NA 개수 (WD 제외) ---")
print(colSums(is.na(data[, na_cols])))






# d. 클러스터 정보를 원본 데이터(data)에 결합
data <- data %>%
  left_join(branch_summary %>% select(branch_ID, cluster), by = "branch_ID")

# 클러스터링 결과 확인 (각 클러스터에 몇 개의 지점이 속했는지)
print("--- 클러스터링 할당 결과 ---")
print(count(branch_summary, cluster))


# 보간 대상: WD를 제외한 숫자형 변수들
imputation_targets <- setdiff(target_cols_no_wd, "WD")

# --- 보간 전 NA 개수 확인 ---
print("--- 보간 전 각 변수별 NA 개수 ---")
print(colSums(is.na(data[, imputation_targets])))


# a. 각 시점(TM) & 클러스터(cluster)별 평균값을 미리 계산하여 임시열로 추가
data <- data %>%
  group_by(cluster, TM) %>%
  mutate(across(
    all_of(imputation_targets),
    ~mean(.x, na.rm = TRUE),
    .names = "{.col}_cluster_mean"
  )) %>%
  ungroup()





# --- 2. '첫 행' 결측치 우선 처리 (RN_HR1 추가) ---
data <- data %>%
  group_by(branch_ID) %>%
  mutate(
    row_num_in_group = row_number(),
    
    TA = if_else(row_num_in_group == 1 & is.na(TA) & TA_na_streak_length >= 3, TA_cluster_mean, TA),
    WS = if_else(row_num_in_group == 1 & is.na(WS) & WS_na_streak_length >= 3, WS_cluster_mean, WS),
    RN_HR1 = if_else(row_num_in_group == 1 & is.na(RN_HR1) & RN_HR1_na_streak_length >= 3, RN_HR1_cluster_mean, RN_HR1), # <-- 로직 추가된 부분
    HM = if_else(row_num_in_group == 1 & is.na(HM) & HM_na_streak_length >= 3, HM_cluster_mean, HM),
    SI = if_else(row_num_in_group == 1 & is.na(SI) & SI_na_streak_length >= 3, SI_cluster_mean, SI),
    ta_chi = if_else(row_num_in_group == 1 & is.na(ta_chi) & ta_chi_na_streak_length >= 3, ta_chi_cluster_mean, ta_chi)
  ) %>%
  select(-row_num_in_group)

# --- 3. 나머지 행 순차적 보간 (기존 for loop) ---
print("--- 순차적 보간 실행 (시간이 소요될 수 있습니다) ---")
for (branch in unique(data$branch_ID)) {
  branch_indices <- which(data$branch_ID == branch)
  for (i in branch_indices[-1]) {
    for (col in imputation_targets) {
      # ... (이하 로직은 이전 답변과 동일)
      streak_col <- paste0(col, "_na_streak_length")
      if (!is.na(data[i, streak_col]) && data[i, streak_col] >= 3 && is.na(data[i, col])) {
        lag_val <- data[[i - 1, col]] 
        mean_col <- paste0(col, "_cluster_mean")
        cluster_mean_val <- data[[i, mean_col]]
        imputed_value <- NA_real_
        if (!is.na(lag_val)) {
          if (!is.nan(cluster_mean_val)) {
            imputed_value <- (lag_val * 0.5) + (cluster_mean_val * 0.5)
          } else {
            imputed_value <- lag_val
          }
        } else {
          if (!is.nan(cluster_mean_val)) {
            imputed_value <- cluster_mean_val
          }
        }
        if (!is.na(imputed_value)) {
          data[[i, col]] <- imputed_value
        }
      }
    }
  }
  print(paste(branch, "지점 처리 완료."))
}

# 4. 임시 열들 제거
data <- data %>%
  select(-ends_with("_na_streak_length"),-ends_with("_cluster_mean"),-cluster)


# --- 보간 후 NA 개수 확인 ---
print("--- 최종 보간 후 각 변수별 남은 NA 개수 ---")
print(colSums(is.na(data[, imputation_targets])))




# 최종적으로 남은 NA를 처리할 변수 목록
final_targets <- c("WS", "RN_HR1") 

print("--- 최종 마무리 전 남은 NA 개수 ---")
print(colSums(is.na(data[, final_targets])))

# 2단계 LOCF를 이용한 최종 결측치 처리
data <- data %>%
  group_by(branch_ID) %>%
  
  # 1단계: Last Observation Carried Forward (중간 NA 및 후미 NA 처리)
  mutate(across(all_of(final_targets), ~na.locf(.x, na.rm = FALSE))) %>%
  
  # 2단계: Next Observation Carried Backward (선두 NA 처리)
  mutate(across(all_of(final_targets), ~na.locf(.x, na.rm = FALSE, fromLast = TRUE))) %>%
  
  ungroup()

print("--- 최종 마무리 후 남은 NA 개수 ---")
print(colSums(is.na(data[, final_targets])))

print(colSums(is.na(data)))

#--- wd 제외 채워진 데이터
write.csv(data,'test_temporary.csv')

data <- read.csv('test_temporary.csv')








