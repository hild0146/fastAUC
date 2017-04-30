#' Calculate and compare AUC estimate for clustered data
#'
#' This function can be used to calculate and compare AUC estimates. If a single
#' set of test measurements is provided, the corresponding AUC and variance
#' values are computed. If two sets of test measurements are provided, a
#' significance test for a difference in AUC estimates is also computed.
#' Optionally, the previous can also be computed in the presence of clustered
#' data per Obuchowski (1997). This function implicity assumes that higher
#' values of test_1 and test_2 are correlated with the positive gold standard
#' (status = TRUE)
#'
#'
#' @param test_1 the set of test measurments for which the AUC should be
#'   computed
#'
#' @param test_2 a second set of test measurments for which the AUC is computed
#'   and compared to the first set of test measurements
#'
#' @param status a logical vector representing the gold standard
#'   negative/positive status of an observaiton
#'
#' @param cluster variable identifying clusters
#'
#' @return a list containing \item{auc}{a vector of AUC estimates}
#'   \item{var}{the individual variance estimate or covariance matrix}
#'   \item{test_stat}{the test statistic for comparing the AUC estimates,
#'   distributed standard normal} \item{p_value}{the p-value associated with the
#'   AUC test comparison (2-sided)}
#'
#' @references Obuchowski N.A. (1997) Nonparametric Analysis of clustered ROC
#'   Curve Data. Biometrics, 53, 567-578.
#'
#'   DeLong E.R., DeLong D.M, Clarke-Pearson D.L. (1988) Comparing the areas
#'   under two or more correlated receiver operating characteristic curves: a
#'   nonparametric approach. Biometrics, 44, 837-845.
#'
#'
#' @examples
#'
#' library(fastAUC)
#'
#' n_obs <- 10000
#' n_clust <- 100
#'
#' cluster <- sample(x = 1:n_clust, size = n_obs, replace = T)
#' status <- as.logical(rbinom(n = n_obs, size = 1,  prob = 0.3))
#' test_1 <- ifelse(status,
#'                  rnorm(n = n_obs, mean = 15, sd = 20),
#'                  rnorm(n = n_obs,mean = 10, sd = 20))
#' test_2 <- ifelse(status,
#'                  rnorm(n = n_obs, mean = 20, sd = 20),
#'                  rnorm(n = n_obs,mean = 10, sd = 20))
#'
#' auc(test_1 = test_1, status = status)
#'
#' auc(test_1 = test_1, test_2 = test_2, status = status, cluster = cluster)
#'
#' @export auc


auc <- function(test_1,
                test_2 = NULL,
                status,
                cluster = NULL){


  n_test <- ifelse(is.null(test_2), 1, 2)

  # -- define dataset of interest

  dta <- data_frame(test_1, status)

  if(!is.null(test_1)) dta$test_2 <- test_2
  if(!is.null(cluster)) dta$cluster <- cluster

  dta <- na.omit(dta)
  dta$obs_num <- 1:nrow(dta)

  n_pos <- sum(dta$status == T)
  n_neg <- sum(dta$status == F)

  if(!is.null(cluster)){
    n_cluster_tot <- n_distinct(dta$cluster)
    n_cluster_pos <- n_distinct(dta$cluster[dta$status == T])
    n_cluster_neg <- n_distinct(dta$cluster[dta$status == F])
  }


  # -- compute auc & variance estimates for each test
  for(test_i in 1:n_test) {
    if(test_i == 1) {
      out <- data_frame('test' = (1:2)[1:n_test],
                        'auc' = NA,
                        'var' = NA)
    }

    if(test_i == 1) dta$test_i <- dta$test_1 else
      dta$test_i <- dta$test_2

    # must be sorted prior to calling C++ code
    dta_i <- arrange(dta, test_i)

    dta_i$kern_sum <- auc_kernel(test = dta_i$test_i,
                                 status = dta_i$status)$out

    auc <- sum(dta_i$kern_sum[dta_i$status == 1]) / n_pos / n_neg

    out$auc[test_i] <- auc


    if(is.null(cluster)){

      v_values <- dta_i %>%
        arrange(obs_num) %>%
        mutate(v = ifelse(status == T,
                          kern_sum / n_neg,
                          kern_sum / n_pos))

      s_values <- v_values %>%
        group_by(status) %>%
        summarize(s = sum( (v - auc) ^ 2))  %>%
        ungroup %>%
        mutate(s = ifelse(status == T,
                          s / (n_pos - 1) / n_pos,
                          s / (n_neg - 1) / n_neg))

      out$var[test_i] <- sum(s_values$s)


    } else{

      v_values <- dta_i %>%
        group_by(status, cluster) %>%
        summarize(cnt = n(),
                  kern_sum = sum(kern_sum)) %>%
        ungroup %>%
        mutate(v = ifelse(status == T,
                          kern_sum / n_neg,
                          kern_sum / n_pos))

      s_values <- v_values %>%
        group_by(status) %>%
        summarize(s = sum( (v - cnt * auc) ^ 2))  %>%
        ungroup %>%
        mutate(s = ifelse(status == T,
                          n_cluster_pos * s / (n_cluster_pos - 1) / n_pos / n_pos,
                          n_cluster_neg * s / (n_cluster_neg - 1) / n_neg / n_neg))

      s_cross <- v_values %>%
        mutate(s = v - cnt * auc)

      s_11 <- inner_join(select(filter(s_cross, status == T), cluster, s),
                         select(filter(s_cross, status == F), cluster, s),
                         by = 'cluster',
                         suffix = c('_pos', '_neg')) %>%
        summarize(s = sum(s_pos * s_neg) * n_cluster_tot / (n_cluster_tot - 1)) %>%
        as.numeric

      out$var[test_i] <- sum(s_values$s) + 2 * s_11 / n_pos / n_neg

      rm(s_values, s_cross)

    }

    # retain values for covariance calculation
    if(n_test == 2){
      if(test_i == 1) v_values_1 <- v_values else
        v_values_2 <- v_values
    }
  }


  # -- covariance between auc estimates
  if(n_test == 2){
    if(is.null(cluster)){

      v_values_1 <- v_values_1 %>%
        mutate(val_auc_1 = v - out$auc[1],
               val_auc_2 = v - out$auc[2])

      v_values_2 <- v_values_2 %>%
        mutate(val_auc_1 = v - out$auc[1],
               val_auc_2 = v - out$auc[2])

      s_10_12 <- 1 / (n_pos - 1) *
        sum(select(filter(v_values_1, status == T), val_auc_1) *
              select(filter(v_values_2, status == T), val_auc_2))

      s_01_12 <- 1 / (n_neg - 1) *
        sum(select(filter(v_values_1, status == F), val_auc_1) *
              select(filter(v_values_2, status == F), val_auc_2))

      cov_auc <- s_10_12 / n_pos + s_01_12 / n_neg

    } else{

      v_values_1 <- v_values_1 %>%
        mutate(val_auc_1 = v - cnt * out$auc[1],
               val_auc_2 = v - cnt * out$auc[2])

      v_values_2 <- v_values_2 %>%
        mutate(val_auc_1 = v - cnt * out$auc[1],
               val_auc_2 = v - cnt * out$auc[2])

      s_10_12 <- n_cluster_pos / ((n_cluster_pos - 1) * n_pos) *
        sum(select(filter(v_values_1, status == T), val_auc_1)  *
              select(filter(v_values_2, status == T), val_auc_2))

      s_01_12 <- n_cluster_neg / ((n_cluster_neg - 1) *  n_neg) *
        sum(select(filter(v_values_1, status == F), val_auc_1) *
              select(filter(v_values_2, status == F), val_auc_2))

      s_11_12 <- n_cluster_tot / (n_cluster_tot - 1) *
        inner_join(select(filter(v_values_1, status == T), cluster, val_auc_1),
                   select(filter(v_values_2, status == F), cluster, val_auc_2),
                   by = 'cluster') %>%
        summarize(sum(val_auc_1 * val_auc_2)) %>%
        as.numeric

      s_11_21 <- n_cluster_tot / (n_cluster_tot - 1) *
        inner_join(select(filter(v_values_2, status == T), cluster, val_auc_2),
                   select(filter(v_values_1, status == F), cluster, val_auc_1),
                   by = 'cluster') %>%
        summarize(sum(val_auc_1 * val_auc_2)) %>%
        as.numeric

      cov_auc <- s_10_12 / n_pos + s_01_12 / n_neg + (s_11_12 + s_11_21) / n_pos / n_neg

    }


    test_stat <- abs(out$auc[1] - out$auc[2]) / sqrt(sum(out$var) - 2 * cov_auc)
    p_value <- 2 * pnorm(test_stat, lower.tail = F)

    out_return <- list(auc = out$auc,
                       var = matrix(c(out$var[1],
                                      rep(cov_auc, 2),
                                      out$var[2]),
                                    nrow = 2),
                       test_stat = test_stat,
                     p_value = p_value)

    return(out_return)

  } else return(list(auc = out$auc, var = out$var))

}
