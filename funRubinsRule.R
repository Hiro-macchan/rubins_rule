funRubinsRule <- function(estimate_list, var_list,df=NULL,names_list =NULL){
  #入力値のエラーチェック
  require(mice)
  if(!is.list(estimate_list)|!is.list(var_list)){stop("An estimate_list,var_list must be list.")}
  if(!is.null(df)){
    if(length(df) != 1 & !is.vector(df)){stop("The df must be one length vector.")}
  }
  
  n_imp <- length(estimate_list)
  n_est_variable <- unique(sapply(estimate_list,length))
  n_var_variable <- unique(sapply(var_list,length))
  if(n_est_variable != n_var_variable){stop("Error1")}
  estimate_matrix <- matrix(,n_imp,n_est_variable)
  for(i in 1:n_imp){
    estimate_matrix[i,] <- estimate_list[[i]]
  }
  var_matrix <- matrix(,n_imp,n_var_variable)
  for(i in 1:n_imp){
    var_matrix[i,] <- var_list[[i]]
  }
  # Rubins Rule で統合
  #estimate の統合
  pooled_estimate <- apply(estimate_matrix,MARGIN = 2,mean)
  conbind_est <- pooled_estimate
  # variance の統合（u_var とb_var を計算して足す。t_var = u_var + (1+1/n_imp)*b_var
  u_var <- apply(var_matrix,MARGIN = 2,mean)
  #var_est <- apply(var_matrix,MARGIN = 2,var)
  #b_var <- (var_est)
  e <- estimate_matrix - matrix(pooled_estimate,nrow = n_imp, ncol = n_est_variable,byrow = T)
  b_var <- diag(((n_imp-1)^(-1))*(t(e) %*% e))
  t_var <- u_var + (1+1/n_imp)*b_var
  conbind_var <- t_var
  # 自由度の計算
  # もしinput される自由度がなければOld Rule で計算
  # df がinput されていれば補正する。
  r <- (1 + 1/n_imp) * (b_var/u_var)
  lambda <- (1 + 1/n_imp) * (b_var/t_var)
  mice.df <- function(m, lambda, dfcom, method) {
    if (is.null(dfcom)) {
      dfcom <- 999999
      warning("Large sample assumed.")
    }
    lambda[lambda < 1e-04] <- 1e-04
    dfold <- (m - 1)/lambda^2
    dfobs <- (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
    df <- dfold * dfobs/(dfold + dfobs)
    if (method != "smallsample") 
      df <- dfold  ## Rubin 1987, 3.1.6, Van Buuren 2012, 2.30, added 31/10/2012
    return(df)
  }
  df_mod <- mice.df(m=n_imp, lambda = lambda,dfcom = df, method = "smallsample")

  t_stat <- (conbind_est-0)/sqrt(conbind_var)
  p_value <- (1-pt(q = abs(t_stat),df = df_mod))*2
  upper <- conbind_est + qt(p = 0.975, df = df_mod) * sqrt(conbind_var)
  lower <- conbind_est + qt(p = 0.025, df = df_mod) * sqrt(conbind_var)
  return_df <- data.frame(est = sprintf(conbind_est,fmt = "%5.3f"),t_stat = sprintf(t_stat,fmt = "%5.2f"), 
                          p = sprintf(p_value,fmt = "%5.4f"),upper=sprintf(upper,fmt="%5.3f"), lower = sprintf(lower,fmt="%5.3f"))
  return_df$CI <- paste0(return_df$est,return_df$conbind_est, " (",return_df$upper,", ", return_df$lower, ")") 
  row.names(return_df) <- names_list
  return(list(return_df=return_df,est = conbind_est, var = conbind_var, t_stat = t_stat, p = p_value#, CI = 
                #paste0( c(conbind_est, " (",upper,", ", lower, ")") )
              ,df=df_mod,u_var = u_var,b_var=b_var,r=r,lambda=lambda))
  
}
