

# Calculate simulation result quantities of interest
post_sims <- function(sims_df, method){
  
  post_df <- sims_df[,c("scenario", "true_diff", "point_estimate")]
  post_df["var"] = sims_df[method]
  post_df["lower"] <- post_df[,"point_estimate"] - qnorm(0.975) * sqrt(sims_df[method])
  post_df["upper"] <- post_df[,"point_estimate"] + qnorm(0.975) * sqrt(sims_df[method])
  post_df["cover"] <- ifelse(post_df[,"lower"] < sims_df[,"true_diff"] & post_df[,"upper"] > sims_df[,"true_diff"], 1, 0)
  post_df["reject"] <- ifelse((post_df[,"lower"] > 0) | (post_df[,"upper"] < 0), 1, 0)
  
  return(cbind(post_df, method = method))
}




## Summarise simulation results
sum_sims <- function(post_sims_df){
  
  df_grouped <- dplyr::group_by(post_sims_df, scenario)
  
  df_sum <- dplyr::summarise(df_grouped,
                             av_se = mean(sqrt(var)),
                             coverage = mean(cover),
                             # power = mean(reject),
                             method = method[1])
  
  return(as.data.frame(df_sum))
}


