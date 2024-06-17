
library(beeca)
library(RobinCar)

beeca_wrapper <- function(glm_fit, method="Ge", type="model-based", mod=F, strata=NULL) {
  
  res <- glm_fit |>
    beeca::get_marginal_effect("trt", method=method, type=type, contrast="diff", reference = "0",
                               mod=mod, strata=strata)
  
  return(list(marginal_est = res$marginal_est[[1]],
              marginal_var = res$marginal_se[[1]] ^ 2))
  
}



robincar_wrapper <- function(glm_fit) {
  
  fit_formula <- glm_fit$formula
  
  res <- robincar_glm(data.frame(glm_fit$data),
                      response_col = "y", 
                      treat_col="trt", 
                      formula = fit_formula,
                      g_family = glm_fit$family,
                      contrast_h = "diff",
                      adj_method = "homogeneous")
  
  return(res$contrast$varcov[[1]])
  
}