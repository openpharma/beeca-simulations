# Simulations supporting paper: *'Estimating the Variance of Covariate-Adjusted Estimators of Average Treatment Effects in Clinical Trials with Binary Endpoints'*, Magirr et al (2024)

This repo contains the reproducible code for simulation results in the paper by [Magirr et al. (2024)](https://osf.io/9mp58/).
Simulations are based on variance estimators as implemented in the [{beeca}](https://github.com/openpharma/beeca) package.


#### Contents

-   `R/extreme_scenario.Rmd`: Results for section 4 ('Extreme scenario').

-   `R/realistic_scenario.Rmd`: Results for section 4 ('Realistic scenario').

-   `R/misspecified_scenario.Rmd`: Results for section 4 ('Model misspecification').

-   `R/plot_ge_vs_ye.Rmd`: Results for section 5 ('Assessing conservativeness').

-   **R/utils**: Helper functions for simulation machinery.

-   `sessionInfo.yaml`: snapshot of R session to aid reproducibility.

#### References

* Ge, Miaomiao, L Kathryn Durham, R Daniel Meyer, Wangang Xie, and Neal Thomas. 2011. "Covariate-Adjusted Difference in Proportions from Clinical Trials Using Logistic Regression and Weighted Risk Differences." *Drug Information Journal: DIJ/Drug Information Association* 45: 481--93. <https://doi.org/10.1177/009286151104500409>

* Magirr, Dominic, Mark Baillie, Craig Wang, and Alexander Przybylski. 2024. “Estimating the Variance of Covariate-Adjusted Estimators of Average Treatment Effects in Clinical Trials with Binary Endpoints.” OSF. May 16. <https://osf.io/9mp58>.

* Przybylski A, Baillie M, Wang C, Magirr D (2024). beeca: Binary Endpoint Estimation with Covariate Adjustment. R package version 0.1.2, <https://openpharma.github.io/beeca/>

* Ye, Ting, Marlena Bannick, Yanyao Yi, and Jun Shao. 2023. "Robust Variance Estimation for Covariate-Adjusted Unconditional Treatment Effect in Randomized Clinical Trials with Binary Outcomes." *Statistical Theory and Related Fields* 7 (2): 159--63. <https://doi.org/10.1080/24754269.2023.2205802>
