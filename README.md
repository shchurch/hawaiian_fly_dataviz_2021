# hawaiian_fly_dataviz_2021

The data visualization is hosted publically at https://shchurch.shinyapps.io/hawaiian_fly_dataviz_2021
At that site you can access and query the dataset of gene expression evolution across Hawaiian flies

The data visualization accompanies our study of gene expression evolution in Hawaiian Drosophilidae flies,
as described in the following preprint: https://doi.org/10.1101/2021.11.30.470652
and located in the repostiory https://github.com/shchurch/hawaiian_drosophilidae_expression_2021.

This data visualization is an `R::shiny` app. It can be accessed locally through R, by following these steps:

## Dependencies
These dependencies must be installed:

```
R
install.packages("shiny")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("ggcorrplot")
install.packages("gridExtra")
install.packages("ape")
install.packages("viridis")
install.pacakges("png")
```

## Run the app
```
library(shiny)
runGitHub("hawaiian_fly_dataviz_2021","shchurch",ref="main")
```



