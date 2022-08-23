[![DOI](https://zenodo.org/badge/528005028.svg)](https://zenodo.org/badge/latestdoi/528005028)

# GcgrKO Shiny App

This app displays RNA sequencing data from livers of GcgrKO mice compared to wildtype littermates. If you wish to view the app online, it is accessible here:
  
https://weweralbrechtsenlab.shinyapps.io/GcgrKO/ 

If you would like to run the app from github, you will need a few R packages:

```{r}
install.packages(c(shiny, EnhancedVolcano, tidyverse, DT, tidyverse, shinyjs, writexl))
```

Once you have all the packages installed, simply run these lines in R. It will download the app and display it in a browser window:

```{r}
library(EnhancedVolcano)
library(tidyverse)
library(DT)
library(shiny)
library(shinyjs)
library(shinythemes) 
library(writexl)

runGitHub(rep = "GcgrKO_ShinyApp", username = "nicwin98", ref = "main")
```
