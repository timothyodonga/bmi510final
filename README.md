# bmi510final
This project includes a package of R functions from the BMI510 final project. The following functios are included in the package.
- logLikBernoulli
- survCurv
- unscale
- pcApprox
- standardizeNames
- minimumN
- downloadRedcapReport

## Installation
To install the package from github, first install the devtools package in your R installation if you haven't already done using the command below
```{r}
install.packages("devtools")
```
To download the BMI510 final package to your R installation, run the command below
```{r}
devtools::install_github("https://github.com/timothyodonga/bmi510final")
```
In order to use functions from the package;
You can load the entire package
```{r}
library(bmi510final)
```
Or, use the individual functions from the package without loading the entire package. For example if you were to use the logLikBernoulli function from the package
```{r}
bmi510final::logLikBernoulli(c(1,0,0,1))
```

### Note
In order to run the downloadRedcapReport function, you need to create a .Renviron file and add the redcap token you have into this file.
```{text}
redcap_token=XXXX
``` 
Then load the environment variables  using the command below in your R console
```{r}
Sys.getenv()
```
If the token isn't visible you might need to restart your RStudio
