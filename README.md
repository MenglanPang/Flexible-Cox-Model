# Extension of the Flexible Cox Model
R code for estimating the baseline hazard from the extended flexible Cox model.

## Content
### Description
The flexible extension of the Cox model is built based on the flexible Cox model developped by Wynant W and Wang Y.


Reference:
"Wynant W, Abrahamowicz M. Impact of the model building strategy on the inference about time-dependent and non-linear covariate effects in survival analysis. *SStatistics in medicin* 2014 Aug;33(19):3318-37."


"Willy Wynant and Michal Abrahamowicz. Flexible estimation of survival curves conditional on non-linear and time-dependent predictor effects. *Statistics in medicine* 2016 ; 35(4):553–565."


The code has been written using R with the following version information:<br/>
- R version 3.6.3 (2020-02-29)<br/> 
- Platform x86_64-apple-darwin15.6.0 (64-bit)<br/> 
- Using R packages:<br/> 
  - survival version 3.1-12
  - splines version v3.6.2
  
#### Code to implement the extended flexible Cox model:
##### `FlexCox.R`
This program is orignally developped by Wynant and Wang, with small modification by Pang. 
It includes functions to provide estimates of:
- TD effects
- NL effects
- Hazard ratio when TD is not present


The program is called by the program `Example.R`. 

#### Code to estimate the baseline hazard function from the extended flexible Cox model:
##### `FlexCox_BH.R`
It includes functions to provide estimates of:
- Hazard function conditional on an arbitrary covariate pattern with their relevant TD and NL effects
- Survival curve conditional on an arbitrary covariate pattern with their relevant TD and NL effects


The program is called by the program `Example.R`. 

#### Code to run the flexible Cox model:
##### `Example.R`
The program read the data in `sepsis.csv`, and generates the results save in `sepsis_flexcox.rda`
 
For questions or comments about the code please contact Menglan Pang (menglan.pang at mail.mcgill.ca).
