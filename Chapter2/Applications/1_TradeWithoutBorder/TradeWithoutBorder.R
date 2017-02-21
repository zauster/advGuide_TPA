##*******************************************************************************
##**************  CHAPTER 2 - GENERAL EQUILIBRIUM TRADE POLICY  ***************** 
##**************		      ANALYSIS WITH STRUCTURAL GRAVITY 	*********
##*******************************************************************************

##******************  APPLICATION 1: TRADE WITHOUT BORDER  ********************** 

## This application applies the methods developed by Anderson et al. (2015) in 
## order to investigates the potential effects of removing all international 
## borders, while preserving the effects of geography.

## Data source: The database reports bilateral trade, including international and
##              intra-national trade, at the aggregated manufacturing level for 
##              69 countries for the period 1986-2006, provided by Thomas Zylkin,
##              based on UN COMTRADE, CEPII TradeProd and UN UNIDO INDSTAT 
##              databases. Standard gravity variables such as distance,  
##              continuous borders, and common language, are taken from the CEPII 
##              GeoDist database.


##***************************** PRELIMINARY STEP ********************************


library(haven)                          # for read_dta
library(sandwich)                       # for vcocHC
library(lmtest)                         # for coeftest

## Set directory path, where "$input" refers to the path of the main folder 
## "Practical Guide to Gravity"	
advGuide.path <- "D:/Reiter/Projects/WTO_Guides/Advanced Guide to TPA/"
setwd(paste0(advGuide.path, "Chapter2"))



##************************ OPEN AND MANAGE THE DATABASE *************************
## Open the database according to the Stata version you are using
twb <- read_dta("Datasets/Chapter2Application1.dta")

## I prefer to use data.tables:
setDT(twb)

## Consider data only for the year == 2006
twb <- twb[year == 2006, ]

## Create the log of distance variable
twb[, ln_DIST := log(DIST)]

## Create the international border dummy variable
twb[, INTL_BRDR := 1]
twb[exporter == importer, INTL_BRDR := 0]

## Create aggregate output by exporter
twb[, Y := sum(trade), by = .(exporter)]

## Create aggregate expenditure by importer
twb[, E := sum(trade), by = .(importer)]

## Chose a country for reference group: GERMANY
## The country code of the reference country is set to "ZZZ" so that the exporter
## and exporter fixed effects of the reference country are always the last ones
## created
twb[exporter == "DEU", exporter := "AAA"]
twb[importer == "DEU", importer := "AAA"]

tmp <- twb[importer == "AAA", unique(E)]
twb[, E_R := tmp]

## Create exporter fixed effects -> do in regression call
## Create importer fixed effects -> do in regression call

## Set the number of exporter fixed effects variables
N <- twb[, uniqueN(exporter)]
N_1 <- N - 1

## Save data
save(twb, file = "Datasets/1_TradeWithoutBorder.RData")

##************************ GENERAL EQUILIBRIUM ANALYSIS *************************

## Step I: Solve the baseline gravity model

## Step I.a. Obtain estimates of trade costs and trade elasticities baseline 
##			indexes
## Estimate the gravity model in the "baseline" scenario with the PPML estimator:
res.Ia <- glm(trade ~ exporter + importer + ln_DIST + CNTG + INTL_BRDR - 1,
              family = "quasipoisson", data = twb)
summary(res.Ia)


## adjust the standard errors to be heteroskedastic
res.Ia.1 <- coeftest(res.Ia, vcov = vcovHC(res.Ia, type = "HC1"))
res.Ia.1

twb[, tradehat_BLN := res.Ia$fitted]

## Step I.b. Construct baseline indexes	
## Based on the estimated exporter and importer fixed effects, create
## the actual set of fixed effects

## extract the coefficients
fixed.effects <- as.data.table(tidy(res.Ia))
fixed.effects.importer <- fixed.effects[term %like% "importer", .(term, estimate)]
fixed.effects.exporter <- fixed.effects[term %like% "exporter", .(term, estimate)]

## importer fixed effects
setnames(fixed.effects.importer, "term", "importer")
fixed.effects.importer[, exp_chi_BLN := exp(estimate)]
fixed.effects.importer[, importer := substring(importer, 9, 100)]
fixed.effects.importer[, estimate := NULL]

## add estimate for DEU back to the table
fixed.effects.importer <- rbind(fixed.effects.importer,
                                list("AAA", 1))

## exporter fixed effects
setnames(fixed.effects.exporter, "term", "exporter")
fixed.effects.exporter[, exp_pi_BLN := exp(estimate)]
fixed.effects.exporter[, exporter := substring(exporter, 9, 100)]
fixed.effects.exporter[, estimate := NULL]

## add both fixed effects to the data.table
setkey(fixed.effects.exporter, exporter)
setkey(twb, exporter)
twb <- twb[fixed.effects.exporter]

setkey(fixed.effects.importer, importer)
setkey(twb, importer)
twb <- twb[fixed.effects.importer]

## Compute the variable of bilateral trade costs
twb[, tij_BLN := exp(res.Ia$coefficients["ln_DIST"] * ln_DIST + res.Ia$coefficients["CNTG"] * CNTG + res.Ia$coefficients["INTL_BRDR"] * INTL_BRDR)]


## Compute the outward and inward multilateral resistances using the 
## additive property of the PPML estimator that links the exporter and  
## importer fixed effects with their respective multilateral resistances
## taking into account the normalisation imposed
twb[, OMR_BLN := Y * E_R / exp_pi_BLN]
twb[, IMR_BLN := E / (exp_chi_BLN * E_R)]

## Compute the estimated level of international trade in the baseline for
## the given level of ouptput and expenditures
twb[exporter != importer, Xi_BLN := sum(tradehat_BLN), by = .(exporter)]
twb[, Y_BLN := Y]
twb[, E_BLN := E]


## Step II: Define a conterfactual scenario
## The counterfactual scenario consists in removing the international borders
## by constraining the parameter associated with the variable INTL_BRDR to be
## zero and assuming the effects of geographic variables (DIST and CNTG) 
## remain the same.

## Constructing the counterfactual bilateral trade costs	by imposing the
## constraints associated with the counterfactual scenario
## Option 1:
generate INTL_BRDR_CFL = 0
generate tij_CFL = exp(_b[ln_DIST]*ln_DIST + _b[CNTG]*CNTG + _b[INTL_BRDR]*INTL_BRDR_CFL)

## Option 2:
##	generate tij_CFL = exp(_b[ln_DIST]*ln_DIST + _b[CNTG]*CNTG + _b[INTL_BRDR]*INTL_BRDR*0)


## Step III: Solve the counterfactual model

## Step III.a.: Obtain conditional general equilibrium effects

## (i):	Estimate the gravity model by imposing the constraints associated 
## 		with the counterfactual scenario. The constraint is defined  
## 		separately by taking the log of the counterfactual bilateral trade 
## 		costs. The parameter of thisexpression will be constrainted to be 
##		equal to 1 in the ppml estimator	

## Specify the constraint in log
generate ln_tij_CFL = log(tij_CFL)	

## Re-create the exporters and imports fixed effects
drop EXPORTER_FE* IMPORTER_FE*
quietly tabulate exporter, generate(EXPORTER_FE)
quietly tabulate importer, generate(IMPORTER_FE)

## Estimate the constrained gravity model and generate predicted trade
## value
ppml trade EXPORTER_FE* IMPORTER_FE1-IMPORTER_FE$N_1 , iter(30) noconst offset(ln_tij_CFL)
predict tradehat_CDL, mu

## (ii):	Construct conditional general equilibrium multilateral resistances

## Based on the estimated exporter and importer fixed effects, create
## the actual set of counterfactual fixed effects	
forvalues i = 1(1)$N_1 {
   quietly replace EXPORTER_FE`i' = EXPORTER_FE`i' * exp(_b[EXPORTER_FE`i'])
   quietly replace IMPORTER_FE`i' = IMPORTER_FE`i' * exp(_b[IMPORTER_FE`i'])
   }

## Create the exporter and importer fixed effects for the country of 
## reference (Germany)
quietly replace EXPORTER_FE$N = EXPORTER_FE$N * exp(_b[EXPORTER_FE$N ])
quietly replace IMPORTER_FE$N = IMPORTER_FE$N * exp(0)

## Create the variables stacking all the non-zero exporter and importer 
## fixed effects, respectively		
egen exp_pi_CDL = rowtotal( EXPORTER_FE1-EXPORTER_FE$N )
egen exp_chi_CDL = rowtotal( IMPORTER_FE1-IMPORTER_FE$N )

## Compute the outward and inward multilateral resistances 				
generate OMR_CDL = Y * E_R / exp_pi_CDL
generate IMR_CDL = E / (exp_chi_CDL * E_R)

## Compute the estimated level of conditional general equilibrium 
## international trade for the given level of ouptput and expenditures		
generate tempXi_CDL = tradehat_CDL if exporter != importer
bysort exporter: egen Xi_CDL = sum(tempXi_CDL)
drop tempXi_CDL


## Step III.b: Obtain full endowment general equilibrium effects

## Create the iterative procedure by specifying the initial variables, 
## where s = 0 stands for the baseline (BLN) value and s = 1 stands for  
## the conditional general equilibrium (CD) value

## The constant elasticity of substitutin is taken from the literature
scalar sigma = 7

## The parameter phi links the value of output with expenditures
generate  phi = E/Y if exporter == importer

## Compute the change in bilateral trade costs resulting from the 
## counterfactual
generate change_tij = tij_CFL / tij_BLN	

## Re-specify the variables in the baseline and conditional scenarios
## Output 
generate Y_0 = Y
generate Y_1 = Y

## Expenditures, including with respect to the reference country   
generate E_0 = E
generate E_R_0 = E_R
generate E_1 = E
generate E_R_1 = E_R			

## Predicted level of trade 
generate tradehat_1 = tradehat_CDL


## (i)	Allow for endogenous factory-gate prices

## Re-specify the factory-gate prices under the baseline and 
## conditional scenarios				
generate exp_pi_0 = exp_pi_BLN
generate tempexp_pi_ii_0 = exp_pi_0 if exporter == importer
bysort importer: egen exp_pi_j_0 = mean(tempexp_pi_ii_0)
generate exp_pi_1 = exp_pi_CDL
generate tempexp_pi_ii_1 = exp_pi_1 if exporter == importer
bysort importer: egen exp_pi_j_1 = mean(tempexp_pi_ii_1)
drop tempexp_pi_ii_*
generate exp_chi_0 = exp_chi_BLN	
generate exp_chi_1 = exp_chi_CDL	

## Compute the first order change in factory-gate prices	in the 
## baseline and conditional scenarios
generate change_pricei_0 = 0				
generate change_pricei_1 = ((exp_pi_1 / exp_pi_0) / (E_R_1 / E_R_0))^(1/(1-sigma))
generate change_pricej_1 = ((exp_pi_j_1 / exp_pi_j_0) / (E_R_1 / E_R_0))^(1/(1-sigma))

## Re-specify the outward and inward multilateral resistances in the
## baseline and conditional scenarios
generate OMR_FULL_0 = Y_0 * E_R_0 / exp_pi_0
generate IMR_FULL_0 = E_0 / (exp_chi_0 * E_R_0)		
generate IMR_FULL_1 = E_1 / (exp_chi_1 * E_R_1)
generate OMR_FULL_1 = Y_1 * E_R_1 / exp_pi_1

## Compute initial change in outward and multilateral resitances, which 
## are set to zero		
generate change_IMR_FULL_1 = exp(0)		
generate change_OMR_FULL_1 = exp(0)
