##*******************************************************************************
##***************  CHAPTER 2 - GENERAL EQUILIBRIUM TRADE POLICY  **************** 
##***************		      ANALYSIS WITH STRUCTURAL GRAVITY 	**************** 
##*******************************************************************************

##************  APPLICATION 2: IMPACT OF REGIONAL TRADE AGREEMENTS  ************* 

## This application applies the methods developed by Anderson et al. (2015) in 
## order to investigates the potential effects of removing the North American
## Free Trade Agreement (NAFTA).

## Data source: The database reports bilateral trade, including international and
##              intra-national trade, at the aggregated manufacturing level for 
##              69 countries for the period 1986-2006, provided by Thomas Zylkin,
##              based on UN COMTRADE, CEPII TradeProd and UN UNIDO INDSTAT 
##              databases. Information on RTAs come from Mario Larch's Regional
##              Trade Agreements Database (http://www.ewf.uni-bayreuth.de/en/
##              research/RTA-data/index.html). Standard gravity variables such 
##              as distance, continuous borders, and common language, are taken
##              from the CEPII GeoDist database.


##****************************** PRELIMINARY STEP *******************************

library(data.table)
library(haven)                          # for read_dta
library(sandwich)                       # for vcocHC (robust standard errors)
library(lmtest)                         # for coeftest
library(broom)                          # for tidying of lm objects
library(ggplot2)                        # for beautiful graphics


## Set directory path, where "$input" refers to the path of the main folder 
## "Practical Guide to Gravity"	
advGuide.path <- "D:/Reiter/Projects/WTO_Guides/Advanced Guide to TPA/"
setwd(paste0(advGuide.path, "Chapter2"))

##************************ OPEN AND MANAGE THE DATABASE *************************
## Open the database according to the Stata version you are using
rtaeff <- read_dta("Datasets/Chapter2Application2.dta")
setDT(rtaeff)

library(ormisc)
rtaeff[, (names(rtaeff)) := lapply(.SD, stripAttr),
       .SDcols = names(rtaeff)]

## Consider panel data with 4 year interval (1986, 1990, ..., 2006)
## rtaeff <- rtaeff[year %in% c(1986 + (0:6) * 4)] # 4 should be 0
rtaeff <- rtaeff[year %in% c(1994, 1998)]

## Create the log of distance variable
rtaeff[, ln_DIST := log(DIST)]

## Create aggregate output
rtaeff[, Y := sum(trade), by = .(exporter, year)]

## Create aggregate expenditure
rtaeff[, E := sum(trade), by = .(importer, year)]

## Chose a country for reference group: GERMANY
## The country code of the reference country is set to "ZZZ" so that the exporter
## and exporter fixed effects of the reference country are always the last ones
## created
rtaeff[exporter == "DEU", exporter := "AAA"]
rtaeff[importer == "DEU", importer := "AAA"]
tmp <- rtaeff[importer == "AAA", unique(E)]
rtaeff[, E_R := tmp]

## Create exporter time fixed effects -> will be done in the regression call
## Create importer time fixed effects -> will be done in the regression call

## Rearrange so that country pairs (e.g. NER-PAN, MWI-MAC, NPL-MWI, PAN-MWI,
## NPL-CMR), which will be dropped due to no trade, are last
rtaeff[, X := sum(trade), by = .(pair_id)]
rtaeff <- rtaeff[X != 0, ]
rtaeff[, X := NULL]

## Set additional exogenous parameters
NT <- length(rtaeff[, table(exporter, year)])
NT_1 <- NT - 1

Nyr <- rtaeff[, uniqueN(year)]
NT_yr <- NT - Nyr

NTij <- rtaeff[, uniqueN(pair_id)]
NTij_1 <- NTij - 1
NTij_8 <- NTij - 8

## Save data
save.image("Datasets/2_RTAsEffects.RData")

##************************ GENERAL EQUILIBRIUM ANALYSIS *************************
## load("Datasets/2_RTAsEffects.RData")

## Step I: Solve the baseline gravity model

## Step I.a. Obtain estimates of trade costs and trade elasticities baseline 
##			indexes

## Implementation of Anderson and Yotov (2016) two-stage procedure to 
## construct the full matrix of trade costs, including when there is no
## trade or zero trade

## Stage 1: Obtain the estimates of pair fixed effects and the effects of RTAs
setorder(rtaeff, pair_id)
rtaeff[, pair_id := as.character(pair_id)]
res.Ia <- glm(trade ~  - 1 + pair_id + exporter:year + importer:year + RTA,
              family = "quasipoisson", data = rtaeff)
summary(res.Ia)

test <- tidy(res.Ia)

RTA_est <- res.Ia[["coefficients"]]["RTA"]

## Construct the trade costs from the pair fixed effects
getFixedEffects <- function(lm.obj, fe.string) {
    require(broom)
    fe <- tidy(lm.obj)
    setDT(fe)
    fe <- fe[term %like% paste0("^", fe.string), .(term, estimate)]
    fe[, term := gsub(fe.string, "", term)]
    ## fe[, exp.estimate := exp(estimate)]
    ## setnames(fe, "exp.estimate", paste0("exp.fe.", fe.string))
    setnames(fe, "estimate", paste0("fe.", fe.string))
    setkey(fe, term)
    fe
}

pair.fe.var <- "pair_id"
pair.fe <- getFixedEffects(res.Ia, pair.fe.var)

setkeyv(rtaeff, pair.fe.var)
rtaeff <- rtaeff[pair.fe]

rtaeff[, tij_bar := exp(fe.pair_id)]
rtaeff[, tij_bln := exp(fe.pair_id) + RTA_est * RTA]

## Stage 2:
## Regress the estimates of pair fixed effects on gravity variables
## and country fixed effects Perform the regression for the baseline
## year
rtaeff <- rtaeff[year == 1994, ]

## Specify the dependent variable as the estimates of pair fixed
## effects
rtaeff[, tij := exp(fe.pair_id)]

## Estimate the standard gravity model 
res.tij <- glm(tij ~ - 1 + exporter + importer + ln_DIST + CNTG + LANG + CLNY,
               family = "quasipoisson", data = rtaeff[exporter != importer, ])



ppml tij EXPORTER_FE* IMPORTER_FE* ln_DIST CNTG LANG CLNY if exporter != importer, cluster(pair_id)
estimates store gravity_est

## Create the predicted values 	
predict tij_noRTA, mu
replace tij_noRTA = 1 if exporter == importer

## Replace the missing trade costs with predictions from the
## standard gravity regression
replace tij_bar = tij_noRTA if tij_bar == . 
replace tij_bln = tij_bar * exp(RTA_est*RTA) if tij_bln == .

## Specify the complete set of bilateral trade costs in log to
## be used as a constraint in the PPML estimation of the 
## structural gravity model
generate ln_tij_bln = log(tij_bln)	

## Set the number of exporter fixed effects variables
quietly ds EXPORTER_FE*
global N = ': word count 'r(varlist)'' 
global N_1 = $N - 1	

## Estimate the gravity model in the "baseline" scenario with the PPML
## estimator constrained with the complete set of bilateral trade costs
ppml trade EXPORTER_FE* IMPORTER_FE1-IMPORTER_FE$N_1 , iter(30) noconst offset(ln_tij_bln)
predict tradehat_BLN, mu

## Step I.b. Construct baseline indexes	
## Based on the estimated exporter and importer fixed effects, create
## the actual set of fixed effects
forvalues i = 1 (1) $N_1 {
   quietly replace EXPORTER_FE'i' = EXPORTER_FE'i' * exp(_b[EXPORTER_FE'i'])
   quietly replace IMPORTER_FE'i' = IMPORTER_FE'i' * exp(_b[IMPORTER_FE'i'])
   }

## Create the exporter and importer fixed effects for the country of 
## reference (Germany)
quietly replace EXPORTER_FE$N = EXPORTER_FE$N * exp(_b[EXPORTER_FE$N ])
quietly replace IMPORTER_FE$N = IMPORTER_FE$N * exp(0)

## Create the variables stacking all the non-zero exporter and importer 
## fixed effects, respectively		
egen exp_pi_BLN = rowtotal(EXPORTER_FE1-EXPORTER_FE$N )
egen exp_chi_BLN = rowtotal(IMPORTER_FE1-IMPORTER_FE$N ) 

## Compute the variable of bilateral trade costs, i.e. the fitted trade
## value by omitting the exporter and importer fixed effects		
generate tij_BLN = tij_bln			

## Compute the outward and inward multilateral resistances using the 
## additive property of the PPML estimator that links the exporter and  
## importer fixed effects with their respective multilateral resistances
## taking into account the normalisation imposed
generate OMR_BLN = Y * E_R / exp_pi_BLN
generate IMR_BLN = E / (exp_chi_BLN * E_R)	

## Compute the estimated level of international trade in the baseline for
## the given level of ouptput and expenditures			
generate tempXi_BLN = tradehat_BLN if exporter != importer
bysort exporter: egen Xi_BLN = sum(tempXi_BLN)
drop tempXi_BLN
generate Y_BLN = Y
generate E_BLN = E


## Step II: Define a conterfactual scenario
## The counterfactual scenario consists in removing the impact of the NAFTA
## by re-specifying the RTA variable with zeros for the country pairs 
## associated with the NAFTA
generate RTA_NO_NAFTA = RTA
replace RTA_NO_NAFTA=0 if (exporter == "CAN" & importer == "USA") | (exporter == "CAN" & importer == "MEX") | (exporter == "MEX" & importer == "USA") | (exporter == "MEX" & importer == "CAN") | (exporter == "USA" & importer == "MEX") | (exporter == "USA" & importer == "CAN")

## Constructing the counterfactual bilateral trade costs	by imposing the
## constraints associated with the counterfactual scenario
generate tij_CFL = tij_bar * exp(RTA_est * RTA_NO_NAFTA) 

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
   quietly replace EXPORTER_FE'i' = EXPORTER_FE'i' * exp(_b[EXPORTER_FE'i'])
   quietly replace IMPORTER_FE'i' = IMPORTER_FE'i' * exp(_b[IMPORTER_FE'i'])
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
bysort year: generate phi = E/Y if exporter == importer

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


##***************************************************************************
##******************* Start of the Iterative Procedure  *********************

## Set the criteria of convergence, namely that either the standard errors or
## maximum of the difference between two iterations of the factory-gate 
## prices are smaller than 0.01, where s is the number of iterations	
local s = 3	
local sd_dif_change_pi = 1
local max_dif_change_pi = 1
while ('sd_dif_change_pi' > 0.01) | ('max_dif_change_pi' > 0.01) {
   local s_1 = 's' - 1
   local s_2 = 's' - 2
   local s_3 = 's' - 3
   
   * (ii)	Allow for endogenous income, expenditures and trade	
   generate trade_'s_1' =  tradehat_'s_2' * change_pricei_'s_2' * change_pricej_'s_2' / (change_OMR_FULL_'s_2'*change_IMR_FULL_'s_2')

   
   * (iii)	Estimation of the structural gravity model
   drop EXPORTER_FE* IMPORTER_FE*
   quietly tabulate exporter, generate (EXPORTER_FE)
   quietly tabulate importer, generate (IMPORTER_FE)
   capture ppml trade_'s_1' EXPORTER_FE* IMPORTER_FE*, offset(ln_tij_CFL) noconst iter(30) 
   predict tradehat_'s_1', mu
   
   * Update output & expenditure			
   bysort exporter: egen Y_'s_1' = total(tradehat_'s_1')
   quietly generate tempE_'s_1' = phi * Y_'s_1' if exporter == importer
   bysort importer: egen E_'s_1' = mean(tempE_'s_1')
   quietly generate tempE_R_'s_1' = E_'s_1' if importer == "ZZZ"
   egen E_R_'s_1' = mean(tempE_R_'s_1')
   
   * Update factory-gate prices 
   forvalues i = 1(1)$N_1 {
      quietly replace EXPORTER_FE'i' = EXPORTER_FE'i' * exp(_b[EXPORTER_FE'i'])
      quietly replace IMPORTER_FE'i' = IMPORTER_FE'i' * exp(_b[IMPORTER_FE'i'])
      }
   quietly replace EXPORTER_FE$N = EXPORTER_FE$N * exp(_b[EXPORTER_FE$N ])
   egen exp_pi_'s_1' = rowtotal(EXPORTER_FE1-EXPORTER_FE$N ) 
   quietly generate tempvar1 = exp_pi_'s_1' if exporter == importer
   bysort importer: egen exp_pi_j_'s_1' = mean(tempvar1) 		
   
   * Update multilateral resistances
   generate change_pricei_'s_1' = ((exp_pi_'s_1' / exp_pi_'s_2') / (E_R_'s_1' / E_R_'s_2'))^(1/(1-sigma))
   generate change_pricej_'s_1' = ((exp_pi_j_'s_1' / exp_pi_j_'s_2') / (E_R_'s_1' / E_R_'s_2'))^(1/(1-sigma))
   generate OMR_FULL_'s_1' = (Y_'s_1' * E_R_'s_1') / exp_pi_'s_1' 
   generate change_OMR_FULL_'s_1' = OMR_FULL_'s_1' / OMR_FULL_'s_2'					
   egen exp_chi_'s_1' = rowtotal(IMPORTER_FE1-IMPORTER_FE$N )	
   generate IMR_FULL_'s_1' = E_'s_1' / (exp_chi_'s_1' * E_R_'s_1')
   generate change_IMR_FULL_'s_1' = IMR_FULL_'s_1' / IMR_FULL_'s_2'
   
   * Iteration until the change in factory-gate prices converges to zero
   generate dif_change_pi_'s_1' = change_pricei_'s_2' - change_pricei_'s_3'
   display "************************* iteration number " 's_2' " *************************"
   summarize dif_change_pi_'s_1', format
   display "**********************************************************************"
   display " "
   local sd_dif_change_pi = r(sd)
   local max_dif_change_pi = abs(r(max))	
   
   local s = 's' + 1
   drop temp* 
	}

##******************** End of the Iterative Procedure  **********************
##***************************************************************************

## (iv)	Construction of the "full endowment general equilibrium" 
##		effects indexes
## Use the result of the latest iteration S
local S = 's' - 2
##	forvalues i = 1 (1) $N_1 {
##   		quietly replace IMPORTER_FE'i' = IMPORTER_FE'i' * exp(_b[IMPORTER_FE'i'])
##   	}		
## Compute the full endowment general equilibrium of factory-gate price
generate change_pricei_FULL = ((exp_pi_'S' / exp_pi_0) / (E_R_'S' / E_R_0))^(1/(1-sigma))		

## Compute the full endowment general equilibrium of the value output
generate Y_FULL = change_pricei_FULL  * Y_BLN

## Compute the full endowment general equilibrium of the value of 
## aggregate expenditures
generate tempE_FULL = phi * Y_FULL if exporter == importer
bysort importer: egen E_FULL = mean(tempE_FULL)
drop tempE_FULL

## Compute the full endowment general equilibrium of the outward and 
## inward multilateral resistances 
generate OMR_FULL = Y_FULL * E_R_'S' / exp_pi_'S'
generate IMR_FULL = E_'S' / (exp_chi_'S' * E_R_'S')	

## Compute the full endowment general equilibrium of the value of 
## bilateral trade 
generate X_FULL = (Y_FULL * E_FULL * tij_CFL) /(IMR_FULL * OMR_FULL)			

## Compute the full endowment general equilibrium of the value of 
## total international trade 
generate tempXi_FULL = X_FULL if exporter != importer
bysort exporter: egen Xi_FULL = sum(tempXi_FULL)
drop tempXi_FULL

## Save the conditional and general equilibrium effects results		
save "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULLGE.dta", replace


## Step IV: Collect, construct, and report indexes of interest
use "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULLGE.dta", clear
collapse(mean) OMR_FULL OMR_CDL OMR_BLN change_pricei_FULL Xi_* Y_BLN Y_FULL, by(exporter)
rename exporter country
replace country = "DEU" if country == "ZZZ"
sort country

## Percent change in full endowment general equilibrium of factory-gate prices
generate change_price_FULL = (1 - change_pricei_FULL) / 1 * 100

## Percent change in full endowment general equilibirum of outward multilateral resistances
generate change_OMR_CDL = (-OMR_CDL^(1/(1-sigma)) + OMR_BLN^(1/(1-sigma))) / OMR_CDL^(1/(1-sigma)) * 100

## Percent change in full endowment general equilibrium of outward multilateral resistances			
generate change_OMR_FULL = (-OMR_FULL^(1/(1-sigma)) + OMR_BLN^(1/(1-sigma))) / OMR_FULL^(1/(1-sigma)) * 100

## Percent change in conditional general equilibrium of bilateral trade
generate change_Xi_CDL = (-Xi_CDL + Xi_BLN) / Xi_CDL * 100	

## Percent change in full endowment general equilibrium of bilateral trade		
generate change_Xi_FULL = (-Xi_FULL + Xi_BLN) / Xi_FULL * 100
save "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULL_PROD.dta", replace


## Construct the percentage changes on import/consumption side
use "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULLGE.dta", clear
collapse(mean) IMR_FULL IMR_CDL IMR_BLN, by(importer)
rename importer country
replace country = "DEU" if country == "ZZZ"
sort country		

## Conditional general equilibrium of inward multilateral resistances
generate change_IMR_CDL = (-IMR_CDL^(1/(1-sigma)) + IMR_BLN^(1/(1-sigma))) / IMR_CDL^(1/(1-sigma)) * 100

## Full endowment general equilibrium of inward multilateral resistances
generate change_IMR_FULL = (-IMR_FULL^(1/(1-sigma)) + IMR_BLN^(1/(1-sigma))) / IMR_FULL^(1/(1-sigma)) * 100
save "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULL_CONS.dta", replace

## Merge the general equilibrium results from the production and consumption
## sides
use "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULL_PROD.dta", clear
joinby country using "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULL_CONS.dta"

## Full endowment general equilibrium of real GDP
generate rGDP_BLN = Y_BLN / (IMR_BLN ^(1 / (1 -sigma)))
generate rGDP_FULL = Y_FULL / (IMR_FULL ^(1 / (1 -sigma)))
generate change_rGDP_FULL = (-rGDP_FULL + rGDP_BLN) / rGDP_FULL * 100

## Keep indexes of interest	
keep country change_Xi_CDL change_Xi_FULL change_price_FULL change_IMR_FULL change_rGDP_FULL Y_BLN
order country change_Xi_CDL change_Xi_FULL change_price_FULL change_IMR_FULL change_rGDP_FULL Y_BLN

## Export the results in Excel
export excel using "Applications/2_RTAsEffects/Results/2_RTAsEffects_FULL.xls", firstrow(variables) replace

##*******************************************************************************
