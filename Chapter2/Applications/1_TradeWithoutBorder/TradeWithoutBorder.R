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

load("Datasets/1_TradeWithoutBorder.RData")

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
twb[, Xi_BLN := tradehat_BLN]
twb[exporter == importer, Xi_BLN := 0]
twb[, Xi_BLN := sum(Xi_BLN),
    by = .(exporter)]

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
twb[, INTL_BRDR_CFL := 0]
twb[, tij_CFL := exp(res.Ia$coefficients["ln_DIST"] * ln_DIST + res.Ia$coefficients["CNTG"] * CNTG + res.Ia$coefficients["INTL_BRDR"] * INTL_BRDR_CFL)]

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
twb[, ln_tij_CFL := log(tij_CFL)]

## Re-create the exporters and imports fixed effects
## not needed in R

## Estimate the constrained gravity model and generate predicted trade
## value
res.IIIa <- glm(trade ~ exporter + importer + offset(ln_tij_CFL) - 1,
              family = "quasipoisson", data = twb)
summary(res.IIIa)

## adjust the standard errors to be heteroskedastic
res.IIIa.1 <- coeftest(res.IIIa, vcov = vcovHC(res.Ia, type = "HC1"))
res.IIIa.1

twb[, tradehat_CDL := res.IIIa$fitted]

## (ii):	Construct conditional general equilibrium multilateral resistances

## Based on the estimated exporter and importer fixed effects, create
## the actual set of counterfactual fixed effects	

## extract the coefficients
fixed.effects <- as.data.table(tidy(res.IIIa))
fixed.effects.importer <- fixed.effects[term %like% "importer", .(term, estimate)]
fixed.effects.exporter <- fixed.effects[term %like% "exporter", .(term, estimate)]

## importer fixed effects
setnames(fixed.effects.importer, "term", "importer")
fixed.effects.importer[, exp_chi_CDL := exp(estimate)]
fixed.effects.importer[, importer := substring(importer, 9, 100)]
fixed.effects.importer[, estimate := NULL]

## add estimate for DEU back to the table
fixed.effects.importer <- rbind(fixed.effects.importer,
                                list("AAA", 1))

## exporter fixed effects
setnames(fixed.effects.exporter, "term", "exporter")
fixed.effects.exporter[, exp_pi_CDL := exp(estimate)]
fixed.effects.exporter[, exporter := substring(exporter, 9, 100)]
fixed.effects.exporter[, estimate := NULL]

## add both fixed effects to the data.table
setkey(fixed.effects.exporter, exporter)
setkey(twb, exporter)
twb <- twb[fixed.effects.exporter]

setkey(fixed.effects.importer, importer)
setkey(twb, importer)
twb <- twb[fixed.effects.importer]


## Compute the outward and inward multilateral resistances
twb[, OMR_CDL := Y * E_R / exp_pi_CDL]
twb[, IMR_CDL := E / (exp_chi_CDL * E_R)]

## Compute the estimated level of conditional general equilibrium 
## international trade for the given level of ouptput and expenditures
twb[, Xi_CDL := tradehat_CDL]
twb[exporter == importer, Xi_CDL := 0]
twb[, Xi_CDL := sum(Xi_CDL),
    by = .(exporter)]



## 
## Step III.b: Obtain full endowment general equilibrium effects
## 

## Create the iterative procedure by specifying the initial variables, 
## where s = 0 stands for the baseline (BLN) value and s = 1 stands for  
## the conditional general equilibrium (CD) value

## The constant elasticity of substitutin is taken from the literature
sigma <- 7

## The parameter phi links the value of output with expenditures
twb[exporter == importer, phi := E/Y]

## Compute the change in bilateral trade costs resulting from the 
## counterfactual
twb[, change_tij := tij_CFL / tij_BLN]



## For the iterative procedure
## "s" will be the current iteration
## "sm1" will be the previous (read "s minus 1") iteration

## Re-specify the variables in the baseline and conditional scenarios
## Output
twb[, Y_sm1 := Y]
## twb[, Y_s := Y]

## Expenditures, including with respect to the reference country
twb[, E_sm1 := E]
twb[, E_R_sm1 := E_R]
## twb[, E_s := E]
## twb[, E_R_s := E_R]

## Predicted level of trade
twb[, tradehat_sm1 := tradehat_CDL]
## twb[, tradehat_s := tradehat_CDL]

## twb.copy <- copy(twb)
## twb <- copy(twb.copy)


## (i)	Allow for endogenous factory-gate prices

## Re-specify the factory-gate prices under the baseline and 
## conditional scenarios
exppi.1 <- twb[exporter == importer, .(exp_pi_j_sm1 = exp_pi_BLN), keyby = .(importer)]
exppi.2 <- twb[exporter == importer, .(exp_pi_j_s = exp_pi_CDL), keyby = .(importer)]
setkey(twb, importer)
twb <- twb[exppi.1][exppi.2]

twb[, exp_chi_sm1 := exp_chi_CDL]
## twb[, exp_chi_sm1 := exp_chi_BLN]
## twb[, exp_chi_s := exp_chi_CDL]

## Compute the first order change in factory-gate prices in the 
## baseline and conditional scenarios
## twb[, change_pricei_sm1 := 0]

twb[, exp_pi_sm1 := exp_pi_CDL]
## twb[, exp_pi_sm1 := exp_pi_BLN]
## twb[, exp_pi_s := exp_pi_CDL]
## twb[, change_pricei_s := (exp_pi_CDL / exp_pi_BLN)^(1/(1 - sigma))]
twb[, change_pricei_sm1 := (exp_pi_CDL / exp_pi_BLN)^(1/(1 - sigma))]
twb[, change_pricej_sm1 := (exp_pi_j_s / exp_pi_j_sm1)^(1/(1 - sigma))]

## Re-specify the outward and inward multilateral resistances in the
## baseline and conditional scenarios
twb[, OMR_FULL_sm1 := Y_sm1 * E_R_sm1 / exp_pi_sm1]
twb[, IMR_FULL_sm1 := E_sm1 / (exp_chi_sm1 * E_R_sm1)]
## twb[, IMR_FULL_s := E_sm1 / (exp_chi_sm1 * E_R_sm1)]
## twb[, OMR_FULL_s := Y_sm1 * E_R_sm1 / exp_pi_sm1]

## Compute initial change in outward and multilateral resitances, which 
## are set to zero		
twb[, change_IMR_FULL_sm1 := exp(0)]
twb[, change_OMR_FULL_sm1 := exp(0)]

twb[, exp_pi_j_sm1 := NULL]
setnames(twb, "exp_pi_j_s", "exp_pi_j_sm1")

##***************************************************************************
##******************* Start of the Iterative Procedure  *********************

## Set the criteria of convergence, namely that either the standard errors or
## maximum of the difference between two iterations of the factory-gate 
## prices are smaller than 0.01, where s is the number of iterations	
sd_dif_change_pi <- 1
max_dif_change_pi <- 1
iter <- 1L


while((sd_dif_change_pi > 0.001) | (max_dif_change_pi > 0.001)) {

    ## (ii)	Allow for endogenous income, expenditures and trade	
    twb[, trade_s := tradehat_sm1 * change_pricei_sm1 * change_pricej_sm1 / (change_OMR_FULL_sm1 * change_IMR_FULL_sm1)]

    ## (iii) Estimation of the structural gravity model
    ppml.tmp <- glm(trade_s ~ exporter + importer + offset(ln_tij_CFL) - 1,
                    family = "quasipoisson", data = twb)
    twb[, tradehat_s := ppml.tmp$fitted]

    ## Update output & expenditure
    twb[, Y_s := sum(tradehat_s), by = .(exporter)]

    tmp.Es <- twb[exporter == importer, .(E_s = phi * Y_s), keyby = .(importer)]
    setkey(twb, importer)
    twb <- twb[tmp.Es]

    tmp <- twb[importer == "AAA" & !is.na(E_s), unique(E_s)]
    twb[, E_R_s := tmp]

    ## Update factory-gate prices
    ## exporter fixed effects
    fixed.effects <- as.data.table(tidy(ppml.tmp))
    fixed.effects.exporter <- fixed.effects[term %like% "exporter", .(term, estimate)]
    fixed.effects.importer <- fixed.effects[term %like% "importer", .(term, estimate)]

    ## importer fixed effects
    setnames(fixed.effects.importer, "term", "importer")
    fixed.effects.importer[, exp_chi_s := exp(estimate)]
    fixed.effects.importer[, importer := substring(importer, 9, 100)]
    fixed.effects.importer[, estimate := NULL]
    fixed.effects.importer <- rbind(fixed.effects.importer,
                                    list("AAA", 1))

    ## exporter fixed effects
    setnames(fixed.effects.exporter, "term", "exporter")
    fixed.effects.exporter[, exp_pi_s := exp(estimate)]
    fixed.effects.exporter[, exporter := substring(exporter, 9, 100)]
    fixed.effects.exporter[, estimate := NULL]

    ## add exporter fixed effects to the data.table
    setkey(fixed.effects.exporter, exporter)
    setkey(twb, exporter)
    twb <- twb[fixed.effects.exporter]

    setkey(fixed.effects.importer, importer)
    setkey(twb, importer)
    twb <- twb[fixed.effects.importer]

    tmp.exppij <- twb[exporter == importer, .(exp_pi_j_s = exp_pi_s), keyby = .(importer)]
    ## twb[, exp_pi_j_sm1 := exp_pi_j_s]
    setkey(twb, importer)
    twb <- twb[tmp.exppij]

    ## Update multilateral resistances
    twb[, change_pricei_s := ((exp_pi_s / exp_pi_sm1) / (E_R_s / E_R_sm1))^(1/(1-sigma))]
    twb[, change_pricej_s := ((exp_pi_j_s / exp_pi_j_sm1) / (E_R_s / E_R_sm1))^(1/(1-sigma))]
    twb[, OMR_FULL_s := (Y_s * E_R_s) / exp_pi_s]
    twb[, change_OMR_FULL_s := OMR_FULL_s / OMR_FULL_sm1]

    twb[, IMR_FULL_s := E_s / (exp_chi_s * E_R_s)]
    twb[, change_IMR_FULL_s := IMR_FULL_s / IMR_FULL_sm1]

    ## Iteration until the change in factory-gate prices converges to zero
    twb[, dif_change_pi_s := change_pricei_s - change_pricei_sm1]

    sd_dif_change_pi <- twb[, sd(dif_change_pi_s)]
    max_dif_change_pi <- twb[, max(abs(dif_change_pi_s))]

    message(" => Iteration: ", iter, "   sd: ", sd_dif_change_pi, "   max: ", max_dif_change_pi)

    ## delete columns from previous iteration
    twb[, `:=`(c("Y_sm1", "E_sm1", "E_R_sm1", "tradehat_sm1", "exp_pi_j_sm1", 
                 "exp_chi_sm1", "exp_pi_sm1", "change_pricej_sm1", "change_pricei_sm1",  
                 "OMR_FULL_sm1", "IMR_FULL_sm1", "change_IMR_FULL_sm1", "change_OMR_FULL_sm1",
                 "dif_change_pi_s"),
               NULL)]

    ## rename all column "s" in "sm1"
    s.columns <- colnames(twb)[grepl("s$", colnames(twb))]
    setnames(twb, s.columns, gsub("s$", "sm1", s.columns))

    iter <- iter + 1L
}


##******************** End of the Iterative Procedure  **********************
##***************************************************************************

## (iv)	Construction of the "full endowment general equilibrium" 
##		effects indexes

## Compute the full endowment general equilibrium of factory-gate price
twb[, change_pricei_FULL := ((exp_pi_sm1 / exp_pi_BLN) / (E_R_sm1 / E_R))^(1/(1-sigma))]

## Compute the full endowment general equilibrium of the value output
twb[, Y_FULL := change_pricei_FULL * Y_BLN]

## Compute the full endowment general equilibrium of the value of 
## aggregate expenditures
tmp.Es <- twb[exporter == importer, .(E_FULL = phi * Y_FULL), keyby = .(importer)]
setkey(twb, importer)
twb <- twb[tmp.Es]

## Compute the full endowment general equilibrium of the outward and 
## inward multilateral resistances
twb[, OMR_FULL := Y_FULL * E_R_sm1 / exp_pi_sm1]
twb[, IMR_FULL := E_sm1 / (exp_chi_sm1 * E_R_sm1)]

## Compute the full endowment general equilibrium of the value of 
## bilateral trade
twb[, X_FULL := (Y_FULL * E_FULL * tij_CFL) /(IMR_FULL * OMR_FULL)]

## Compute the full endowment general equilibrium of the value of 
## total international trade
twb[, Xi_FULL := X_FULL]
twb[exporter == importer, Xi_FULL := 0]
twb[, Xi_FULL := sum(Xi_FULL),
    by = .(exporter)]

## Save the conditional and general equilibrium effects results		
save(twb, file = "Applications/1_TradeWithoutBorder/Results/1_TradeWithoutBorder_FULLGE.RData")


## Step IV: Collect, construct, and report indexes of interest
load("Applications/1_TradeWithoutBorder/Results/1_TradeWithoutBorder_FULLGE.RData")
indexes <- c("OMR_FULL",  "OMR_CDL", "OMR_BLN", "change_pricei_FULL",
             "Xi_BLN", "Xi_CDL", "Xi_FULL", "Y_BLN", "Y_FULL")

twb[exporter == "AAA", exporter := "DEU"]
twb[importer == "AAA", importer := "DEU"]

twb.ex <- twb[, lapply(.SD, mean), keyby = .(exporter), .SD = indexes]
setnames(twb.ex, "exporter", "country")

## Percent change in full endowment general equilibrium of factory-gate prices
twb.ex[, change_price_FULL := (change_pricei_FULL - 1) * 100]

## Percent change in full endowment general equilibirum of outward multilateral resistances
twb.ex[, change_OMR_CDL := (OMR_CDL^(1/(1-sigma)) - OMR_BLN^(1/(1-sigma))) / OMR_BLN^(1/(1-sigma)) * 100]

## Percent change in full endowment general equilibrium of outward multilateral resistances
twb.ex[, change_OMR_FULL := (OMR_FULL^(1/(1-sigma)) - OMR_BLN^(1/(1-sigma))) / OMR_BLN^(1/(1-sigma)) * 100]

## Percent change in conditional general equilibrium of bilateral trade
twb.ex[, change_Xi_CDL := (Xi_CDL - Xi_BLN) / Xi_BLN * 100]

## Percent change in full endowment general equilibrium of bilateral trade		
twb.ex[, change_Xi_FULL := (Xi_FULL - Xi_BLN) / Xi_BLN * 100]



## Construct the percentage changes on import/consumption side
indexes <- c("IMR_FULL", "IMR_CDL", "IMR_BLN")

twb.im <- twb[, lapply(.SD, mean), keyby = .(importer), .SD = indexes]
setnames(twb.im, "importer", "country")

## Conditional general equilibrium of inward multilateral resistances
twb.im[, change_IMR_CDL := (IMR_CDL^(1/(1-sigma)) - IMR_BLN^(1/(1-sigma))) / IMR_BLN^(1/(1-sigma)) * 100]

## Full endowment general equilibrium of inward multilateral resistances
twb.im[, change_IMR_FULL := (IMR_FULL^(1/(1-sigma)) - IMR_BLN^(1/(1-sigma))) / IMR_BLN^(1/(1-sigma)) * 100]

## Merge the general equilibrium results from the production and consumption
## sides
twb.res <- twb.im[twb.ex]

## Full endowment general equilibrium of real GDP
twb.res[, rGDP_BLN := Y_BLN / (IMR_BLN^(1 / (1 - sigma)))]
twb.res[, rGDP_FULL := Y_FULL / (IMR_FULL^(1 / (1 - sigma)))]
twb.res[, change_rGDP_FULL := (rGDP_FULL - rGDP_BLN) / rGDP_BLN * 100]

## Keep indexes of interest	
indexes <- c("country", "change_Xi_CDL", "change_Xi_FULL", "change_price_FULL",
             "change_IMR_FULL", "change_rGDP_FULL", "Y_BLN")
twb.res <- twb.res[, indexes, with = FALSE]


twb.res[, ln_Y := log(Y_BLN)]
twb.res[, change_IMR_FULL := -1 * change_IMR_FULL]

## Create a graphic showing the conditional and full endowment general
## equilibrium on exports

plot1 <- melt(twb.res, id = "country", measure = patterns("^change_Xi"),
              variable.name = "scenario", variable.factor = FALSE)
plot1[scenario == "change_Xi_CDL", scenario := "Conditional general equilibrium"]
plot1[scenario == "change_Xi_FULL", scenario := "Full endowment general equilibrium"]
setkey(plot1, country)
plot1 <- plot1[twb.res[, .(country, ln_Y)]]

p <- ggplot(plot1[country != "HKG", ], aes(x = ln_Y, y = value))
p <- p + geom_point(aes(colour = scenario), size = 3)
p <- p + ylab("Percentage change of exports") + xlab("Log value of output")
p <- p + theme(legend.pos = "bottom")
p

ggsave(filename = "Applications/1_TradeWithoutBorder/Results/scatter_trade_output.png")

## Create a graphic showing the impact on real GDP, factory-gate prices and
## -1* inward multilateral resistances
plot2 <- melt(twb.res, id = "country",
              measure = c("change_price_FULL", "change_IMR_FULL", "change_rGDP_FULL"), 
              variable.factor = FALSE)
plot2[variable == "change_price_FULL", variable := "Factory-gate price"]
plot2[variable == "change_rGDP_FULL", variable := "real GDP"]
plot2[variable == "change_IMR_FULL", variable := "-(inward multilateral resistances)"]
setkey(plot2, country)
plot2 <- plot2[twb.res[, .(country, ln_Y)]]

p <- ggplot(plot2, aes(x = ln_Y, y = value))
p <- p + geom_point(aes(colour = variable), size = 3)
p <- p + ylab("Percentage changes") + xlab("Log value of output")
p <- p + theme(legend.pos = "bottom")
p

ggsave(filename = "Applications/1_TradeWithoutBorder/Results/scatter_rGDP_output.png")
