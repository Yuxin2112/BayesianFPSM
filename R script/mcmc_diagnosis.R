#########################################################################
## This R scrip contains all necessary code for 
## creating trace plot & density plots, as well as Geweke test

## created 2023-3-31 by Yuxin Huang 
## updated 2023-5-24 by Yuxin Huang
##########################################################################
setwd("C:/Users/n11117761/Work/Paper2")
## load package

library(coda)
library(mcmcplots)

library(bayesplot)
library(ggplot2)
#library(gridExtra)

## load .rds file save from previous analysis

WinBUGS_output <- readRDS("mcmc1.rds")

## convert Large bugs object to mcmc object for diagnosis test

MCMC1 <- as.mcmc.bugs(WinBUGS_output)

## save the trace plot & density plot

png(filename = "mcmc_beta.png",
    width = 520, height = 380)
plot(MCMC1[,c("beta[1]","beta[2]")],
     trace = TRUE, density = TRUE, smooth = TRUE, auto.layout = TRUE, 
     type = "l",
     col=c(rgb(red = 0, green = 0, blue = 1, alpha = 0.5),rgb(red = 1, green = 0.2, blue = 0.4, alpha = 0.5)),
     xlab = "Iterations", ylab = "value",
)
dev.off()

png(filename = "mcmc_gamma_1.png",
    width = 520, height = 520)
plot(MCMC1[,c("gamma[1]","gamma[2]","gamma[3]")],
     trace = TRUE, density = TRUE, smooth = TRUE, auto.layout = TRUE, 
     type = "l",
     col=c(rgb(red = 0, green = 0, blue = 1, alpha = 0.5),rgb(red = 1, green = 0.2, blue = 0.4, alpha = 0.5)),
     xlab = "Iterations", ylab = "value")
dev.off()


png(filename = "mcmc_gamma_2.png",
    width = 520, height = 520)
plot(MCMC1[,c("gamma[1]","gamma[2]","gamma[3]")],
     trace = TRUE, density = TRUE, smooth = TRUE, auto.layout = TRUE, 
     type = "l",
     col=c(rgb(red = 0, green = 0, blue = 1, alpha = 0.5),rgb(red = 1, green = 0.2, blue = 0.4, alpha = 0.5)),
     xlab = "Iterations", ylab = "value")
dev.off()

# Diagnosis plot for spatial random effects
# I am not sure which grids I should look at, therefore I just randomly pick 3 grids

png(filename = "mcmc_spatial.png",
    width = 520, height = 520)
plot(MCMC1[,c("ranuv[58]","ranuv[200]","ranuv[400]")],
     trace = TRUE, density = TRUE, smooth = TRUE, auto.layout = TRUE, 
     type = "l",
     col=c(rgb(red = 0, green = 0, blue = 1, alpha = 0.5),rgb(red = 1, green = 0.2, blue = 0.4, alpha = 0.5)),
     xlab = "Iterations", ylab = "value")
dev.off()


#it can also be done by traceplot() and densplot()

################################################################################

## Convergence test
# As we only do one chain for now, therefore I conduct Geweke test only

## Geweke test
geweke_results <- geweke.diag(MCMC1)
## convert as datafame
geweke_results_df <- data.frame(geweke_results[[1]])

## Gelman test
#gelman.diag(post_Gamma5[,"ranuv[58]"],
#            confidence = 0.95)

## autocorrelation plot (if needed)
#autocorr.plot(post_Gamma5[,"gamma[1]"],lag.max = 4000,
#              main="gamma1")

## save the results of test as .csv file

write.csv(geweke_results_df,
          file = "geweke_results.csv",row.names = TRUE)

################################################################

## compute 2.5%, 50% and 97.5% of posterior distribution for each parameter

## save the posterior distribution to a dataframe

pars_mcmc <- as.data.frame(WinBUGS_output$sims.list)

## create another dataframe for quantile of parameters

pars_name <- colnames(pars_mcmc)
npars <- ncol(pars_mcmc)
pars_quantile <- structure(list(para.Name = pars_name, 
                                `2.5` = rep(0,npars), `50` = rep(0,npars), `97.5` = rep(0,npars)), 
                           .Names = c("variable", "q2", "q50", "q97"), 
                           class = "data.frame", row.names = pars_name)


for (i in 1:npars) {
  pars_quantile[i,2] = quantile(pars_mcmc[,i],.025)
  pars_quantile[i,3] = quantile(pars_mcmc[,i],.50)
  pars_quantile[i,4] = quantile(pars_mcmc[,i],.975)
}

# consistent with Stata .do file, switch the position of beta and gamma (ignore this step for now)
#tmp_row_beta <- para_quantile[c(7,8),]
#tep_row <- para_quantile[-c(7,8),]
#pars_quantile <- rbind(tmp_row_beta,tep_row)


# compute v
#for (i in 1:478) {
#  para_quantile[(i+486),2] <- log(para_quantile[i+486,2]) - para_quantile[(i+8),2]
#  para_quantile[(i+486),3] <- log(para_quantile[i+486,3]) - para_quantile[(i+8),3]
#  para_quantile[(i+486),4] <- log(para_quantile[i+486,4]) - para_quantile[(i+8),4]
#}

## save the final result in .dta for Stata
#write.dta(para_quantile,"C:/Users/n11117761/Work/Data/WinBUGS_output/quantiles_est_gamma.dta")

## save the final result in a .csv file
write.csv(pars_quantile,
          file = "pars_quantile.csv",row.names = FALSE)

