# Perform validation of predictions of precipitation occurrence

library(h5)
library(verification)

ds <- h5file('/storage/home/jwo118/scratch/ncar_statmos/data/ncar_statmos_cheswx_20110101_20151231.nc')

pred <- as.vector(ds['prcp'][] > 0)
obs <- as.vector(ds['prcp_obs'][] > 0)

mask_fin <- (is.finite(pred) & is.finite(obs))
pred <- pred[mask_fin]
obs <- obs[mask_fin]

tstats <- table.stats(obs,pred)

# Confusion matrix with predicted as columns
print(tstats$tab)

# Percent correct
print(sprintf("Percent Correct: %.2f",tstats$PC))

# Specificity: true negative rate
print(sprintf("Specificity (true negative rate): %.2f",tstats$tab[1,1]/sum(tstats$tab[1,])))

# Sensitivity: true positive rate
print(sprintf("Sensitivity (true positive rate): %.2f",tstats$tab[2,2]/sum(tstats$tab[2,])))

# Bias
print(sprintf("Bias (1 = unbiased): %.2f",tstats$BIAS))

#Heidke Skill Score
print(sprintf("Heidke Skill Score: %.2f",tstats$HSS))
