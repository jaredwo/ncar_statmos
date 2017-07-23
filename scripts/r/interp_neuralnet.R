library(neuralnet)

load('/storage/home/jwo118/scratch/ncar_statmos/data/prcp_stns_cheasapeake_20060101_20151231.rdata')

mask_time <- df_prcp$time< as.POSIXct('2015-01-01')
dftrain <- df_prcp[mask_time,]
dftest <- df_prcp[!mask_time,]

#Training data
binresptr <- ifelse(dftrain$prcp > 0, 1, 0)
lattr <- dftrain$latitude
lontr <- dftrain$longitude
elevtr <- dftrain$elevation
prcptr <- dftrain$prcp
ysintr <- dftrain$yday_sin
ycostr <- dftrain$yday_cos
dlon_nn01tr <- dftrain$dlon_nn01
dlat_nn01tr <- dftrain$dlat_nn01
delev_nn01tr <- dftrain$delev_nn01
prcp_nn01tr <- dftrain$prcp_nn01
dlon_nn02tr <- dftrain$dlon_nn02
dlat_nn02tr <- dftrain$dlat_nn02
delev_nn02tr <- dftrain$delev_nn02
prcp_nn02tr <- dftrain$prcp_nn02

trainset1 <- cbind(binresptr,lattr,lontr,elevtr,ysintr,ycostr,
                   dlon_nn01tr, dlat_nn01tr, delev_nn01tr, prcp_nn01tr,
                   dlon_nn02tr, dlat_nn02tr, delev_nn02tr, prcp_nn02tr)

# Testing data

binrespte <- ifelse(dftest$prcp > 0, 1, 0)
latte <- dftest$latitude
lonte <- dftest$longitude
elevte <- dftest$elevation
prcpte <- dftest$prcp
ysinte <- dftest$yday_sin
ycoste <- dftest$yday_cos
dlon_nn01te <- dftest$dlon_nn01
dlat_nn01te <- dftest$dlat_nn01
delev_nn01te <- dftest$delev_nn01
prcp_nn01te <- dftest$prcp_nn01
dlon_nn02te <- dftest$dlon_nn02
dlat_nn02te <- dftest$dlat_nn02
delev_nn02te <- dftest$delev_nn02
prcp_nn02te <- dftest$prcp_nn02

testset1 <- cbind(binrespte,latte,lonte,elevte,ysinte,ycoste,
                  dlon_nn01te, dlat_nn01te, delev_nn01te, prcp_nn01te,
                  dlon_nn02te, dlat_nn02te, delev_nn02te, prcp_nn02te)

#lifesign = "minimal", 

neonet2 <- neuralnet(binresptr ~ lattr + lontr + elevtr + ysintr + ycostr+dlon_nn01tr+ dlat_nn01tr+ delev_nn01tr+ prcp_nn01tr+
                       dlon_nn02tr+ dlat_nn02tr+ delev_nn02tr+ prcp_nn02tr, trainset1, 
                     hidden = 4,  act.fct = "logistic",
                     stepmax = 1e+05,linear.output = FALSE, threshold = 0.1)
