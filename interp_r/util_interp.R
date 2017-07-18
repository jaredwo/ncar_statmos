# Utility functions for interpolating precipitation occurence and amount.
###############################################################################

library(caret)
library(randomForest)
library(mgcv)
library(verification)
library(quantregForest)
library(gstat)
library(pROC)

.ecdf_linear <- function(x) {
	
	x <- sort(x)
	n <- length(x)
	q_step <- 1/n #1/(n-1)
	probs <- seq(q_step,1,q_step)#seq(0,1,q_step)
	
	return(approxfun(x, probs, method = "linear", rule=2, yleft=0, yright=1,ties=max))
	
}

.quantile_linear <- function(x) {
	
	x <- sort(x)
	return(approxfun(.ecdf_linear(x)(x), x, method = "linear", rule=2, ties=max))
	
}

optimal_threshold <- function(o,p) {
	
	# Determine global optimal wet probability threshold
	tseq <- seq(0,1,.01)
	bias <- sapply(tseq,function(a_t){table.stats(o,p > a_t)$BIAS})
	othres <- tseq[which.min(abs(1-bias))]
	
	return(othres)
	
}

# K-fold cross validation for precipitation occurrence
kfold_xval <- function(a_obs, vname, stns_xval, folds=10, pred_func, a_seed=NULL, ...)
{
	# Subset xval stations to only those that have observations
	stns_xval <- stns_xval[stns_xval$station_id %in% a_obs$station_id,]
	stns_xval@data[,vname] <- a_obs[stns_xval$station_id,]@data[,vname]
	
	n <- nrow(stns_xval)
	
	predicts <- rep(NA, n)
	names(predicts) <- stns_xval$station_id
	
	if (!is.null(a_seed))
	{
		set.seed(a_seed)
	}
	
	rand <- createFolds(stns_xval@data[,vname], k=folds, list=FALSE)
	foldnum <- sort(unique(rand))
	
	for (i in foldnum) 
	{
		
		rows.in <- rand != i
		rows.out <- rand == i
		
		in_stns_xval <- stns_xval[rows.in,]
		out_stns_xval <- stns_xval[rows.out,]
		
		in_obs <- a_obs[! a_obs$station_id %in% out_stns_xval$station_id,]
		
		a_predicts <- tryCatch(pred_func(in_obs, newdata=out_stns_xval, ...), error=function(e) e)
		if (is(a_predicts,"simpleError")) {
			
			if (a_predicts$message == "Need at least two classes to do classification.")
				a_predicts <- rep(unique(in_obs[[vname]]),nrow(out_stns_xval))
			else
				stop(a_predicts)
			
		}
		
		predicts[rows.out] <- a_predicts      
		
	}
	
	stns_xval@data[,paste(vname,"p",sep='_')] <- predicts
	
	# reset seed if set
	if (!is.null(a_seed))
	{
		set.seed(NULL)
	}
	
	return(stns_xval)	
}

func_rf_wet <- function(obs_pts, a_model=NULL, newdata=NULL, ntree=2000, a_formula=as.formula('wetf~elevation+longitude1+latitude1')) {
	
	if (is.null(a_model)) {
		
		obs_pts$wetf <- factor(obs_pts$wet, levels=c(TRUE,FALSE),ordered=TRUE)
		a_model <- randomForest(a_formula, data=obs_pts,ntree=ntree)
		
	}
	
	same_data <- tryCatch(all(coordinates(obs_pts) == coordinates(newdata)),
			error=function(e) FALSE)
	
	if (is.null(newdata)) {
		
		return(a_model)
		
	} else if (class(newdata)=='RasterStack') {
		
		predicts <- predict(newdata,a_model,type='prob')# > attr(a_model,'wett')  
		return(predicts)
		
	} else if (same_data) {
		
		predicts <- predict(a_model,type='prob')[,1]# > attr(a_model,'wett') 
		
	} else {
		
		predicts <- predict(a_model,newdata = newdata,type='prob')[,1]# > attr(a_model,'wett') 
		#predicts <- predicts > attr(a_model,'wett')
		return(predicts)
	}
	
}

func_gam_wet <- function(obs_pts, a_model=NULL, newdata=NULL,
		a_formula=as.formula('wet~s(elevation,longitude1,latitude1,k=1)+s(longitude1,latitude1)')) {
	
	if (is.null(a_model)) {
		
		options(warn=-1)
		a_model <- gam(a_formula, data=obs_pts, family=binomial)
		options(warn=0)
		
		wett <- optimal_threshold(obs_pts$wet, predict(a_model,type='response'))
		attr(a_model,'wett') <- wett
		
	}
	
	if (is.null(newdata)) {
		
		return(a_model)
		
	} else if (class(newdata)=='RasterStack') {
		
		predicts <- predict(newdata, a_model, type='response')
		#predicts <- predicts > attr(a_model,'wett')    
		return(predicts)
		
	} else {
		
		predicts <- as.numeric(predict(a_model, newdata = newdata, type='response'))
		#predicts <- predicts > attr(a_model,'wett')
		return(predicts)
	}
	
	
}

func_rfq_amt <- function(a_obs_wet, a_pts_wet, ntree=2000) {
	
	df_obs_wet <- as.data.frame(a_obs_wet)[,c('elevation','longitude1','latitude1')]
	
	n <- nrow(df_obs_wet)
	
	cond_ecdf <- list()
	
	rand <- createFolds(a_obs_wet@data[,'prcp'], k=10, list=FALSE)
	foldnum <- sort(unique(rand))
	
	for (i in foldnum) 
	{
		
		rows.in <- rand != i
		rows.out <- rand == i
		
		in_stns <- df_obs_wet[rows.in,]
		in_prcp <- a_obs_wet$prcp[rows.in]
		out_stns <- df_obs_wet[rows.out,]
		
		a_qrf_model <- quantregForest(in_stns, in_prcp, ntree=ntree)
		a_cond_ecdf <- predict(a_qrf_model, out_stns, what=.ecdf_linear)
		cond_ecdf[which(rows.out)] <- a_cond_ecdf    
		
	}
	
	a_obs_wet$prcp_q <- sapply(seq(nrow(a_obs_wet)),function(i) cond_ecdf[[i]](a_obs_wet$prcp[i]))
	
	###########################################################################
	evgm <- variogram(prcp~1,a_obs_wet,cutoff=1000,width=15)
	if (sum(evgm$np)==1) stop("Error: Not enough obs to build vgm")
	
	options(warn=-1)
	vgm_global <- fit.variogram(evgm, vgm(c('Gau','Sph','Exp','Mat')), fit.kappa=TRUE)
	options(warn=0)
	sill <- sum(vgm_global[,'psill'])
	vgm_global[1,'psill'] <- (vgm_global[1,'psill']/sill)*var(a_obs_wet$prcp_q)
	vgm_global[2,'psill'] <- (vgm_global[2,'psill']/sill)*var(a_obs_wet$prcp_q)
	
	nngh_krig <- 10
	krig_prcp_q <- krige(prcp_q~1,a_obs_wet,a_pts_wet, vgm_global,nmax=nngh_krig,
			             nsim=1, debug.level=0)
	krig_prcp_q$sim1[krig_prcp_q$sim1 < 0] = 0			 
	krig_prcp_q$sim1[krig_prcp_q$sim1 > 1] = 1			 
#	krig_prcp_q$sim1[] <- inv.logit(krig_prcp_q$sim1)
	
	##############################################################################
	
	# Check that all points fall within data range of neighboring points
	# Get nngh_krig nearest stns to each grid point
	gpts <- krig_prcp_q
	gpts_m <- gpts@coords
	opts_m <- a_obs_wet@coords
	if (nngh_krig > nrow(a_obs_wet)) nngh_wet <- nrow(a_obs_wet) else nngh_wet <- nngh_krig  
	knn_stns <- apply(gpts_m,1,function(x) order(spDistsN1(opts_m,x,longlat=TRUE))[1:nngh_wet] )
	knn_stns <- t(knn_stns)
	gpts$zmin <- apply(knn_stns,1,function(x)  min(a_obs_wet$prcp_q[x]))
	gpts$zmax <- apply(knn_stns,1,function(x)  max(a_obs_wet$prcp_q[x]))
	
	mask_too_low <- gpts$sim1 < gpts$zmin
	mask_too_high <- gpts$sim1 > gpts$zmax
	gpts$sim1[mask_too_low] <- gpts$zmin[mask_too_low]
	gpts$sim1[mask_too_high] <- gpts$zmax[mask_too_high]
	
	#gpts$sim1 <- inv.logit(gpts$sim1)
	a_pts_wet$prcp_q <- gpts$sim1
	##############################################################################                                  
		
	qrf_model <- quantregForest(df_obs_wet,a_obs_wet$prcp,keep.inbag=TRUE,ntree=ntree)
	
	cond_q_gpts <- predict(qrf_model,as.data.frame(a_pts_wet)[,c('elevation','longitude1','latitude1')],
			               what=.quantile_linear)
	gpts_pred_prcp <- sapply(seq(length(cond_q_gpts)),function(i) cond_q_gpts[[i]](a_pts_wet$prcp_q[i]))                              
		
	return(gpts_pred_prcp)
	
}