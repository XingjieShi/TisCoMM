# simulation 2
# Type I effor rate evluation under NULL models
library(mvtnorm)
library(glmnet)
library(Rcpp)
library(CoMM)
sourceCpp("vb_aux.cpp") 
source("functions.r")  
source("multiEstB.R")
sourceCpp("PXem_part.cpp")
#sourceCpp("PXem.cpp")

PartType1Err <- function(h_z, h_c, s, rhoX, rhoE, rhoW, nz_ti, Ti, a0, batch, Rep)
{
	PX <- T
	tic <- Sys.time()
	pvalue_comm <- matrix(0, ncol=Ti, nrow=Rep)
	pvalue_mammotpart <- matrix(0, ncol=Ti, nrow=Rep)
	pvalue_predixcan <- matrix(0, ncol=Ti, nrow=Rep)
	pvalue_TWAS <- matrix(0, ncol=Ti, nrow=Rep)
	alpha_mammotpart <- matrix(0, ncol=Ti, nrow=Rep)
	alpha_comm <- matrix(0, ncol=Ti, nrow=Rep)

 	# mammot, mammot_ss, multixcan, multixcan_ss, utmost.
	for(i in 1:Rep)	{
		m <- 400
		n <- 5000
		p <- 400
		# genotypes in eQTL
		x1 <- x2snp(rmvnorm(m, mean=rep(0, p), sigma = AR(rhoX, p)))
		x1_snp <- x1
		x1mean <- matrix(rep(colMeans(x1), m), ncol=p, byrow = T)
		x1sd <- matrix(rep(apply(x1, 2, sd), m), ncol=p, byrow = T)
		x1 <- (x1 - x1mean)/x1sd/sqrt(m)


		# genotypes in GWAS
		x2 <- x2snp(rmvnorm(n, mean=rep(0, p), sigma = AR(rhoX, p)))
		x2_snp <- x2
		x2mean <- matrix(rep(colMeans(x2), n), ncol=p, byrow = T)
		x2sd <- matrix(rep(apply(x2, 2, sd), n), ncol=p, byrow = T)
		x2 <- (x2 - x2mean)/x2sd/sqrt(n)

 
 		# eQTL coefficient matrix, across Ti tissues.
		b <- rnorm(p)
		W <- matrix(0, nrow = p, ncol=Ti)
		idx <- sample(p, s*p, replace=FALSE)
		W[idx, ] <- rmvnorm(s*p, mean=rep(0, Ti), sigma = AR(rhoW, Ti))
  
		B <-  diag(b) %*% W
		lp <- x1 %*% B

		# error matrix in eQTL.
		if (Ti != 1){
			sd.lp = diag(sqrt(diag(var(lp)) * (1/h_c - 1)))
			Ve = sd.lp %*% AR(rhoE, Ti) %*% sd.lp
			E <- rmvnorm(m, mean = rep(0, Ti), sigma = Ve)
		} else {
			sd.lp = sqrt(diag(var(lp)* (1/h_c - 1)))
			E <- rnorm(m, sd = sd.lp)
		}
		y <- lp + E

		if (h_z == 0) z <- rnorm(n, 0, sqrt(3)) else {
			a <- rep(0, Ti)
			#a <- runif(Ti, -1, 1)
			a[1:nz_ti] <- a0
			lp_z <- x2 %*% B %*% a
			sig2 <- c(var(lp_z)) * ((1-h_z)/h_z)
			z <- lp_z + rnorm(n, 0, sqrt(sig2))
  		}

 		maxIter = 3000
 		y <- resid(lm(y~rep(1, m)))
       	z <- resid(lm(z~rep(1, n)))
    
 		# mammot (individual data)
 		# adaptive weight for mammot and mammot_ss.
		w <- matrix(1, nrow=p, ncol=1)
		if(Ti != 1) w = ADW(x1, y)  
 		for (ti in 1:Ti) {
       		constrFactor <- numeric(Ti)
			constrFactor[ti] <- 1
			fit <- mammot_part_test(x1, y, x2, z, w, constrFactor, TRUE, maxIter=2000)
			pvalue_mammotpart[i, ti] <- pchisq(fit$chisq, 1, lower.tail=F)
		}	

		# comm
  		w1 = matrix(rep(1, m), ncol=1)		
		w2 = matrix(rep(1, n), ncol=1)
		for (ti in 1:Ti)	{
			fmHa = CoMM_covar_pxem(y[, ti], z, x1, x2, w1, w2, constr = 0)
			fmH0 = CoMM_covar_pxem(y[, ti], z, x1, x2, w1, w2, constr = 1)
			loglikHa = max(fmHa$loglik, na.rm=T)
			loglikH0 = max(fmH0$loglik, na.rm=T)
			tstat = 2 * (loglikHa - loglikH0)
 			pvalue_comm[i, ti] <- pchisq(tstat, 1, lower.tail=F)
 			alpha_comm[i, ti] <- fmHa$alpha
 		} 

 		# prediXcan
 		for (ti in 1:Ti)	{
 			fit_elasnet <- cv.glmnet(x1, y[, ti], type.measure="mse", alpha=0.5, family="gaussian")
			elasnet_M <- predict(fit_elasnet, s=fit_elasnet$lambda.min, newx=x2)
			if (var(elasnet_M) == 0) {
				pvalue_predixcan[i, ti] <- 1
			} else {
				bet_lm <- coefficients(summary(lm(z~elasnet_M)))
				pvalue_predixcan[i, ti] <- bet_lm[2,4]	
			}
		}


		# TWAS
		for (ti in 1:Ti)	{
			setwd(paste0("./","twasTemp"))
			path <- paste0("part_hc", h_c*10, "_rhoX", rhoX*10, "_rhoW", rhoW*10, "nz_ti", nz_ti, "_batch-", batch)
			base <-  paste(path, ti, i, sep="_")
			dir.create(base)
			setwd(base)
       	 	ad2bed(x1_snp, y[, ti], stringname = base, plink_exe="/mnt/home/xuliu/Xingjie/plink/plink")
       	 	system("rm -rf *.tfam"); 
        	system("rm -rf *.tped");
        	system(paste0("/mnt/home/xuliu/Xingjie/gemma/bin/gemma -bfile ", base, " -bslmm 1 -w 5000 -s 5000 -o BSLMM"))
        	param = read.table("./output/BSLMM.param.txt", header = T)
        	setwd("..")
        	system(paste0("rm -r ", base))
        	
        	beta_hat_TWAS = param[,5] + param[,6]*param[,7]
        	impute_M <- x2 %*%  (beta_hat_TWAS*apply(x1_snp, 2, sd)) # scale back
        	if (var(impute_M) == 0) {
        		pvalue_TWAS[i, ti] <- 1
        	   	} else {
        		bet_lm <- coefficients(summary(lm(z~impute_M)))
				pvalue_TWAS[i, ti] <- bet_lm[2,4]	
			}
			setwd("..")
		}

 	}
 	toc <- Sys.time()
 	print(toc - tic)
 	saveRDS(list(pvalue_mammotpart = pvalue_mammotpart,
 				 pvalue_comm = pvalue_comm,
 				 pvalue_predixcan = pvalue_predixcan,
 				 pvalue_TWAS = pvalue_TWAS),
 			paste0(path, ".rds")
 	        )

	 

}
