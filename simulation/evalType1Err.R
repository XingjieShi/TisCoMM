# simulation I
# Type I effor rate evluation under NULL models
library(mvtnorm)
library(GBJ)
library(Rcpp)
source("functions.r")  
source("multiEstB.R")
sourceCpp("PXem_ss.cpp")
sourceCpp("PXem.cpp")

evalType1Err <- function(h_z, h_c, s, rhoX, rhoE, Ti, batch, Rep = 1000)
{
	PX <- T
	tic <- Sys.time()
	pvalue <- matrix(0, ncol=5, nrow=Rep)
 	# mammot, mammot_ss, multixcan, multixcan_ss, utmost.
	for(i in 1:Rep)	{
		m = 400
		n = 5000
		p = 300

		# genotypes in eQTL
		x1 <- x2snp(rmvnorm(m, mean=rep(0, p), sigma = AR(rhoX, p)))
		x1mean <- matrix(rep(colMeans(x1), m), ncol=p, byrow = T)
		x1 <- (x1 - x1mean)

		# genotypes in GWAS
		x2 <- x2snp(rmvnorm(n, mean=rep(0, p), sigma = AR(rhoX, p)))
		x2mean <- matrix(rep(colMeans(x2), n), ncol=p, byrow = T)
		x2 <- (x2 - x2mean)
 
 		# eQTL coefficient matrix, across Ti tissues.
		b <- rnorm(p)
		W <- matrix(0, nrow = p, ncol=Ti)
		for (j in 1:Ti)	
			W[sample(p, s*p, replace=FALSE),j] <- 1
  
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
			a <- runif(Ti, -1, 1)
			lp_z <- x2 %*% B %*% a
			#h_z <- .01
			sig2 <- c(var(lp_z)) * ((1-h_z)/h_z)
			z <- lp_z + rnorm(n, 0, sqrt(sig2))
  		}

 		maxIter = 3000
 		# mammot (individual data)
 		# adaptive weight for mammot and mammot_ss.
		w <- matrix(1, nrow=p, ncol=1)
		if(Ti != 1) w = adw(x1, y)  
 		fit_PX <- mammot_PXem(x1, y, x2, z, w, FALSE, PX, maxIter=maxIter)
		fit_PX_0 <- mammot_PXem(x1, y, x2, z, w, TRUE, PX, maxIter=maxIter)
		chisq <- 2*(max(fit_PX$loglik) - max(fit_PX_0$loglik))
		pvalue[i, 1] <- pchisq(chisq, Ti, lower.tail=F)

  		# reference panal
  		n_p <- 400
		x3 <- rmvnorm(n_p, mean=rep(0, p), sigma = AR(rhoX, p))
		x3mean <- matrix(rep(colMeans(x3), n_p), ncol=p, byrow = T)
		#x3sd <- matrix(rep(apply(x3, 2, sd), n_p), ncol=p, byrow = T)
		x3 <- (x3 - x3mean)#/x3sd/sqrt(n_p)

		 

		lam = 0.95
 		sumx3 = apply(x3*x3, 2, sum)
		RR = matrix(0, p, p);
		for (i1 in 1:p){
    		for (j1 in 1:p){
        		RR[i1, j1] = t(x3[,i1])%*%x3[,j1]/sqrt(sumx3[i1]*sumx3[j1])                        
    		}
		}
 
		R = RR*lam + (1 - lam)*diag(p)

	
		# summary statisitcs for GWAS	
 		hatmu = matrix(0, p, 1)
		hats  = matrix(0, p, 1)
		for (j in 1:p){
    		fm <- lm(z ~ 1 + x2[, j]);
    		hatmu[j] <- summary(fm)$coefficients[2,1]
    		hats[j]  <- summary(fm)$coefficients[2,2];
		}

		# mammot_ss (GWAS summary data)
		S <- diag(c(hats))
		SR <- S %*% R 
		V <- SR %*% S
		L <- chol(V)
		Ltinv <- solve(t(L))
		gam_L <- Ltinv %*% hatmu
		U <- Ltinv %*% SR %*% solve(S)

		fit_ss   <- mammot_PXem_ss(x1, y, U, gam_L, w, FALSE, PX, maxIter=maxIter)
		fit_ss_0 <- mammot_PXem_ss(x1, y, U, gam_L, w, TRUE,  PX, maxIter=maxIter)
		chisq  <- 2*(max(fit_ss$loglik) - max(fit_ss_0$loglik))
		pvalue[i, 2] <- pchisq(chisq, Ti, lower.tail=F)
	

		# multiXcan (individual data)
		cond <- 30
 		pvalue[i, 3] <- MulTiXcan(x1, y, x2, z, cond)$pval	

		# utmost (summary data)
		B_hat = multiEstB(x1, y)
		if (sum(abs(B_hat)) > 1e-6) {
			B_hat <- as.matrix(B_hat[, apply(B_hat, 2, sd)!= 0])
			Ge_impute <- x3 %*% B_hat
			eta <- apply(Ge_impute, 2, sd)
			#sig_z <- sd(z)
			sig_x <- sqrt(diag(var(x3)))
			Ztilde <- hatmu/hats
			Lambda <-  diag(sig_x) %*% sweep(B_hat, 2, eta, "/")
			Lambda_sub <- Lambda[, which(eta != 0)]
			test_stats <- t(Lambda_sub) %*% Ztilde
			cov_Z <- t(Lambda_sub) %*% R %*% Lambda_sub
			pvalue[i, 5] <- GBJ(test_stats=test_stats, cor_mat=cov_Z)$GBJ_pvalue

			# multiXcan (summary data)
 			rhoGE <- cor(Ge_impute)
			rhoGE_svd <- svd(rhoGE)
			ind_top <- which(rhoGE_svd$d[1]/rhoGE_svd$d < cond)
			if(length(ind_top) == 0) ind_top <- 1
			u <- rhoGE_svd$u
			v <- rhoGE_svd$v
			us <- as.matrix(u[, ind_top])
			vs <- as.matrix(v[, ind_top])
			d <- rhoGE_svd$d[ind_top]
			if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
			rhoGE_ginv <- vs %*% ds %*% t(us)
			chisq <- c(t(test_stats) %*% rhoGE_ginv %*% test_stats)
			pvalue[i, 4] <- pchisq(chisq,  length(ind_top), lower.tail=F)
			
		} else { 
			pvalue[i, 5] <- NA
			pvalue[i, 4] <- NA
		}	
		

 		
 	}
 	toc <- Sys.time()
 	print(toc - tic)

 	colnames(pvalue) <- c("mammot", "mammot_ss", "multixan", "multixcan_ss", "utmost")
	write.table(pvalue, paste0("pvalue_hz", h_z*10, "_hc", h_c*10, "_rhoX", rhoX*10, "_s", s*10, "_batch-", batch, ".txt"), row.names = F, col.names = F)

}
