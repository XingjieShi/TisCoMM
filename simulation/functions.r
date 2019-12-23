
library(glmnet)
predixcan <- function(x1, y, x2, z)	{
	Ti <- ncol(y)
	p <- ncol(x1)
	pvalue <- numeric(Ti)
	for (k in 1:Ti) {
		fit.elasnet.cv <- cv.glmnet(x1, y[,k], type.measure="mse", alpha=0.5, family="gaussian")
		elasnet_M <- predict(fit.elasnet.cv, s=fit.elasnet.cv$lambda.min, newx=x2)
		hat <- coefficients(summary(lm(z~elasnet_M)))
		pvalue[k] <- ifelse(var(elasnet_M) < 1E-6, 1, hat[2,4])
	}
	pvalue
}





# transform Normal x to SNPs.
x2snp = function (X)
{
    n <- nrow(X)
    p <- ncol(X)
    maf <- runif(p, 0.05, 0.5)
    AAprob = maf^2;
    Aaprob = 2 * maf * (1-maf);
    quanti = cbind(1-Aaprob-AAprob, 1- AAprob)  ## attention
    snp = matrix(0, n, p);
    for (j in 1:p){
        cutoff = qnorm(quanti[j,]);
        snp[X[,j] <  cutoff[1], j] = 0;
        snp[X[,j] >= cutoff[1]  & X[,j] < cutoff[2],j] = 1;  ## attention
        snp[X[,j] >= cutoff[2], j] = 2;
    }
    snp
}



# check some equations for block matrix
AR = function(rho, p) {
    outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))) 
}

# adaptive weighting based on marginal regression.
# SNP matrix, GE matrix of all tissues.
adw <- function(x, y){
	w <- matrix(0, nrow = ncol(x), ncol = ncol(y))
	for (i in 1:ncol(x))
		{
			for (j in 1:ncol(y))
			{
				w[i, j] <- coefficients(lm(y[,j]~x[,i]))[2]
			}
		}
	w
}



###################
####################
# multiple predixcan (individual data).	
 	MulTiXcan <- function(x1, y, x2, z, cond = 30)
 	{
        n <- nrow(x2)
        Ti <- ncol(y)
  		yhat2 <- matrix(0, nrow = n, ncol = Ti)	
  		for (ti in 1:Ti)	{
  		  	cvfit = cv.glmnet(x1, y[, ti], alpha = 0.5)
            beta_hat = coef(cvfit, s = "lambda.min")
              
            yhat2[, ti] = x2%*%beta_hat[-1] + beta_hat[1]
        }
        # remove zero columns.
        yhat2_var <- apply(yhat2, 2, var)
        if (all(yhat2_var == 0))	{
        	alpha <- 0
        	pval <- 1
        }	else  if (sum(yhat2_var != 0) == 1) {
        	twas <- summary(lm(z ~ yhat2))
            alpha <- twas$coefficients[2,1]
			pval <- twas$coefficients[2,4]
        } 	else  {
        	yhat2 <- yhat2[, yhat2_var != 0]

        	yhat2_pca <- summary(prcomp(yhat2, center = TRUE, scale. = TRUE))
        	ind_top <- which(yhat2_pca$importance[2, 1]/yhat2_pca$importance[2, ] < cond)
        	pc_top <- yhat2_pca$x
        	twas <- summary(lm(z ~ pc_top))
        	alpha <- twas$coefficients[-1]
        	pval <- pf(twas$fstatistic[1],twas$fstatistic[2],twas$fstatistic[3],lower.tail=FALSE)
    	}
    	list(alpha = alpha, 
    		 pval  =  unname(pval))
    }
      

# prepare binary ped files for GEMMA/Plink.
# imput:
#       x          --- simulated design matrix, enties are 
#                      counts of minor alleles.
#       y          --- mulvariate outcomes
#       stringname --- prefix for filenames
#       plink_exe  --- location of plink.exe
# output:
#       write files fro analysis in Gemma/Plink.
ad2bed = function(x, y, stringname, plink_exe="/home/sysadmin/Software/plink/plink")
{   
    n <- nrow(x)
    p <- ncol(x)
    snp <- ad2mm(x) # its body is in vb_aux.cpp.
    tped <- cbind(rep(1, p), 
                  paste0("snp", 1:p), 
                  rep(0, p),
                  rep(0, p),
                  snp
                 )
    write.table(tped, quote = FALSE, 
                row.names = FALSE, col.names = FALSE,
                paste0(stringname, ".tped"))
    tfam = data.frame(paste0("ID", 1:n), 
                      paste0("ID", 1:n), 
                      paste0("PD", 1:n),  
                      paste0("MD", 1:n),  
                      rep(1, n),
                      rep(-9,n))
    write.table(tfam, quote = FALSE, 
                row.names = FALSE, col.names = FALSE,
                paste0(stringname, ".tfam"))
    cmd = paste0(plink_exe, " --tfile ", stringname, " --noweb --make-bed --out ", stringname)
    system(cmd)

    # fill y in ".fam" files
    famfile <- read.table(paste0(stringname, ".fam"), 
                          stringsAsFactors=FALSE, header=FALSE)
    famfile$V6 <- NULL
    famfileK <- cbind(famfile, y)
    write.table(famfileK, quote = FALSE, 
                row.names = FALSE, col.names = FALSE,
                paste0(stringname, ".fam"))

}

