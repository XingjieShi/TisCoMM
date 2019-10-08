//  Created by Xingjie Shi, 24/06/2019.
//  Copyright © 2019 Xingjie Shi. All rights reserved.
// NOTE: the standardization methods of X1 for mammot and mammotSS are different.
//
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "paral.hpp"

using namespace std;
using namespace arma;

void parGene::loop_by_gene(int i){

	uvec idx = find(snp_info.col(0) < expr_info(i,2) + bw && snp_info.col(0) > expr_info(i,1) - bw
        && snp_info.col(1) == expr_info(i,0));

	vec out_param_par = vec(out_param.colptr(i), 4, false); 
	out_param_par(3) = idx.n_elem;

    if (idx.n_elem > 1){
        //cout << "Error: no matching SNPs for " << i << "-th gene ... " << endl;
    //} else {

		mat y = expr.subcube(span(i), span::all, span::all);
 		if (y.n_rows == 1) y = y.t();
		y = y - Prjy * y;
		int T = y.n_cols;

		vec Alpha_par = vec(Alpha.colptr(i), T, false);
 
		mat X1tmp = conv_to<mat>::from(X1.cols(idx)); 
		mat X2tmp = conv_to<mat>::from(X2.cols(idx));


		rowvec meanX1tmp = mean(X1tmp, 0);
		rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
		rowvec meanX2tmp = mean(X2tmp, 0);
		rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual

    
		X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1)) / repmat(sdX1tmp, X1tmp.n_rows, 1)/ sqrt(X1tmp.n_cols);

		X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);

		rowvec X1row = X1tmp.row(0);
    
   		uvec idx3 = find_finite(X1row);
		X1tmp = X1tmp.cols(idx3);
		X2tmp = X2tmp.cols(idx3);
    	if (idx3.n_elem == 0){
     	   cout << "Error: the number of SNPs are 0 for " << i << "-th gene ... " << endl;
    	}

		int p = X1tmp.n_cols, maxIter = 2000;

		//H0: constr = true; produce initial values for Ha, and also hertitability estimate in EQTL. 
    	int IterationH0 = 1;
    	vec LogLikH0 = zeros<vec>(maxIter);

    	mat w = ones<mat>(X1tmp.n_cols, 1);
		if (T != 1) {
			w = cor(X1tmp, y);
			rowvec sdy = stddev(y, 0, 0); // see manual
			colvec sdx = vectorise(stddev(X1tmp, 0, 0)); 
			w.each_row() %= sdy; 
			w.each_col() /= sdx;	
		}

		mat B_hat = zeros<mat>(p, T);
        vec mu = zeros<vec>(p);
	    mat Sigma = eye<mat>(p, p);
	    mat Ve = eye<mat>(T, T);
		Ve.diag() = vectorise(var(y)); 

 		double sigmab = 1.0;
 		double sigmaz = var(z);
		vec alpha = zeros<vec>(T);

 		PXem(X1tmp, y,  X2tmp, z, w, B_hat, mu, Sigma, sigmab, Ve, alpha, sigmaz, LogLikH0, true, true, IterationH0, maxIter);
 

		//H1: constr = false; 
		int IterationH1 = 1;
		vec LogLikH1 = zeros<vec>(maxIter);
    	PXem(X1tmp, y,  X2tmp, z, w, B_hat, mu, Sigma, sigmab, Ve, alpha, sigmaz, LogLikH1, false, true, IterationH1, maxIter);
    	Alpha_par = alpha;

    	out_param_par(0) = sigmaz;
		out_param_par(1) = 2 * (LogLikH1(IterationH1) - LogLikH0(IterationH0));
		out_param_par(2) = sigmab;
				
	}
	if ( i % 100 == 0 && i != 0)	{
					cout << i + 1 << "-th Gene starts working ..." << endl;
	}

}

std::mutex _mtx2;
int parGene::next(){
    std::lock_guard<std::mutex> lockGuard(_mtx2);
    if(current_idx >= (int)Ngene){
        return -1;
    }
    current_idx++;
    return current_idx-1;
}

void parGene::update_by_thread(int thread_id){
    while(true){
        int idx = next();
        if(idx == -1){
            break;
        }
        // cout << idx << endl;
        loop_by_gene(idx);
    }
}


void parGeneSS::loop_by_gene(int i){

	uvec idx = find(snp_info.col(0) < expr_info(i,2) + bw && snp_info.col(0) > expr_info(i,1) - bw
        && snp_info.col(1) == expr_info(i,0));

	vec out_param_par = vec(out_param.colptr(i), 4, false); 
	out_param_par(3) = idx.n_elem;

    if (idx.n_elem > 1){
        //cout << "Error: no matching SNPs for " << i << "-th gene ... " << endl;
    //} else {

		mat y = expr.subcube(span(i), span::all, span::all);
		//mat y = expr_active;
		if (y.n_rows == 1) y = y.t();
		y = y - Prjy * y;
		int T = y.n_cols;

		vec Alpha_par = vec(Alpha.colptr(i), T, false);
 
		mat X1tmp = conv_to<mat>::from(X1.cols(idx)); 
		mat Xpaneltmp = conv_to<mat>::from(Xpanel.cols(idx));
		mat GWAS_SStmp = GWAS_SS.rows(idx);

		rowvec meanX1tmp = mean(X1tmp, 0);
		rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
		rowvec meanXpaneltmp = mean(Xpaneltmp, 0);
		//rowvec sdXpaneltmp = stddev(Xpaneltmp, 0, 0); // see manual

    	// centering
		X1tmp = X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1);
		Xpaneltmp = Xpaneltmp - repmat(meanXpaneltmp, Xpaneltmp.n_rows, 1);

		rowvec X1row = 1/sdX1tmp;
    
   		uvec idx3 = find_finite(X1row);
		X1tmp = X1tmp.cols(idx3);
		Xpaneltmp = Xpaneltmp.cols(idx3);
		GWAS_SStmp = GWAS_SStmp.rows(idx3);


    	if (idx3.n_elem == 0){
     	   cout << "Error: the number of SNPs are 0 for " << i << "-th gene ... " << endl;
    	}


		
		// prepared summary statistics.
		int p = Xpaneltmp.n_cols;
		vec hatmu = GWAS_SStmp.cols(zeros<uvec>(1));
		vec hats  = GWAS_SStmp.cols(ones<uvec>(1));
		mat R1 = cor(Xpaneltmp);

		// positive definite sparse correlation matrix
		mat R = lam * R1 + (1 - lam) * eye(p, p);
		mat S = diagmat(hats);
		mat SR = S * R;
		mat V = SR * S;
		mat L = chol(V);
		mat Ltinv = inv(L.t());
		mat gam_L = Ltinv * hatmu;
		mat U = Ltinv * SR * diagmat(1.0/hats);

		//H0: constr = true; produce initial values for Ha, and also heritability estimate in EQTL. 
    	int IterationH0 = 1, maxIter = 2000;
    	vec LogLikH0 = zeros<vec>(maxIter);

    	mat w = ones<mat>(X1tmp.n_cols, 1);
		if (T != 1) {
			w = cor(X1tmp, y);
			rowvec sdy = stddev(y, 0, 0); // see manual
			colvec sdx = vectorise(stddev(X1tmp, 0, 0)); 
			w.each_row() %= sdy; 
			w.each_col() /= sdx;	
		}

		mat B_hat = zeros<mat>(p, T);
        vec mu = zeros<vec>(p);
	    mat Sigma = eye<mat>(p, p);
	    mat Ve = eye<mat>(T, T);
		Ve.diag() = vectorise(var(y)); 

 		double sigmab = 1.0;
 		double sigmaz = 1.0;
		vec alpha = zeros<vec>(T);

 		PXem_ss(X1tmp, y,  U, gam_L, w, B_hat, mu, Sigma, sigmab, Ve, alpha, sigmaz, LogLikH0, true, true, IterationH0, 2000);
   		//h_y2_par = 1.0/(1.0+Ve.diag()/Sigma.diag());


		//H1: constr = true; 
		int IterationH1 = 1;
		vec LogLikH1 = zeros<vec>(maxIter);
    	PXem_ss(X1tmp, y,  U, gam_L, w, B_hat, mu, Sigma, sigmab, Ve, alpha, sigmaz, LogLikH1, false, true, IterationH1, 2000);
    	Alpha_par = alpha;

    	out_param_par(0) = sigmaz;
		out_param_par(1) = 2 * (LogLikH1(IterationH1) - LogLikH0(IterationH0)); 
		out_param_par(2) = sigmab;
	}
	if ( i % 100 == 0 && i != 0)	{
					cout << i + 1 << "-th Gene starts working ..." << endl;
	}

}

//std::mutex _mtx3;
int parGeneSS::next(){
    std::lock_guard<std::mutex> lockGuard(_mtx2);
    if(current_idx >= (int)Ngene){
        return -1;
    }
    current_idx++;
    return current_idx-1;
}

void parGeneSS::update_by_thread(int thread_id){
    while(true){
        int idx = next();
        if(idx == -1){
            break;
        }
        // cout << idx << endl;
        loop_by_gene(idx);
    }
}



void parGene_part::loop_by_gene(int i){

	uvec idx = find(snp_info.col(0) < expr_info(i,2) + bw && snp_info.col(0) > expr_info(i,1) - bw
        && snp_info.col(1) == expr_info(i,0));

	vec out_param_par = vec(out_param.colptr(i), 4, false); 
	out_param_par(3) = idx.n_elem;

    if (idx.n_elem > 1){
        //cout << "Error: no matching SNPs for " << i << "-th gene ... " << endl;
    //} else {

		mat y = expr.subcube(span(i), span::all, span::all);
		//mat y = expr_active;
		if (y.n_rows == 1) y = y.t();
		y = y - Prjy * y;
		int T = y.n_cols;

		vec Alpha_par = vec(Alpha.colptr(i), T, false);
		vec chisq_par = vec(chisq.colptr(i), T, false);
 		vec h_y2_par = vec(h_y2.colptr(i), T, false);

		mat X1tmp = conv_to<mat>::from(X1.cols(idx)); 
		mat X2tmp = conv_to<mat>::from(X2.cols(idx));


		rowvec meanX1tmp = mean(X1tmp, 0);
		rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
		rowvec meanX2tmp = mean(X2tmp, 0);
		rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual

    
		X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1)) / repmat(sdX1tmp, X1tmp.n_rows, 1)/ sqrt(X1tmp.n_cols);

		X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);

		rowvec X1row = X1tmp.row(0);
    
   		uvec idx3 = find_finite(X1row);
		X1tmp = X1tmp.cols(idx3);
		X2tmp = X2tmp.cols(idx3);
    	if (idx3.n_elem == 0){
     	   cout << "Error: the number of SNPs are 0 for " << i << "-th gene ... " << endl;
    	}

		int p = X1tmp.n_cols, maxIter = 2000;
    	mat w = ones<mat>(X1tmp.n_cols, 1);
		if (T != 1) {
			w = cor(X1tmp, y);
			rowvec sdy = stddev(y, 0, 0); // see manual
			colvec sdx = vectorise(stddev(X1tmp, 0, 0)); 
			w.each_row() %= sdy; 
			w.each_col() /= sdx;	
		}

	
		//1. H0_t: part test
		//2. H1, using output of H0_t as initial values (more robust); 
 		for (int t = 0; t < T; t++)
 		{
 			vec mu0_t = zeros<vec>(p);
			mat B_hat0_t = zeros<mat>(p, T), Sigma0_t = eye<mat>(p, p); 
			mat Ve0_t = eye<mat>(T, T);
			Ve0_t.diag() = vectorise(var(y));
			double sigmab0_t = 1.0, sigmaz = 1.0;
			vec constrFactor = zeros<vec>(T);
			constrFactor(t) = 1;
		    int IterationH0_t = 1, IterationH1 = 1;
			vec LogLikH0_t = zeros<vec>(maxIter), LogLikH1 = zeros<vec>(maxIter);
			vec alpha0_t = zeros<vec>(T);
	
    		PXem_part(X1tmp, y,  X2tmp, z, w, 
    					B_hat0_t, mu0_t, Sigma0_t, sigmab0_t, Ve0_t, alpha0_t, 
    					sigmaz, LogLikH0_t, constrFactor, true, IterationH0_t, maxIter);

    		constrFactor = zeros<vec>(T);
    		PXem_part(X1tmp, y,  X2tmp, z, w, 
    					B_hat0_t, mu0_t, Sigma0_t, sigmab0_t, Ve0_t, alpha0_t, 
    					sigmaz, LogLikH1, constrFactor, true, IterationH1, maxIter);
 			chisq_par(t) = 2 * (LogLikH1(IterationH1) - LogLikH0_t(IterationH0_t));

			//estimate heritability.
			vec covar = ones<vec>(X1tmp.n_rows);
 			List lmm = lmm_pxem(y.col(t), covar, X1tmp, maxIter);
 			double sigma2y = lmm["sigma2y"], sigma2beta= lmm["sigma2beta"];
 			h_y2_par(t) = 1 / (1 + sigma2y/sigma2beta );
 		}	
	}
	if ( i % 10 == 0 && i != 0)	{
					cout << i + 1 << "-th Gene starts working ..." << endl;
	}

}

int parGene_part::next(){
    std::lock_guard<std::mutex> lockGuard(_mtx2);
    if(current_idx >= (int)Ngene){
        return -1;
    }
    current_idx++;
    return current_idx-1;
}

void parGene_part::update_by_thread(int thread_id){
    while(true){
        int idx = next();
        if(idx == -1){
            break;
        }
        // cout << idx << endl;
        loop_by_gene(idx);
    }
}


void parGeneSS_part::loop_by_gene(int i){

	uvec idx = find(snp_info.col(0) < expr_info(i,2) + bw && snp_info.col(0) > expr_info(i,1) - bw
        && snp_info.col(1) == expr_info(i,0));

	vec out_param_par = vec(out_param.colptr(i), 4, false); 
	out_param_par(3) = idx.n_elem;

    if (idx.n_elem > 1){
        //cout << "Error: no matching SNPs for " << i << "-th gene ... " << endl;
    //} else {

		mat y = expr.subcube(span(i), span::all, span::all);
		//mat y = expr_active;
		if (y.n_rows == 1) y = y.t();
		y = y - Prjy * y;
		int T = y.n_cols;

		vec chisq_par = vec(chisq.colptr(i), T, false);
		vec Alpha_par = vec(Alpha.colptr(i), T, false);
		vec h_y2_par = vec(h_y2.colptr(i), T, false);

		mat X1tmp = conv_to<mat>::from(X1.cols(idx)); 
		mat Xpaneltmp = conv_to<mat>::from(Xpanel.cols(idx));
		mat GWAS_SStmp = GWAS_SS.rows(idx);

		rowvec meanX1tmp = mean(X1tmp, 0);
		rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
		rowvec meanXpaneltmp = mean(Xpaneltmp, 0);
		//rowvec sdXpaneltmp = stddev(Xpaneltmp, 0, 0); // see manual

    	// centering
		X1tmp = X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1);
		Xpaneltmp = Xpaneltmp - repmat(meanXpaneltmp, Xpaneltmp.n_rows, 1);
		mat X1_norm = X1tmp / repmat(sdX1tmp, X1tmp.n_rows, 1)/ sqrt(X1tmp.n_cols);//for heritability estimation.

		rowvec X1row = 1/sdX1tmp;
    
   		uvec idx3 = find_finite(X1row);
		X1tmp = X1tmp.cols(idx3);
		X1_norm = X1_norm.cols(idx3);
		Xpaneltmp = Xpaneltmp.cols(idx3);
		GWAS_SStmp = GWAS_SStmp.rows(idx3);


    	if (idx3.n_elem == 0){
     	   cout << "Error: the number of SNPs are 0 for " << i << "-th gene ... " << endl;
    	}


		
		// prepared summary statistics.
		int p = Xpaneltmp.n_cols, maxIter = 2000;
		vec hatmu = GWAS_SStmp.cols(zeros<uvec>(1));
		vec hats  = GWAS_SStmp.cols(ones<uvec>(1));
		mat R1 = cor(Xpaneltmp);

		// positive definite sparse correlation matrix
		mat R = lam * R1 + (1 - lam) * eye(p, p);
		mat S = diagmat(hats);
		mat SR = S * R;
		mat V = SR * S;
		mat L = chol(V);
		mat Ltinv = inv(L.t());
		mat gam_L = Ltinv * hatmu;
		mat U = Ltinv * SR * diagmat(1.0/hats);
    	mat w = ones<mat>(X1tmp.n_cols, 1);
		if (T != 1) {
			w = cor(X1tmp, y);
			rowvec sdy = stddev(y, 0, 0); // see manual
			colvec sdx = vectorise(stddev(X1tmp, 0, 0)); 
			w.each_row() %= sdy; 
			w.each_col() /= sdx;	
		}


		//1. H0_t: part test
		//2. H1, using output of H0_t as initial values (more robust); 
 		for (int t = 0; t < T; t++)
 		{
 			vec mu0_t = zeros<vec>(p);
			mat B_hat0_t = zeros<mat>(p, T), Sigma0_t = eye<mat>(p, p); 
			mat Ve0_t = eye<mat>(T, T);
			Ve0_t.diag() = vectorise(var(y));
			double sigmab0_t = 1.0, sigmaz = 1.0;;
			vec constrFactor = zeros<vec>(T);
			constrFactor(t) = 1;
		    int IterationH0_t = 1, IterationH1 = 1;
			vec LogLikH0_t = zeros<vec>(maxIter), LogLikH1 = zeros<vec>(maxIter);
			vec alpha0_t = zeros<vec>(T);
	
    		PXem_ss_part(X1tmp, y,  U, gam_L, w, 
    					B_hat0_t, mu0_t, Sigma0_t, sigmab0_t, Ve0_t, alpha0_t, 
    					sigmaz, LogLikH0_t, constrFactor, true, IterationH0_t, maxIter);

    		constrFactor = zeros<vec>(T);
    		PXem_ss_part(X1tmp, y,  U, gam_L, w, 
    					B_hat0_t, mu0_t, Sigma0_t, sigmab0_t, Ve0_t, alpha0_t, 
    					sigmaz, LogLikH1, constrFactor, true, IterationH1, maxIter);
 			chisq_par(t) = 2 * (LogLikH1(IterationH1) - LogLikH0_t(IterationH0_t));

 			//estimate heritability.
			vec covar = ones<vec>(X1tmp.n_rows);
 			List lmm = lmm_pxem(y.col(t), covar, X1_norm, maxIter);
 			double sigma2y = lmm["sigma2y"], sigma2beta= lmm["sigma2beta"];
 			h_y2_par(t) = 1 / (1 + sigma2y/sigma2beta );
 		}	
	}
	if ( i % 100 == 0 && i != 0)	{
					cout << i + 1 << "-th Gene starts working ..." << endl;
	}

}

int parGeneSS_part::next(){
    std::lock_guard<std::mutex> lockGuard(_mtx2);
    if(current_idx >= (int)Ngene){
        return -1;
    }
    current_idx++;
    return current_idx-1;
}

void parGeneSS_part::update_by_thread(int thread_id){
    while(true){
        int idx = next();
        if(idx == -1){
            break;
        }
        // cout << idx << endl;
        loop_by_gene(idx);
    }
}
