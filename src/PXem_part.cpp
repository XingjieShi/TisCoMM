#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat ADW(const arma::mat& x, const arma::mat& y)	{
		mat w = cor(x, y);
				
		rowvec sdy = stddev(y, 0, 0); // see manual
		colvec sdx = vectorise(stddev(x, 0, 0)); 
		w.each_row() %= sdy; 
		w.each_col() /= sdx;		
		return w;
}

void	PXem_part(const arma::mat& x1, const arma::mat& y, const arma::mat& x2, const arma::vec& z, const arma::mat& w, arma::mat& B_hat, arma::vec& mu, arma::mat& Sigma,
	 	  double& sigmab, arma::mat& Ve, arma::vec& alpha, double& sigmaz,
	 	  arma::vec& LogLik, const arma::vec& constrFactor, const bool& PX, int& Iteration, const int& maxIter)
{
	int m = y.n_rows, T = y.n_cols, p = x1.n_cols, n = z.n_elem;
	mat I_p = eye<mat>(p, p);
	// Parameter Expansion
	double lam = 1.0;

	// precomputation
	mat x1tx1 = x1.t() * x1;
	mat x2tx2 = x2.t() * x2;
	mat ytx1 = y.t() * x1;
	mat x2tz = x2.t() * z;
 	mat Sigma_i = I_p;

	mat yhat = x1 * B_hat;
	mat residual_y = y - yhat;

	mat ztilde = x2 * B_hat;
	vec residual_z = z - ztilde * alpha;

	//mat wv_i = w * v_i;
	mat v_iwt = solve(Ve, w.t());
	mat wv_i = v_iwt.t();
	mat wv_iwt = wv_i * w.t();
	mat wa = w * alpha;

	mat Lv = Ve;
	mat Ls = I_p;
   	LogLik(0) = -INFINITY;
   	int iter;
   	uvec idx = find(constrFactor == 0);
	//--------------------------------------------------------------------------------	
	// EM algrithm
	//--------------------------------------------------------------------------------
	
	for (iter = 1; iter < maxIter; iter ++ ) {

		// E-step
		Sigma_i = lam *lam * (x1tx1 % wv_iwt) + x2tx2 % (wa * wa.t())/sigmaz + I_p/sigmab;		
 	    Sigma = solve(Sigma_i, I_p);
	    mat x1tx1S = x1tx1 % Sigma;

	    mu = Sigma * (lam * diagvec(wv_i * ytx1) + diagmat(x2tz) * wa /sigmaz);
		B_hat = diagmat(mu) * w;
	    //cout << "diagmat(mu)  = " << mu  << endl;
		// M-step 
		// residual
		yhat = lam * x1 * B_hat;
		residual_y = y - yhat;
        ztilde = x2 * B_hat;
		residual_z = z - ztilde * alpha;

		// v
		Ve = residual_y.t()*residual_y/m + lam*lam/m * w.t() * x1tx1S * w; 
		//wv_i = w * v_i;
		v_iwt = solve(Ve, w.t());
	    wv_i = v_iwt.t();
		wv_iwt = wv_i * w.t();

		// sigmab
	 	sigmab = (sum(mu % mu) + trace(Sigma))/p;

		// sigmaz
		mat sigmaz_2 = w.t() * (x2tx2 % Sigma) * w;
        sigmaz = sum(residual_z % residual_z)/n + trace(sigmaz_2 * alpha * alpha.t())/n;

		// alpha 
		mat w_j = w.cols(idx);
        mat alpha_d_1 = w_j.t() * (x2tx2 % (mu * mu.t())) * w_j;
        mat alpha_d_2 = w_j.t() * (x2tx2 % Sigma) * w_j;
        mat B_hat_j = B_hat.cols(idx);
        vec alpha_j = solve(alpha_d_1 + alpha_d_2, B_hat_j.t()) * x2tz;

        alpha(idx) = alpha_j;
        //cout << alpha << endl;

		vec residual_z_new = z - ztilde * alpha;
		double n_sigmaz_new = sum(residual_z_new % residual_z_new) + trace(sigmaz_2 * alpha * alpha.t());
		wa = w * alpha;
        
		// Gamma
		if (PX)	{
			double lam2 = sum(mu.t() * (x1tx1 % wv_iwt) * mu) + trace(wv_i.t() * x1tx1S * w);
			lam = sum(mu % diagvec(wv_i * ytx1))/lam2;
		} 

		// current log-likelihood
		Lv = chol(Ve, "lower");
		Ls = chol(Sigma, "lower");

		LogLik(iter) = -(m*T + n) * 0.5 *  log(2*datum::pi) - m * sum(log(Lv.diag())) - m*T/2 - 
	 						n * 0.5 * log(sigmaz) - n_sigmaz_new/2/sigmaz -
	 						p * 0.5 * log(sigmab) +
	 						 sum(log(Ls.diag()));

		if ( LogLik(iter) - LogLik(iter - 1) < -1e-8 ){
			perror("The likelihood failed to increase!");
			break;
		}

		if (abs(LogLik(iter) - LogLik(iter - 1)) < 1e-5 || rcond(Sigma_i) < 1e-6) {
			break;
		}
	}

	if (iter == maxIter) {
		Iteration = iter - 1;
	} else {
		Iteration = iter;
	}
	// Reduce-step
	sigmab = sigmab * lam * lam;
	alpha = alpha/lam;
}


// PXem for summary statisitcs
void	PXem_ss_part(const arma::mat& x1, const arma::mat& y, const arma::mat& x2, const arma::vec& z, const arma::mat& w, arma::mat& B_hat, arma::vec& mu, arma::mat& Sigma,
	 	  double& sigmab, arma::mat& Ve, arma::vec& alpha, double& sigmaz,
	 	  arma::vec& LogLik, const arma::vec& constrFactor, const bool& PX, int& Iteration, const int& maxIter)
{
	int m = y.n_rows, T = y.n_cols, p = x1.n_cols, n = z.n_elem;
	mat I_p = eye<mat>(p, p);
	// Parameter Expansion
	double lam = 1.0;

	// precomputation
	mat x1tx1 = x1.t() * x1;
	mat x2tx2 = x2.t() * x2;
	mat ytx1 = y.t() * x1;
	mat x2tz = x2.t() * z;
 	mat Sigma_i = I_p;

	mat yhat = x1 * B_hat;
	mat residual_y = y - yhat;

	mat ztilde = x2 * B_hat;
	vec residual_z = z - ztilde * alpha;

	//mat wv_i = w * v_i;
	mat v_iwt = solve(Ve, w.t());
	mat wv_i = v_iwt.t();
	mat wv_iwt = wv_i * w.t();
	mat wa = w * alpha;

	mat Lv = Ve;
	mat Ls = I_p;
   	LogLik(0) = -INFINITY;
   	int iter;
   	uvec idx = find(constrFactor == 0);

	//--------------------------------------------------------------------------------	
	// EM algrithm
	//--------------------------------------------------------------------------------
	
	for (iter = 1; iter < maxIter; iter ++ ) {

		// E-step
		Sigma_i = lam *lam * (x1tx1 % wv_iwt) + x2tx2 % (wa * wa.t())/sigmaz + I_p/sigmab;		
 	    Sigma = solve(Sigma_i, I_p);
	    mat x1tx1S = x1tx1 % Sigma;

	    mu = Sigma * (lam * diagvec(wv_i * ytx1) + diagmat(x2tz) * wa /sigmaz);
		B_hat = diagmat(mu) * w;
	    //cout << "diagmat(mu)  = " << mu  << endl;
		// M-step 
		// residual
		yhat = lam * x1 * B_hat;
		residual_y = y - yhat;
        ztilde = x2 * B_hat;
		residual_z = z - ztilde * alpha;

		// v
		Ve = residual_y.t()*residual_y/m + lam*lam/m * w.t() * x1tx1S * w; 
		//wv_i = w * v_i;
		v_iwt = solve(Ve, w.t());
	    wv_i = v_iwt.t();
		wv_iwt = wv_i * w.t();

		// sigmab
	 	sigmab = (sum(mu % mu) + trace(Sigma))/p;

		// sigmaz
		mat sigmaz_2 = w.t() * (x2tx2 % Sigma) * w;
        //sigmaz = sum(residual_z % residual_z)/n + trace(sigmaz_2 * alpha * alpha.t())/n;
		sigmaz = 1;

		// alpha 
		mat w_j = w.cols(idx);
        mat alpha_d_1 = w_j.t() * (x2tx2 % (mu * mu.t())) * w_j;
        mat alpha_d_2 = w_j.t() * (x2tx2 % Sigma) * w_j;
        mat B_hat_j = B_hat.cols(idx);
        vec alpha_j = solve(alpha_d_1 + alpha_d_2, B_hat_j.t()) * x2tz;

        alpha(idx) = alpha_j;
        //cout << alpha << endl;

		vec residual_z_new = z - ztilde * alpha;
		double n_sigmaz_new = sum(residual_z_new % residual_z_new) + trace(sigmaz_2 * alpha * alpha.t());
		wa = w * alpha;
        
		// Gamma
		if (PX)	{
			double lam2 = sum(mu.t() * (x1tx1 % wv_iwt) * mu) + trace(wv_i.t() * x1tx1S * w);
			lam = sum(mu % diagvec(wv_i * ytx1))/lam2;
		} 

		// current log-likelihood
		Lv = chol(Ve, "lower");
		Ls = chol(Sigma, "lower");

		LogLik(iter) = -(m*T + n) * 0.5 *  log(2*datum::pi) - m * sum(log(Lv.diag())) - m*T/2 - 
	 						n * 0.5 * log(sigmaz) - n_sigmaz_new/2/sigmaz -
	 						p * 0.5 * log(sigmab) +
	 						 sum(log(Ls.diag()));

		if ( LogLik(iter) - LogLik(iter - 1) < -1e-8 ){
			perror("The likelihood failed to increase!");
			break;
		}

		if (abs(LogLik(iter) - LogLik(iter - 1)) < 1e-5 || rcond(Sigma_i) < 1e-6) {
			break;
		}
	}

	if (iter == maxIter) {
		Iteration = iter - 1;
	} else {
		Iteration = iter;
	}
	// Reduce-step
	sigmab = sigmab * lam * lam;
	alpha = alpha/lam;
}


// [[Rcpp::export]]
List	mammot_part_test(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter)	{

	if (y.n_rows != x1.n_rows || x1.n_cols != x2.n_cols || z.n_elem != x2.n_rows)	{
		perror("Dimensions of inpute (x1, y, x2, z) are not matched.");
	}
	int T = y.n_cols, p = x1.n_cols;

	mat Ve = eye<mat>(T, T);
	Ve.diag() = vectorise(var(y)); 
	double sigmab = 1;

 	double sigmaz = var(z);

 	vec alpha = zeros<vec>(T); 	
	vec LogLik0 = zeros<vec>(maxIter);
	vec LogLik1 = zeros<vec>(maxIter);

	int Iteration0 = 1, Iteration1 = 1;

	mat B_hat = zeros<mat>(x1.n_cols, y.n_cols);
	vec mu = zeros<vec>(p);
	mat Sigma = eye<mat>(p, p);

	PXem_part(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik0, constrFactor, PX, Iteration0, maxIter);
	vec alpha0 = alpha;

	vec Noconstr = zeros<vec>(T);
	PXem_part(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik1, Noconstr, PX, Iteration1, maxIter);
	vec alpha1 = alpha;
	
	List output = List::create(
					Rcpp::Named("alpha0") = alpha0,
					Rcpp::Named("alpha1") = alpha1,
					Rcpp::Named("chisq") = 2 * (LogLik1(Iteration1) - LogLik0(Iteration0)));

	return output; 
}



// [[Rcpp::export]]
List	mammotSS_part_test(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter)	{

	if (y.n_rows != x1.n_rows || x1.n_cols != x2.n_cols || z.n_elem != x2.n_rows)	{
		perror("Dimensions of inpute (x1, y, x2, z) are not matched.");
	}
	int T = y.n_cols, p = x1.n_cols;

	mat Ve = eye<mat>(T, T);
	Ve.diag() = vectorise(var(y)); 
	double sigmab = 1;
 	double sigmaz = 1;

 	vec alpha = zeros<vec>(T); 	
	vec LogLik0 = zeros<vec>(maxIter);
	vec LogLik1 = zeros<vec>(maxIter);

	int Iteration0 = 1, Iteration1 = 1;

	mat B_hat = zeros<mat>(x1.n_cols, y.n_cols);
	vec mu = zeros<vec>(p);
	mat Sigma = eye<mat>(p, p);

	PXem_ss_part(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik0, constrFactor, PX, Iteration0, maxIter);
	vec alpha0 = alpha;

	vec Noconstr = zeros<vec>(T);
	PXem_ss_part(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik1, Noconstr, PX, Iteration1, maxIter);
	vec alpha1 = alpha;
	

	List output = List::create(
					Rcpp::Named("alpha0") = alpha0,
					Rcpp::Named("alpha1") = alpha1,
					Rcpp::Named("mu") = mu,
					Rcpp::Named("Sigma") = Sigma,
					Rcpp::Named("B_hat") = B_hat,
					Rcpp::Named("Ve") = Ve,
					Rcpp::Named("chisq") = 2 * (LogLik1(Iteration1) - LogLik0(Iteration0)));

	return output; 
}



// [[Rcpp::export]]
List	mammot_PXem_part(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter)	{

	if (y.n_rows != x1.n_rows || x1.n_cols != x2.n_cols || z.n_elem != x2.n_rows)	{
		perror("Dimensions of inpute (x1, y, x2, z) are not matched.");
	}
	int T = y.n_cols, p = x1.n_cols;

	mat Ve = eye<mat>(T, T);
	Ve.diag() = vectorise(var(y)); 
	double sigmab = 1;

 	double sigmaz = var(z);

 	vec alpha = zeros<vec>(T); 	
	vec LogLik = zeros<vec>(maxIter);

	int Iteration = 1;

	mat B_hat = zeros<mat>(x1.n_cols, y.n_cols);
	vec mu = zeros<vec>(p);
	mat Sigma = eye<mat>(p, p);

	PXem_part(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik, constrFactor, PX, Iteration, maxIter);
	
	vec loglik;
	loglik = LogLik.subvec(1, Iteration);

	List output = List::create(
		Rcpp::Named("x1") = x1,
		Rcpp::Named("y") = y,
		Rcpp::Named("x2") = x2,
		Rcpp::Named("z") = z,
		Rcpp::Named("W") = w,
		Rcpp::Named("mu") = mu,
		Rcpp::Named("Sigma") = Sigma,
		Rcpp::Named("B_hat") = B_hat,
		Rcpp::Named("sigmab") = sigmab,
		Rcpp::Named("Ve") = Ve,
		Rcpp::Named("alpha") = alpha,
		Rcpp::Named("sigmaz") = sigmaz,
		Rcpp::Named("loglik") = loglik);

	return output; 
}



// [[Rcpp::export]]
List	mammotSS_part_est(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::mat B_hat, arma::vec mu, arma::mat Sigma, double sigmab, arma::mat Ve, arma::vec alpha, 
					arma::vec constrFactor, bool PX, int maxIter)	{
	double sigmaz = 1;
	vec LogLik0 = zeros<vec>(maxIter);
	int Iteration0 = 1;

	PXem_ss_part(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik0, constrFactor, PX, Iteration0, maxIter);
	
	List output = List::create(
					Rcpp::Named("alpha") = alpha,
					Rcpp::Named("loglik") = LogLik0(Iteration0));

	return output; 
}



// [[Rcpp::export]]
List	mammotSS_PXem_part(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter)	{

	if (y.n_rows != x1.n_rows || x1.n_cols != x2.n_cols || z.n_elem != x2.n_rows)	{
		perror("Dimensions of inpute (x1, y, x2, z) are not matched.");
	}
	int T = y.n_cols, p = x1.n_cols;

	mat Ve = eye<mat>(T, T);
	Ve.diag() = vectorise(var(y)); 
	double sigmab = 1;
 	double sigmaz = 1;

 	vec alpha = zeros<vec>(T); 	
	vec LogLik = zeros<vec>(maxIter);

	int Iteration = 1;

	mat B_hat = zeros<mat>(x1.n_cols, y.n_cols);
	vec mu = zeros<vec>(p);
	mat Sigma = eye<mat>(p, p);

	PXem_ss_part(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik, constrFactor, PX, Iteration, maxIter);

	vec loglik;
	loglik = LogLik.subvec(1, Iteration);

	List output = List::create(
		Rcpp::Named("x1") = x1,
		Rcpp::Named("y") = y,
		Rcpp::Named("x2") = x2,
		Rcpp::Named("z") = z,
		Rcpp::Named("W") = w,
		Rcpp::Named("mu") = mu,
		Rcpp::Named("Sigma") = Sigma,
		Rcpp::Named("B_hat") = B_hat,
		Rcpp::Named("sigmab") = sigmab,
		Rcpp::Named("Ve") = Ve,
		Rcpp::Named("alpha") = alpha,
		Rcpp::Named("sigmaz") = sigmaz,
		Rcpp::Named("loglik") = loglik);

	return output; 
}



// a simple linear mixed model, PXEM algorithm.
void lmm_pxem_ptr2(const arma::vec& y, const arma::mat& W, const arma::mat& X,  const int& maxIter,
              double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
              int& iteration, arma::mat& Sigb, arma::vec& mub){

  int n = y.n_elem, p = X.n_cols;

  if (y.n_elem != X.n_rows || X.n_rows != W.n_rows){
    perror("The dimensions in outcome and covariates (X and W) are not matched");
  }

  if (beta0.n_elem != W.n_cols){
    perror("The dimensions in covariates are not matched in W and beta0");
  }

  if (p != mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }

  if (p != Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }

  mat XtX = X.t()*X, WtW = W.t()*W, WtX = W.t()*X;
  vec Xty = X.t()*y, Wty = W.t()*y;

  vec SWy;
  mat SWX;

  if(W.n_cols==1){
    SWy = mean(y);
    SWX = mean(X,0);
  } else{
    SWy = solve(WtW, Wty);
    SWX = solve(WtW, WtX);
  }


  double gam, gam2;  // parameter expansion

  vec eVal;
  mat eVec;

  eig_sym(eVal, eVec, XtX);

  // initialize
  sigma2y = var(y);
  sigma2beta = sigma2y/p;
  beta0 = SWy - SWX * mub;
  vec loglik(maxIter);
  loglik(0) = -datum::inf;

  vec D;
  vec Xmu;
  vec y_bar = y - W * beta0;
  double y_Xmu2, E;

  iteration = maxIter-1;
  for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
    D = 1 / sigma2beta +  eVal / sigma2y;
    mub = 1/sigma2y * eVec * (eVec.t() * (X.t() * y_bar) / D);
    Xmu = X * mub;
    y_Xmu2 = sum(pow(y_bar-Xmu,2));

    // Evaluate loglik
    E = y_Xmu2/(2*sigma2y) + accu(pow(mub,2))/(2*sigma2beta);
    loglik(iter) = - p*log(sigma2beta)/2 - n*log(sigma2y)/2 - E - sum(log(D))/2 - n/2*log(2*datum::pi);

    if ( loglik(iter) - loglik(iter - 1) < 0 ){
      perror("The likelihood failed to increase!");
    }

    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }

    // M-step
    gam = sum(y_bar % Xmu) / (accu(pow(Xmu,2)) + sum(eVal/D));
    gam2 = pow(gam , 2);

    beta0 = SWy - (SWX * mub) * gam;
    y_bar = y - W * beta0;;

    sigma2y = sum(pow(y_bar-Xmu*gam,2))/n + gam2 * sum(eVal/D)/n;
    sigma2beta = accu(pow(mub,2))/p + sum(1/D)/p;

    // Reduction step
    sigma2beta = gam2 * sigma2beta;
    // gam = 1;
    // gam2 = pow(gam , 2);
  }

  vec loglik_out;
  loglik_out = loglik.subvec(0, iteration);

  loglik_max = loglik(iteration);
}

// [[Rcpp::export]]
Rcpp::List lmm_pxem(const arma::vec y, const arma::mat w, const arma::mat x, const int maxIter){

    double sigma2y = var(y)/2, sigma2beta = var(y)/2, loglik;
    vec beta0 =zeros<vec>(w.n_cols);
    int iter;
    mat Sigb = zeros<mat>(x.n_cols,x.n_cols);
    vec mub  = zeros<vec>(x.n_cols);

    lmm_pxem_ptr2(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);

    List output = List::create(Rcpp::Named("sigma2y") = sigma2y,
                               Rcpp::Named("sigma2beta") = sigma2beta,
                               Rcpp::Named("beta0") = beta0,
                               Rcpp::Named("loglik") = loglik,
                               Rcpp::Named("iteration") = iter,
                               Rcpp::Named("Sigb") = Sigb,
                               Rcpp::Named("mub") = mub);

    return output;

}
