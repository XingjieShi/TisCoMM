#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

 
void	PXem_ss(const arma::mat& x1, const arma::mat& y, const arma::mat& x2, const arma::vec& z, const arma::mat& w, arma::mat& B_hat, arma::vec& mu, arma::mat& Sigma,
	 	  double& sigmab, arma::mat& Ve, arma::vec& alpha, double& sigmaz,
	 	  arma::vec& LogLik, const bool& constr, const bool& PX, int& Iteration, const int& maxIter)
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

	mat v_iwt = solve(Ve, w.t());
	mat wv_i = v_iwt.t();
	mat wv_iwt = wv_i * w.t();
	mat wa = w * alpha;

	mat Lv = Ve;
	mat Ls = I_p;
   	LogLik(0) = -INFINITY;
   	int iter;
	//--------------------------------------------------------------------------------	
	// EM algrithm
	//--------------------------------------------------------------------------------
	
	for (iter = 1; iter < maxIter; iter ++ ) {

		// E-step
		Sigma_i = lam *lam * (x1tx1 % wv_iwt) + x2tx2 % (wa * wa.t())/sigmaz + I_p/sigmab;
 	    Sigma = inv_sympd(Sigma_i);
	    mat x1tx1S = x1tx1 % Sigma;

	    mu = Sigma * (lam * diagvec(wv_i * ytx1) + diagmat(x2tz) * wa /sigmaz);
		B_hat = diagmat(mu) * w;
		// M-step 
		// residual
		yhat = lam * x1 * B_hat;
		residual_y = y - yhat;
        ztilde = x2 * B_hat;
		residual_z = z - ztilde * alpha;

		// v
		Ve = residual_y.t()*residual_y/m + lam*lam/m * w.t() * x1tx1S * w; 
		v_iwt = solve(Ve, w.t());
	    wv_i = v_iwt.t();
		wv_iwt = wv_i * w.t();

		// sigmab
	 	sigmab = (sum(mu % mu) + trace(Sigma))/p;

		// sigmaz
		mat alpha_d_2 = w.t() * (x2tx2 % Sigma) * w;
        //sigmaz = sum(residual_z % residual_z)/n + trace(alpha_d_2 * alpha * alpha.t())/n;
		sigmaz = 1.0;
		// alpha 
        mat alpha_d_1 = w.t() * (x2tx2 % (mu * mu.t())) * w;
        if (!constr) alpha = solve(alpha_d_1 + alpha_d_2, B_hat.t()) * x2tz;
		vec residual_z_new = z - ztilde * alpha;
		double n_sigmaz_new = sum(residual_z_new % residual_z_new) + trace(alpha_d_2 * alpha * alpha.t());
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
		}

		if (abs(LogLik(iter) - LogLik(iter - 1)) < 1e-5 || rcond(Sigma_i) < 1e-7) {
			
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
List	mammot_PXem_ss(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					bool constr, bool PX, int maxIter)	{

	if (y.n_rows != x1.n_rows || x1.n_cols != x2.n_cols || z.n_elem != x2.n_rows)	{
		perror("Dimensions of inpute (x1, y, x2, z) are not matched.");
	}
	int T = y.n_cols, p = x1.n_cols;

	mat Ve = eye<mat>(T, T);
	Ve.diag() = vectorise(var(y)); 
	double sigmab = 1;

 	double sigmaz = 1.0;

 	vec alpha = zeros<vec>(T);

	vec LogLik = zeros<vec>(maxIter);
	int Iteration;

	mat B_hat = zeros<mat>(x1.n_cols, y.n_cols);
	vec mu = zeros<vec>(p);
	mat Sigma = eye<mat>(p, p);

	PXem_ss(x1, y, x2, z, w, B_hat, mu, Sigma,
		 sigmab, Ve, alpha, sigmaz,
		 LogLik, constr, PX, Iteration, maxIter);


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