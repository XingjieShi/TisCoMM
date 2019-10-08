#ifndef mammot_part_hpp
#define mammot_part_hpp

#include <RcppArmadillo.h>
 
//#include <Rcpp.h>
//#include <omp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat ADW(arma::mat& x, arma::mat& y);

void	PXem_part(const arma::mat& x1, const arma::mat& y, const arma::mat& x2, const arma::vec& z, const arma::mat& w, arma::mat& B_hat, arma::vec& mu, arma::mat& Sigma,
	 	  double& sigmab, arma::mat& Ve, arma::vec& alpha, double& sigmaz,
	 	  arma::vec& LogLik, const arma::vec& constrFactor, const bool& PX, int& Iteration, const int& maxIter);
void  PXem_ss_part(const arma::mat& x1, const arma::mat& y, const arma::mat& x2, const arma::vec& z, const arma::mat& w, arma::mat& B_hat, arma::vec& mu, arma::mat& Sigma,
	 	  double& sigmab, arma::mat& Ve, arma::vec& alpha, double& sigmaz,
	 	  arma::vec& LogLik, const arma::vec& constrFactor, const bool& PX, int& Iteration, const int& maxIter);

List	mammot_PXem_part(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter);
List  mammotSS_PXem_part(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter);

List	mammot_part_test(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter);
List  mammotSS_part_test(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					arma::vec constrFactor, bool PX, int maxIter);

void lmm_pxem_ptr2(const arma::vec& y, const arma::mat& W, const arma::mat& X,  const int& maxIter,
              double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
              int& iteration, arma::mat& Sigb, arma::vec& mub);


List lmm_pxem(const arma::vec y, const arma::mat w, const arma::mat x, const int maxIter);

#endif /* mammot_part_hpp */
