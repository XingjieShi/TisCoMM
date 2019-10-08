#ifndef PXem_ss_hpp
#define PXem_ss_hpp

#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>
#include <math.h>


using namespace Rcpp;
using namespace arma;
using namespace std;



void	PXem_ss(const arma::mat& x1, const arma::mat& y, const arma::mat& x2, const arma::vec& z, const arma::mat& w, arma::mat& B_hat, arma::vec& mu, arma::mat& Sigma,
	 	  double& sigmab, arma::mat& Ve, arma::vec& alpha, double& sigmaz,
	 	  arma::vec& LogLik, const bool& constr, const bool& PX, int& Iteration, const int& maxIter);

List	mammot_PXem_ss(arma::mat x1, arma::mat y, arma::mat x2, arma::vec z, arma::mat w, 
					bool constr, bool PX, int maxIter);


#endif /* PXem_ss_hpp */