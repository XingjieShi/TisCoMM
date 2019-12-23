# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
using namespace Rcpp;
using namespace arma;


// convert a Additive and dominance components
// back to minor and major alleles
// fixed convert habit: minor is 1, major is 2
// [[Rcpp::export()]]
mat ad2mm(mat x)  {
    int n = x.n_rows;
    int p = x.n_cols;
    mat snp(p, 2*n, fill::zeros);
    for (int i = 0; i < n; i++)  
    {
    	for (int j = 0; j < p; j++)    
        {
            if (x(i, j) == 0) { 
                snp(j, 2*i) = 2;
                snp(j, 2*i+1) = 2;
            }   else if (x(i, j) == 1) { 
                snp(j, 2*i) = 1;
                snp(j, 2*i+1) = 2;
            }   else    {
                snp(j, 2*i) = 1;
                snp(j, 2*i+1) = 1;
            }
        }
    }
    return snp;
} 



// [[Rcpp::export()]]
mat vec2mat(vec V, uword n, uword p) {
    mat M(n, p, fill::zeros);
	for (int j=0; j < p ; j++) 
		M.col(j) = V.subvec(j*n, (j+1)*n-1);
	return M;
}

//convert a matrix of local fdr to a matrix of Global FDR
// [[Rcpp::export()]]
mat fdr2Fdr(mat fdr){
	vec pool =  vectorise(fdr);
    uword M = pool.size();
    uvec indices = sort_index(pool);
    vec sort_pool = pool(indices);
    vec Pool = cumsum(sort_pool) / linspace<vec>(1, M, M);
    Pool.elem(find(Pool  > 1)).ones();
    Pool(indices) = Pool;
    return vec2mat(Pool, fdr.n_rows, fdr.n_cols);
}


//convert a vector of local fdr to a vector of Global FDR
// [[Rcpp::export()]]
vec fdr2FDR(vec fdr){
    uword M = fdr.size();
    uvec indices = sort_index(fdr);
    vec sort_fdr = fdr(indices);
    vec FDR = cumsum(sort_fdr) / linspace<vec>(1, M, M);
    FDR.elem(find(FDR  > 1)).ones();
    FDR(indices) = FDR;
    return FDR;
}

// [[Rcpp::export()]]
double calAUC(vec label, vec pred){
    double auc = 0;
    double m = mean(label);
    vec label2 = label;
    label2(find(label >= m)).fill(1);
    label2(find(label < m)).fill(0);
    label = label2;
    uword N = pred.size();
    uword N_pos = sum(label);
    uvec  idx = sort_index(pred,"descend");
    vec above = linspace<vec>(1, N, N) - cumsum(label(idx));
    auc = (1 - sum(above % label(idx)) / (N_pos * (N-N_pos)));
    auc = auc > 0.5?auc:(1-auc);
    return auc;
}