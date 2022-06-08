#include <RcppArmadillo.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <string>
#include <vector>


#include "data_loader.hpp"
#include "paral.hpp"
 
using namespace std;
using namespace Rcpp;
using namespace arma;

// There are four functions here:
// the first two functions conduct MAMMOT based on individual eQTL data and GWAS individual genotypes.
// 1. mammot
// 2. mammot_paral
// the above functions have exactly the same input and output, except that mammot_paral has an extra
// argument "coreNum", which specify the number of cores you may want to use.
//
// the last two functions conduct MAMMOT based on individual eQTL data and GWAS summary statisitcs based
// on margianl linear regression analysis:
// 3. mammotSS
// 4. mammotSS_paral
// Note that we assume the resulting coefficients are not standardized,
// which is the default settings in plink, (output of "--assoc" or "--linear"), 



//' @title
//' CoMM
//' @description
//' CoMM to dissecting genetic contributions to complex traits by leveraging regulatory information.
//'
//' @param stringname1  prefix for eQTL genotype file with plink format (bim, bed).
//' @param stringname2  prefix for GWAS genotype and phenotype file with plink format (bim, bed, fam).
//' @param stringname3  gene expression file with full name.
//' @param stringname4  covariates file for eQTL data.
//' @param stringname5  covariates file for GWAS data, e.g. top 10 PCs.
//' @param whCol  specify which phenotype is used in fam. For example, when whCol = 2, the seven-th column of fam file will be used as phenotype.
//' @param bw  the number of downstream and upstream SNPs that are considered as cis-SNP within a gene.
//'
//' @return List of model parameters
//'
//' @examples
//' ##Working with no summary statistics, no covariates and options
//' file1 = "1000G.EUR.QC.1";
//' file2 = "NFBC_filter_mph10";
//' file3 = "Geuvadis_gene_expression_qn.txt";
//' file4 = "";
//' file5 = "pc5_NFBC_filter_mph10.txt";
//' whichPheno = 1;
//' bw = 500000;
//'
//' fm = CoMM_testing_run(file1,file2,file3, file4,file5, whichPheno, bw);
//'
//' @details
//' \code{CoMM} fits the CoMM model. It requires to provide plink binary eQTL genotype file (bim, bed)
//' the GWAS plink binary file (bim, bed, fam), gene expression file for eQTL.
//' @export
// [[Rcpp::export]]
Rcpp::List mammot(std::string stringname1, std::string stringname2, 
							std::vector<std::string> stringname3, std::string stringname4,
							std::string stringname5, int whCol, int bw){ 
							//int normalize_option = 1, int pred_option = 0){
							//, char* A21, char* A22){
	// normalize_option: 1. normalize each separately, 2. normalize both plink files together
	// match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
	// plink file 1: stringname1; plink file 2: stringname2; expression file: stringname3
	// covariates file for file 1: stringname4; covariates file for file 2: stringname5
	// pred_option :0 (no calculation for prediction) 1 (calcuation for prediction)

	List tmp = dataLoader(stringname1, stringname2, stringname3, stringname4, stringname5, whCol);
	Mat<unsigned> X1 = tmp["X1"], X2 = tmp["X2"];
	vec z = tmp["y"];
	mat w1 = tmp["covar1"], w2 = tmp["covar2"];
	cube expr = tmp["expr_used"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], 
	targetID = tmp["targetID"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];

	uword Ngene = lower.size();
	//Ngene = 30;
	int T = stringname3.size();
	umat snp_info = zeros<umat>(chr_4use_r.n_elem, 2);
	snp_info.col(0) = chr_4use_r;
	snp_info.col(1) = bp_4use_r;

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	uvec idx;
	mat out_param = datum::nan * ones<mat>(Ngene, 4);
	mat Alpha(Ngene, T); Alpha.fill(datum::nan);
	uvec idx_all = zeros<uvec>(0);

	// pre-compute matrix for intercepts.
	mat Prj1 = w1 * inv_sympd(w1.t() * w1) * w1.t();
	mat Prj2 = w2 * inv_sympd(w2.t() * w2) * w2.t();
	z = z - Prj2 * z;
	// conduct testing
	List fmHa;

	uword g = 0;
	for (g = 0; g < Ngene; g++){

		idx = find(conv_to<vec>::from(bp_4use_r) < upper(g) + bw && conv_to<vec>::from(bp_4use_r) > lower(g) - bw
			&& conv_to<vec>::from(chr_4use_r) == chr_expr(g));
		out_param(g, 3) = idx.n_elem;

		if ( g % 100 == 0 && g != 0){
					cout << g << "-th Gene starts working ..." << endl ;
					//cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}
		//if (g == 262) break;

		if (idx.is_empty() == false){
			if (idx.n_elem > 1){
				idx_all = join_cols(idx_all, idx);

				mat X1tmp = conv_to<mat>::from((&X1) ->cols(idx));
				mat X2tmp = conv_to<mat>::from((&X2) ->cols(idx));

				mat y = expr.subcube(span(g), span::all, span::all);
				if (y.n_rows == 1) y = y.t();
				y = y - Prj1 * y;

				rowvec meanX1tmp = mean(X1tmp, 0);
				rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
				X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1))/ repmat(sdX1tmp, X1tmp.n_rows, 1) / sqrt(X1tmp.n_cols);


				rowvec meanX2tmp = mean(X2tmp, 0);
				rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
				X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);

				rowvec X1row = X1tmp.row(0);
				X1tmp = X1tmp.cols(find_finite(X1row));
				X2tmp = X2tmp.cols(find_finite(X1row));

				mat w = ones<mat>(X1tmp.n_cols, 1);
				if (T != 1) {
					w = cor(X1tmp, y);
					rowvec sdy = stddev(y, 0, 0); // see manual
					colvec sdx = vectorise(stddev(X1tmp, 0, 0)); 
					w.each_row() %= sdy; 
					w.each_col() /= sdx;	
				}
				// fit model under Ha and H0
				fmHa = mammot_PXem(X1tmp, y, X2tmp, z, w, false, false, 2000);
				List fmHo = mammot_PXem(X1tmp, y, X2tmp, z, w, true, false, 2000);

			
				vec alphatmp = fmHa["alpha"];
				vec loglikHa = fmHa["loglik"], loglikHo = fmHo["loglik"];
				mat Ve = fmHa["Ve"];
				double sigmab = fmHa["sigmab"];
				mat Ve0 = fmHo["Ve"];

				double sigmaz = fmHa["sigmaz"];
				out_param(g, 0) = sigmaz;
				out_param(g, 1) = 2 * (max(loglikHa) - max(loglikHo));
				out_param(g, 2) = sigmab;
				Alpha.row(g) = alphatmp.t(); 

			}
		}

		
	}



	List out = List::create(
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("snp_info_4use_r") = snp_info,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("expr_info") = expr_info,
		Rcpp::Named("out_param") = out_param,
		Rcpp::Named("Alpha") = Alpha);	 
	

	return out;
}






// [[Rcpp::export]]
Rcpp::List mammot_paral(std::string stringname1, std::string stringname2, 
							std::vector<std::string> stringname3, std::string stringname4,
							std::string stringname5, int whCol, int bw, int coreNum){ 
						 
	List tmp = dataLoader(stringname1, stringname2, stringname3, stringname4, stringname5, whCol);
	Mat<unsigned> X1 = tmp["X1"], X2 = tmp["X2"];
	vec z = tmp["y"];
	mat w1 = tmp["covar1"], w2 = tmp["covar2"];
	cube expr = tmp["expr_used"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], 
	targetID = tmp["targetID"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];

	uword Ngene = lower.size();
	//Ngene = 30;
	int T = stringname3.size();
	mat snp_info = zeros<mat>(chr_4use_r.n_elem, 2);
	snp_info.col(0) = conv_to<vec>::from(bp_4use_r);
	snp_info.col(1) = conv_to<vec>::from(chr_4use_r);

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	 
	// pre-compute matrix for intercepts.
	mat Prj2 = w2 * inv_sympd(w2.t() * w2) * w2.t();
	z = z - Prj2 * z;

	mat prj1 = w1 * inv_sympd(w1.t() * w1) * w1.t();
	mat out_param = datum::nan * ones<mat>(4, Ngene);
	mat Alpha = datum::nan * ones<mat>(T, Ngene);
 
	//set parallele structure object
    parGene parObj(X1, X2, expr, z, prj1, out_param, Alpha, expr_info, snp_info, Ngene, bw);

    //set parallel computation
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
        threads[i_thread] = std::thread(&parGene::update_by_thread, &parObj, i_thread);
    }

    for(int i = 0; i < n_thread; i++){
        threads[i].join();
    }
 


	List out = List::create( 
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("snp_info_4use_r") = snp_info,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("expr_info") = expr_info,
		Rcpp::Named("out_param") = trans(parObj.out_param),
		Rcpp::Named("Alpha") = trans(parObj.Alpha));	 

	return out;
}










//' @title
//' mammot
//' @description
//'  
//'
//' @param stringname1  prefix for eQTL genotype file with plink format (bim, bed).
//' @param stringname2  prefix for reference panal GWAS genotype file with plink format (bim, bed, fam).
//' @param stringname3  gene expression file with full name.
//' @param stringname4  covariates file for eQTL data.
//' @param stringname5  GWAS summary statisitcs, which has specific form.
//' @param whCol  specify which phenotype is used in fam. For example, when whCol = 2, the seven-th column of fam file will be used as phenotype.
//' @param bw  the number of downstream and upstream SNPs that are considered as cis-SNP within a gene.
//'
//' @return List of model parameters
//'
//' @examples
//' ##Working with no summary statistics, no covariates and options
//' file1 = "1000G.EUR.QC.1";
//' file2 = "NFBC_filter_mph10";
//' file3 = "Geuvadis_gene_expression_qn.txt";
//' file4 = "";
//' file5 = "pc5_NFBC_filter_mph10.txt";
//' whichPheno = 1;
//' bw = 500000;
//'
//' fm = mammotSS_testing_run(file1,file2,file3, file4,file5, whichPheno, bw);
//'
//' @details
//' \code{mammotSS_testing_run} fits the mammot model. It requires to provide plink binary eQTL genotype file (bim, bed)
//' the GWAS summary statisitcs file (txt), gene expression file for eQTL.
//' @export

// [[Rcpp::export]]
Rcpp::List mammotSS(std::string stringname1, std::string stringname2, 
							std::vector<std::string> stringname3, std::string stringname4,
							std::string stringname5, double lam, int bw){ 
	List tmp = dataLoaderSS(stringname1, stringname2, stringname3, stringname4, stringname5);
	Mat<unsigned> X1 = tmp["X1"], Xpanel = tmp["Xpanel"];
	mat GWAS_SS = tmp["GWAS_SS"];
 	mat w1 = tmp["covar1"];
	cube expr = tmp["expr_used"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], 
	targetID = tmp["targetID"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];

	uword Ngene = lower.size();
	//Ngene = 30;
	int T = stringname3.size();
	umat snp_info = zeros<umat>(chr_4use_r.n_elem, 2);
	snp_info.col(0) = chr_4use_r;
	snp_info.col(1) = bp_4use_r;

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	uvec idx;

	mat out_param = datum::nan * ones<mat>(Ngene, 4);
	mat Alpha(Ngene, T); Alpha.fill(datum::nan);
	uvec idx_all = zeros<uvec>(0);

	// pre-compute matrix for intercepts.
	mat Prj1 = w1 * inv_sympd(w1.t() * w1) * w1.t();

	// conduct testing
	List fmHa;

	uword g = 0;
	for (g = 0; g < Ngene; g++){

		idx = find(conv_to<vec>::from(bp_4use_r) < upper(g) + bw && conv_to<vec>::from(bp_4use_r) > lower(g) - bw
			&& conv_to<vec>::from(chr_4use_r) == chr_expr(g));


		out_param(g, 3) = idx.n_elem;

		if ( g % 100 == 0 && g != 0){
					cout << g << "-th Gene starts working ..." << endl ;
		}

		if (idx.is_empty() == false){
			if (idx.n_elem > 1){
				idx_all = join_cols(idx_all, idx);

				mat X1tmp = conv_to<mat>::from((&X1) ->cols(idx));
				mat Xpaneltmp = conv_to<mat>::from((&Xpanel) ->cols(idx));
				mat GWAS_SStmp = GWAS_SS.rows(idx);

				// deal with One tissue files.
				mat y = expr.subcube(span(g), span::all, span::all);
				if (y.n_rows == 1) y = y.t();
				y = y - Prj1 * y;

				rowvec meanX1tmp = mean(X1tmp, 0);
				rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
				X1tmp = X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1); // just centering

				rowvec meanXpaneltmp = mean(Xpaneltmp, 0);
				Xpaneltmp = Xpaneltmp - repmat(meanXpaneltmp, Xpaneltmp.n_rows, 1); // just centering



				rowvec X1row = 1/sdX1tmp;
				X1tmp = X1tmp.cols(find_finite(X1row));
				Xpaneltmp = Xpaneltmp.cols(find_finite(X1row));
				GWAS_SStmp = GWAS_SStmp.rows(find_finite(X1row));

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

				mat w = ones<mat>(X1tmp.n_cols, 1);
				if (T != 1) {
					w = cor(X1tmp, y);
					rowvec sdy = stddev(y, 0, 0); // see manual
					colvec sdx = vectorise(stddev(X1tmp, 0, 0)); 
					w.each_row() %= sdy; 
					w.each_col() /= sdx;	
				}
				// fit model under Ha and H0
				fmHa = mammot_PXem_ss(X1tmp, y, U, gam_L, w, false, true, 2000);
				List fmHo = mammot_PXem_ss(X1tmp, y, U, gam_L, w, true, true, 2000);

				//cout <<"break 2 ..." <<endl;
 				vec alphatmp = fmHa["alpha"];
				vec loglikHa = fmHa["loglik"], loglikHo = fmHo["loglik"];
				mat Ve = fmHa["Ve"];
				double sigmab = fmHa["sigmab"];
				mat Ve0 = fmHo["Ve"];

				double sigmaz = fmHa["sigmaz"];
				out_param(g, 0) = sigmaz;
				out_param(g, 1) = 2 * (max(loglikHa) - max(loglikHo));
				out_param(g, 2) = sigmab;
				Alpha.row(g) = alphatmp.t(); 
			}			 
		}
	}

	List out = List::create( 
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("snp_info_4use_r") = snp_info,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("expr_info") = expr_info,
		Rcpp::Named("out_param") = out_param,
		Rcpp::Named("Alpha") = Alpha);
		

	return out;
}



 // [[Rcpp::export]]
Rcpp::List mammotSS_paral(std::string stringname1, std::string stringname2, 
							std::vector<std::string> stringname3, std::string stringname4,
							std::string stringname5, double lam, int bw, int coreNum){ 
						 
	List tmp = dataLoaderSS(stringname1, stringname2, stringname3, stringname4, stringname5);
	Mat<unsigned> X1 = tmp["X1"], Xpanel = tmp["Xpanel"];
	mat GWAS_SS = tmp["GWAS_SS"];
	mat w1 = tmp["covar1"];
	cube expr = tmp["expr_used"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], 
	targetID = tmp["targetID"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];

	uword Ngene = lower.size();
	//Ngene = 30;
	int T = stringname3.size();
	mat snp_info = zeros<mat>(chr_4use_r.n_elem, 2);
	snp_info.col(0) = conv_to<vec>::from(bp_4use_r);
	snp_info.col(1) = conv_to<vec>::from(chr_4use_r);

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	 
	// pre-compute matrix for intercepts.
	mat prj1 = w1 * inv_sympd(w1.t() * w1) * w1.t();
	mat out_param = datum::nan * ones<mat>(4, Ngene);
	mat Alpha = datum::nan * ones<mat>(T, Ngene);
 
	//set paralel structure object
    parGeneSS parObjSS(X1, Xpanel, expr, GWAS_SS, prj1, out_param, Alpha, expr_info, snp_info, Ngene, lam, bw);

    //set paralel computation
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
        threads[i_thread] = std::thread(&parGeneSS::update_by_thread, &parObjSS, i_thread);
    }

    for(int i = 0; i < n_thread; i++){
        threads[i].join();
    }
 


	List out = List::create( 
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("snp_info_4use_r") = snp_info,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("expr_info") = expr_info,
		Rcpp::Named("out_param") = trans(parObjSS.out_param),
		Rcpp::Named("Alpha") = trans(parObjSS.Alpha));	 

	return out;
}



// [[Rcpp::export]]
Rcpp::List mammot_part_paral(std::string stringname1, std::string stringname2, 
							std::vector<std::string> stringname3, std::string stringname4,
							std::string stringname5, CharacterVector targetList, int whCol, int bw, int coreNum){ 
	int T = stringname3.size();
	List tmp = dataLoader(stringname1, stringname2, stringname3, stringname4, stringname5, whCol);
	Mat<unsigned> X1 = tmp["X1"], X2 = tmp["X2"];
	vec z = tmp["y"];
	mat w1 = tmp["covar1"], w2 = tmp["covar2"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];

	cube expr = tmp["expr_used"];
	CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], targetID = tmp["targetID"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];

	uword Ngene = targetList.size();
	IntegerVector idx = match(targetList, targetID) - 1;
 	uvec index_gene = as<uvec>(idx);
 	cube expr_idx(Ngene, expr.n_cols, T);
 	for (int j = 0; j < Ngene; j++)
 	{
 		expr_idx.subcube(span(j), span::all, span::all) = expr.subcube(span(index_gene(j)), span::all, span::all);
 	}
	//expr = expr(span(index_gene), span::all, span::all);
	genetype1 = genetype1[idx], genetype2 = genetype2[idx], targetID = targetID[idx];
	lower = lower(index_gene), upper = upper(index_gene), chr_expr = chr_expr(index_gene);

	//Ngene = 2;
	mat snp_info = zeros<mat>(chr_4use_r.n_elem, 2);
	snp_info.col(0) = conv_to<vec>::from(bp_4use_r);
	snp_info.col(1) = conv_to<vec>::from(chr_4use_r);

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	 
	// pre-compute matrix for intercepts.
	mat Prj2 = w2 * inv_sympd(w2.t() * w2) * w2.t();
	z = z - Prj2 * z;

	mat prj1 = w1 * inv_sympd(w1.t() * w1) * w1.t();
	mat out_param = datum::nan * ones<mat>(4, Ngene);
	mat Alpha = datum::nan * ones<mat>(T, Ngene);
 	mat chisq = datum::nan * ones<mat>(T, Ngene);
 	mat h_y2  = datum::nan * ones<mat>(T, Ngene);

	//set paralel structure object
    parGene_part parObj_part(X1, X2, expr_idx, z, prj1, out_param, Alpha, chisq, h_y2, expr_info, snp_info, Ngene, bw);

    //set paralel computation
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
        threads[i_thread] = std::thread(&parGene_part::update_by_thread, &parObj_part, i_thread);
    }

    for(int i = 0; i < n_thread; i++){
        threads[i].join();
    }
 


	List out = List::create( 
		Rcpp::Named("idx") = idx,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("expr_info") = expr_info,
		//Rcpp::Named("out_param") = trans(parObj_part.out_param),
		//Rcpp::Named("Alpha") = trans(parObj_part.Alpha),
		Rcpp::Named("chisq") = trans(parObj_part.chisq),
		Rcpp::Named("h_y2") = trans(parObj_part.h_y2));	 

	return out;
}




 // [[Rcpp::export]]
Rcpp::List mammotSS_part_paral(std::string stringname1, std::string stringname2, 
							std::vector<std::string> stringname3, std::string stringname4,
							std::string stringname5, CharacterVector targetList, double lam, int bw, int coreNum){ 
						 
	List tmp = dataLoaderSS(stringname1, stringname2, stringname3, stringname4, stringname5);
	int T = stringname3.size();
	Mat<unsigned> X1 = tmp["X1"], Xpanel = tmp["Xpanel"];
	mat GWAS_SS = tmp["GWAS_SS"];
	mat w1 = tmp["covar1"];

	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];

	cube expr = tmp["expr_used"];
	CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], targetID = tmp["targetID"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];

	uword Ngene = targetList.size();
	IntegerVector idx = match(targetList, targetID) - 1;
 	uvec index_gene = as<uvec>(idx);
 	cube expr_idx(Ngene, expr.n_cols, T);
 	for (int j = 0; j < Ngene; j++)
 	{
 		expr_idx.subcube(span(j), span::all, span::all) = expr.subcube(span(index_gene(j)), span::all, span::all);
 	}
	//expr = expr(span(index_gene), span::all, span::all);
	genetype1 = genetype1[idx], genetype2 = genetype2[idx], targetID = targetID[idx];
	lower = lower(index_gene), upper = upper(index_gene), chr_expr = chr_expr(index_gene);

	//Ngene = 30;
	mat snp_info = zeros<mat>(chr_4use_r.n_elem, 2);
	snp_info.col(0) = conv_to<vec>::from(bp_4use_r);
	snp_info.col(1) = conv_to<vec>::from(chr_4use_r);

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	 
	// pre-compute matrix for intercepts.
	mat prj1 = w1 * inv_sympd(w1.t() * w1) * w1.t();
	mat out_param = datum::nan * ones<mat>(4, Ngene);
	mat Alpha = datum::nan * ones<mat>(T, Ngene);
 	mat chisq = datum::nan * ones<mat>(T, Ngene);
 	mat h_y2  = datum::nan * ones<mat>(T, Ngene);

	//set paralel structure object
    parGeneSS_part parObjSS_part(X1, Xpanel, expr_idx, GWAS_SS, prj1, out_param, Alpha, chisq, h_y2, expr_info, snp_info, Ngene, lam, bw);

    //set paralel computation
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
        threads[i_thread] = std::thread(&parGeneSS_part::update_by_thread, &parObjSS_part, i_thread);
    }

    for(int i = 0; i < n_thread; i++){
        threads[i].join();
    }
 


	List out = List::create( 
		Rcpp::Named("idx") = idx,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("expr_info") = expr_info,
		//Rcpp::Named("out_param") = trans(parObjSS_part.out_param),
		//Rcpp::Named("Alpha") = trans(parObjSS_part.Alpha),
		Rcpp::Named("chisq") = trans(parObjSS_part.chisq),
		Rcpp::Named("h_y2") = trans(parObjSS_part.h_y2));	 

	return out;
}
