#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <omp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <armadillo_bits/config.hpp>
#include "plinkfun.hpp"
#include "readExprFile.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List getColNum_Header(std::string filename, char delimiter){
	//cout << "Count columns in the file: ";

	std::ifstream myfile (filename.c_str());
	std::string line, temp;
	
	if (myfile.is_open()){
		getline(myfile, line);
	}

	stringstream ss(line); // Turn the string into a stream.
	string tok;
	CharacterVector header;

	while(getline(ss, tok, delimiter)) {
		header.push_back(tok);
	}

	//cout << "line 1" << line;
	std::istringstream iss(line);
	int columns = 0;
	do{
		std::string sub;
		iss >> sub;
		if (sub.length())
			++columns;
	}
	while(iss);

	//cout << columns << endl;
	
	List output = List::create(//Rcpp::Named("NumGene") = NumGene,
		Rcpp::Named("columns") = columns,
		Rcpp::Named("header") = header);

	return output;
}

CharacterVector charv_subset2(CharacterVector x, uvec idx){
	CharacterVector v(idx.n_elem);
	for (unsigned int i = 0; i < idx.n_elem; i++){
		v(i) = x(idx(i));
	}
	return v;
}

// format of covariate file is "FID IID, covar1, covar2, covar3 ..."
Rcpp::List getCovarFile(std::string filename, char delimiter, int ncols, int nrows){
  std::ifstream myfile (filename.c_str());
  mat covar(nrows, ncols-2);
  std::string line;
  CharacterVector IID(nrows), FID(nrows);

  //clock_t t1 = clock();

  int nrow_ind = 0;
  vector <string> tmp;

  if (myfile.is_open()){
	  while (nrow_ind < nrows ){

		  getline(myfile, line);

		  boost::split(tmp, line, boost::is_any_of(" \t *"));

		  IID(nrow_ind) = tmp[1];
		  FID(nrow_ind) = tmp[0];

		  //   cout << ncols-ncols_omit << ";" << nrow_ind << endl;


		  for (int j = 0; j < ncols - 2; j++){
			  if (tmp[j + 2].compare("NA") == 0){
				  covar(nrow_ind, j) = log(-10);
			  }
			  else {
				  covar(nrow_ind, j) = atof(tmp[j + 2].c_str());
			  }
		  }

		  nrow_ind++;

	  }
  }

  List output = List::create(Rcpp::Named("IID") = IID,
	  Rcpp::Named("FID") = FID,
      Rcpp::Named("covar") = covar);
  return output;
}

// [[Rcpp::export]]
Rcpp::List dataLoader(std::string stringname1, std::string stringname2, std::vector<std::string> stringname3, 
	std::string stringname4, std::string stringname5, int whCol){
	//match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
	// plink file 1: stringname1; plink file 2: stringname2; expression file: stringname3
	// covariates file for file 1: stringname4; covariates file for file 2: stringname5
	// NOTE:
	//      1. overlapped individuals in expression files should be a subset of plink file1,
	// and also should be a subset of covariate file1.
	//		2. individuals in (row order of ) covariate file2 should be exactly the same as 
	// those in the fam file of GWAS (plink file 2).
	// whCol: which pheno is used (1 stands for first pheno (6th column of fam file and so on )
	cout << "## Start matching SNPs in plink files (1, 2). " << endl;

	List tmp =  match_SNPs(stringname1, stringname2);
	uvec idxinFile1 = tmp["idxinFile1"];
	uvec idxinFile2 = tmp["idxinFile2"];
	vec ind = tmp["indicator"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	uvec A1_1_r = tmp["A1_1_r"], A1_2_r = tmp["A1_2_r"], A2_1_r = tmp["A2_1_r"], A2_2_r = tmp["A2_2_r"];

	// load IID in file 1 (fam file 1)
	cout << "## Start loading fam files (1, 2). " << endl;
	string famfile1 = stringname1;
	famfile1 += ".fam";
	int N1 = getLineNum(famfile1);

	IntegerVector tmp1(N1);
	CharacterVector FID_1(N1), IID_1(N1);
	NumericVector tmp2(N1);
	ReadPlinkFamFile(famfile1, FID_1, IID_1, tmp1, tmp2, N1);

	// load pheno in file 2 (fam file 2)
	string famfile2 = stringname2;
	famfile2 += ".fam";
	int N2 = getLineNum(famfile2);

	IntegerVector sex_2(N2);
	NumericVector pheno_2(N2);
	CharacterVector FID_2(N2), IID_2(N2);
	// ReadPlinkFamFile(famfile2, FID_2, IID_2, sex_2, pheno_2, N2);
	ReadPlinkFamFile2(famfile2, FID_2, IID_2, pheno_2, N2, whCol);

	vec y = as<vec>(pheno_2);




// read multiple expression files.
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
	// get intersections of sample IDs, gene IDs for all expression files and
	// TWAS files (plink file 1).
	cout << "## Start loading multiple expression files. " << endl;

	CharacterVector indiv_inter = IID_1;
	char delimiter = '\t';   // for all the expression files.
	int t = 0, T = stringname3.size();

	List IID_expr(T);
	ivec Nindiv = ones<ivec>(T);
	ivec Ngene  = ones<ivec>(T);
	for (t = 0; t < T; t++)	{	// overlapped individuals.
		List tmp_e = getColNum_Header(stringname3[t], delimiter);

		Nindiv(t) = tmp_e["columns"];
		Ngene(t) = getLineNum(stringname3[t]) - 1; // header is first line;
		CharacterVector header = tmp_e["header"];
		uvec idxtmp = linspace<uvec>(6, Nindiv(t)-1, Nindiv(t) - 6);
		CharacterVector IID_expr_tmp = charv_subset2(header, idxtmp);
		IID_expr[t] = IID_expr_tmp;

		indiv_inter = intersect(indiv_inter, header);
	}
//cout << " indiv_inter: " << indiv_inter.size() << endl;


	// read expression files and 
 	List expr_tmp(T);
 	vec  lower, upper, chr_expr;
 	CharacterVector genetype1, genetype2;
 	CharacterVector gene_inter;
 	for (t = 0; t < T; t++)	{	 
 		List expr_t = getExprFile2(stringname3[t], delimiter, Nindiv(t), Ngene(t) + 1, 6);
 		expr_tmp[t] = expr_t;

 		// find overlapped genes.
 		CharacterVector targetID = expr_t["targetID"];
 		if (t == 0) gene_inter = targetID;
 		else gene_inter = intersect(gene_inter, targetID);

 		// save information of overlapped genes.
 		if (t == (T-1))	{
 			IntegerVector ind_gene = match(gene_inter, targetID) - 1;
 			uvec index_gene = as<uvec>(ind_gene);
 			CharacterVector genetype1_t = expr_t["genetype1"];
			CharacterVector genetype2_t = expr_t["genetype2"];
			vec chr_expr_t = expr_t["chr_expr"];
 			vec lower_t = expr_t["lower"];
			vec upper_t = expr_t["upper"];

			genetype1 = genetype1_t[ind_gene];
			genetype2 = genetype2_t[ind_gene];
			chr_expr = chr_expr_t(index_gene);
 			lower = lower_t(index_gene);
 			upper = upper_t(index_gene);
 		}
 	} 

 	// draw overlapped samples and genes acorss all expression files.
 	cube expr_used(gene_inter.size(), indiv_inter.size(), T, fill::zeros);
 	for (t =0; t < T; t++)	{
 		List expr_t = expr_tmp[t];

 		CharacterVector IID_expr_tmp = IID_expr[t];
 		IntegerVector ind_indiv = match(indiv_inter, IID_expr_tmp) - 1;
 		uvec index_indiv = as<uvec>(ind_indiv);

        CharacterVector targetID = expr_t["targetID"];
 		IntegerVector ind_gene = match(gene_inter, targetID) - 1;
 		uvec index_gene = as<uvec>(ind_gene);

 		mat expr_all = expr_t["expr"];
 		expr_used.slice(t) = expr_all(index_gene, index_indiv);	

 	}


// read two plink files.
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
 	clock_t t1 = clock();
	cout << "## Start loading genotype files 1, " ;
	//read file 1 (plink bim bed fam files)
	string bimfile1 = stringname1;
	bimfile1 += ".bim";
	long int P1 =  getLineNum(bimfile1);
	cout << N1 << "-by-" << P1;
	arma::Mat<unsigned> X1(N1,P1);
	//cout << " break 1: " << N1 << "X" << P1 << endl;
	readPlink(stringname1, N1, P1, X1.memptr());

	// subset overlapped samples and SNPs.
	IntegerVector indiv_idx_in1 = match(indiv_inter, IID_1) -1;
	X1 = X1.cols(idxinFile1 - 1);
	X1 = X1.rows(as<uvec>(indiv_idx_in1));
	// replace NA with 0 dosage
	X1.replace(3, 0);

	cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;



	t1 = clock();
	cout << "## Start loading genotype files 2, ";
	//read file 2 (plink bim bed fam files)
	//cout << endl;
	string bimfile2 = stringname2;
	bimfile2 += ".bim";
	//cout << "break 1... " << endl;
	long int P2 =  getLineNum(bimfile2);
	cout << N2 << "-by-" << P2;
	arma::Mat<unsigned> X2(N2,P2);
	//cout << "break 3..." << endl;
	readPlink(stringname2, N2, P2, X2.memptr());

	// subset overlapped SNPs
	X2 = X2.cols(idxinFile2 -1);
	// replace NA with 0 dosage
	X2.replace(3, 0);

	// if needed, change the reference allele in file 2
	uvec idxtmp = linspace<uvec>(1, idxinFile1.size(), idxinFile1.size());
	uvec idx = find(ind == -1);
	uvec idxc = idxtmp.elem(idx);

	X2.cols(idxc-1) = 2 - X2.cols(idxc-1);
 	cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;



	cout << "## Start loading covariates files ... " << endl;
	// load covariates file w1
	mat covar1;
	//cout << "break ..." << covar2.n_rows << ";" << covar2.n_cols << ";" << N1 << ";" << !stringname4.empty() << endl;
	CharacterVector IID_w1;
	if (!stringname4.empty()){
		int Nincovar1 = getLineNum(stringname4);
		tmp = getColNum_Header(stringname4, delimiter);
		int Ncovar = tmp["columns"];
		tmp = getCovarFile(stringname4, delimiter, Ncovar, Nincovar1);
		mat covartmp = tmp["covar"];
		
		IID_w1 = tmp["IID"];

		mat w1one = ones<mat>(X1.n_rows, 1);
		IntegerVector  indiv_idx_covar1 = match(indiv_inter, IID_w1) - 1;

		covar1 = join_rows(w1one, covartmp.rows(as<uvec>(indiv_idx_covar1)));

	}
	else {
		covar1 = ones<mat>(X1.n_rows, 1);
	}
	
	// load covariates file w2
	mat covar2;
	if (!stringname5.empty()){
		tmp = getColNum_Header(stringname5, delimiter);
		int Ncovar = tmp["columns"];
		tmp = getCovarFile(stringname5, delimiter, Ncovar, N2);
		mat covartmp = tmp["covar"];
		
		mat w2one = ones<mat>(N2, 1);
		//cout << endl << "break ..." << covartmp.n_rows << ";" << covartmp.n_cols << ";" << N2 << endl;
		covar2 = join_rows(w2one, covartmp);
	}
	else {
		covar2 = ones<mat>(N2, 1);
	}

	cout << "## End loading files ... " << endl;


 	List output = List::create(
 							   Rcpp::Named("covar1") = covar1,
 							   Rcpp::Named("X1") = X1,
 							   Rcpp::Named("rsname_4use_r") = rsname_4use_r,
 							   Rcpp::Named("chr_4use_r") = chr_4use_r,
							   Rcpp::Named("bp_4use_r") = bp_4use_r,
							   
							   Rcpp::Named("lower") = lower,
							   Rcpp::Named("upper") = upper,
							   Rcpp::Named("genetype1") = genetype1,
							   Rcpp::Named("genetype2") = genetype2,
							   Rcpp::Named("targetID") = gene_inter,
							   Rcpp::Named("chr_expr") = chr_expr,					   	
 							   Rcpp::Named("indiv_inter") = indiv_inter,
							   Rcpp::Named("expr_used") = expr_used,

							   Rcpp::Named("covar2") = covar2,							   
 							   Rcpp::Named("X2") = X2,
							   Rcpp::Named("y") = y);

	return output;
}




// [[Rcpp::export]]
Rcpp::List dataLoaderSS(std::string stringname1, std::string stringname2, std::vector<std::string> stringname3, 
	std::string stringname4, std::string stringnameSS){
	//match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
	// plink file 1: stringname1; plink file 2: stringname2; expression file: stringname3
	// covariates file for file 1: stringname4; covariates file for file 2: stringname5
	// NOTE:
	//      1. overlapped individuals in expression files should be a subset of plink file1,
	// and also should be a subset of covariate file1.
	//		2. individuals in (row order of ) covariate file2 should be exactly the same as 
	// those in the fam file of GWAS (plink file 2).
	// whCol: which pheno is used (1 stands for first pheno (6th column of fam file and so on )
	cout << "## Start matching SNPs in plink files (1, 2). " << endl;

	List tmp =  match_SNPs_SS(stringname1, stringname2, stringnameSS);
	uvec idxinFile1 = tmp["idxinFile1"];
	uvec idxinFile2 = tmp["idxinFile2"];
	uvec idxinFileSS = tmp["idxinFileSS"];
	vec ind = tmp["indicator"];
	vec indSS = tmp["indicatorSS"];

	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	uvec A1_1_r = tmp["A1_1_r"], A1_2_r = tmp["A1_2_r"], A2_1_r = tmp["A2_1_r"], A2_2_r = tmp["A2_2_r"];

	// load IID in file 1 (fam file 1)
	cout << "## Start loading fam files (1, 2). " << endl;
	string famfile1 = stringname1;
	famfile1 += ".fam";
	int N1 = getLineNum(famfile1);

	IntegerVector tmp1(N1);
	CharacterVector FID_1(N1), IID_1(N1);
	NumericVector tmp2(N1);
	ReadPlinkFamFile(famfile1, FID_1, IID_1, tmp1, tmp2, N1);

	

// read multiple expression files.
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
	// get intersections of sample IDs, gene IDs for all expression files and
	// TWAS files (plink file 1).
	cout << "## Start loading multiple expression files. " << endl;

	CharacterVector indiv_inter = IID_1;
	char delimiter = '\t';   // for all the expression files.
	int t = 0, T = stringname3.size();

	List IID_expr(T);
	ivec Nindiv = ones<ivec>(T);
	ivec Ngene  = ones<ivec>(T);
	for (t = 0; t < T; t++)	{	// overlapped individuals.
		List tmp_e = getColNum_Header(stringname3[t], delimiter);

		Nindiv(t) = tmp_e["columns"];
		Ngene(t) = getLineNum(stringname3[t]) - 1; // header is first line;
		CharacterVector header = tmp_e["header"];
		uvec idxtmp = linspace<uvec>(6, Nindiv(t)-1, Nindiv(t) - 6);
		CharacterVector IID_expr_tmp = charv_subset2(header, idxtmp);
		IID_expr[t] = IID_expr_tmp;

		indiv_inter = intersect(indiv_inter, header);
	}
//cout << " indiv_inter: " << indiv_inter.size() << endl;


	// read expression files and 
 	List expr_tmp(T);
 	vec  lower, upper, chr_expr;
 	CharacterVector genetype1, genetype2;
 	CharacterVector gene_inter;
 	for (t = 0; t < T; t++)	{	 
 		List expr_t = getExprFile2(stringname3[t], delimiter, Nindiv(t), Ngene(t) + 1, 6);
 		expr_tmp[t] = expr_t;

 		// find overlapped genes.
 		CharacterVector targetID = expr_t["targetID"];
 		if (t == 0) gene_inter = targetID;
 		else gene_inter = intersect(gene_inter, targetID);

 		// save information of overlapped genes.
 		if (t == (T-1))	{
 			IntegerVector ind_gene = match(gene_inter, targetID) - 1;
 			uvec index_gene = as<uvec>(ind_gene);
 			CharacterVector genetype1_t = expr_t["genetype1"];
			CharacterVector genetype2_t = expr_t["genetype2"];
			vec chr_expr_t = expr_t["chr_expr"];
 			vec lower_t = expr_t["lower"];
			vec upper_t = expr_t["upper"];

			genetype1 = genetype1_t[ind_gene];
			genetype2 = genetype2_t[ind_gene];
			chr_expr = chr_expr_t(index_gene);
 			lower = lower_t(index_gene);
 			upper = upper_t(index_gene);
 		}
 	} 
    cout << "just completed" << endl;
 	// draw overlapped samples and genes acorss all expression files.
 	cube expr_used(gene_inter.size(), indiv_inter.size(), T, fill::zeros);
 	for (t =0; t < T; t++)	{
 		List expr_t = expr_tmp[t];

 		CharacterVector IID_expr_tmp = IID_expr[t];
 		IntegerVector ind_indiv = match(indiv_inter, IID_expr_tmp) - 1;
 		uvec index_indiv = as<uvec>(ind_indiv);

        CharacterVector targetID = expr_t["targetID"];
 		IntegerVector ind_gene = match(gene_inter, targetID) - 1;
 		uvec index_gene = as<uvec>(ind_gene);

 		mat expr_all = expr_t["expr"];
 		expr_used.slice(t) = expr_all(index_gene, index_indiv);	

 	}


// read two plink files.
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
 	clock_t t1 = clock();
	cout << "## Start loading genotype files 1, " ;
	//read file 1 (plink bim bed fam files)
	string bimfile1 = stringname1;
	bimfile1 += ".bim";
	long int P1 =  getLineNum(bimfile1);
	cout << N1 << "-by-" << P1;
	arma::Mat<unsigned> X1(N1,P1);
	//cout << " break 1: " << N1 << "X" << P1 << endl;
	readPlink(stringname1, N1, P1, X1.memptr());

	// subset overlapped samples and SNPs.
	IntegerVector indiv_idx_in1 = match(indiv_inter, IID_1) -1;
	X1 = X1.cols(idxinFile1 - 1);
	X1 = X1.rows(as<uvec>(indiv_idx_in1));
	// replace NA with 0 dosage
	X1.replace(3, 0);

	cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;



	t1 = clock();
	cout << "## Start loading genotype files of reference panel " << endl;
	// load pheno in file 2 (fam file 2)
	string famfile2 = stringname2;
	famfile2 += ".fam";
	int N2 = getLineNum(famfile2);
	//cout << endl;
	string bimfile2 = stringname2;
	bimfile2 += ".bim";
	//cout << "break 1... " << endl;
	long int P2 =  getLineNum(bimfile2);
	cout << N2 << "-by-" << P2;
	arma::Mat<unsigned> X2(N2,P2);
	//cout << "break 3..." << endl;
	readPlink(stringname2, N2, P2, X2.memptr());

	// subset overlapped SNPs
	X2 = X2.cols(idxinFile2 -1);
	// replace NA with 0 dosage
	X2.replace(3, 0);

	// if needed, change the reference allele in file 2
	uvec idxtmp = linspace<uvec>(1, idxinFile1.size(), idxinFile1.size());
	uvec idx = find(ind == -1);
	uvec idxc = idxtmp.elem(idx);

	X2.cols(idxc-1) = 2 - X2.cols(idxc-1);
 	cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;



	cout << "## Start loading covariates files ... " << endl;
	// load covariates file w1
	mat covar1;
	//cout << "break ..." << covar2.n_rows << ";" << covar2.n_cols << ";" << N1 << ";" << !stringname4.empty() << endl;
	CharacterVector IID_w1;
	if (!stringname4.empty()){
		int Nincovar1 = getLineNum(stringname4);
		tmp = getColNum_Header(stringname4, delimiter);
		int Ncovar = tmp["columns"];
		tmp = getCovarFile(stringname4, delimiter, Ncovar, Nincovar1);
		mat covartmp = tmp["covar"];
		
		IID_w1 = tmp["IID"];

		mat w1one = ones<mat>(X1.n_rows, 1);
		IntegerVector  indiv_idx_covar1 = match(indiv_inter, IID_w1) - 1;

		covar1 = join_rows(w1one, covartmp.rows(as<uvec>(indiv_idx_covar1)));

	}
	else {
		covar1 = ones<mat>(X1.n_rows, 1);
	}
	
	// load Summary Statisitc files.
	int P3 = getLineNum(stringnameSS);
	List GWAS = read_GWAS(stringnameSS, P3);
	mat GWAS_sum = GWAS["GWAS_sum"];
	GWAS_sum = GWAS_sum.rows(idxinFileSS - 1);

	// if needed, the signs of marginal coefficients in Summary Statistics file are changed.
	GWAS_sum.cols(zeros<uvec>(1)) = GWAS_sum.cols(zeros<uvec>(1)) % indSS;

	cout << "## End loading files ... " << endl;


 	List output = List::create(
 							   Rcpp::Named("covar1") = covar1,
 							   Rcpp::Named("X1") = X1,
 							   Rcpp::Named("rsname_4use_r") = rsname_4use_r,
 							   Rcpp::Named("chr_4use_r") = chr_4use_r,
							   Rcpp::Named("bp_4use_r") = bp_4use_r,
							   
							   Rcpp::Named("lower") = lower,
							   Rcpp::Named("upper") = upper,
							   Rcpp::Named("genetype1") = genetype1,
							   Rcpp::Named("genetype2") = genetype2,
							   Rcpp::Named("targetID") = gene_inter,
							   Rcpp::Named("chr_expr") = chr_expr,					   	
 							   Rcpp::Named("indiv_inter") = indiv_inter,
							   Rcpp::Named("expr_used") = expr_used,	   
 							   Rcpp::Named("Xpanel") = X2,
 							   Rcpp::Named("GWAS_SS") = GWAS_sum);

	return output;
}