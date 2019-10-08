
#ifndef paral_hpp
#define paral_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include "PXem.hpp"
#include "PXem_ss.hpp"
#include "PXem_part.hpp"
 
using namespace std;
using namespace arma;

class parGene{
  public:
    int current_idx=0;

    // std::mutex mtx;

    
    Mat<unsigned> X1, X2;
    cube expr;
    vec z;
    mat Prjy;
    mat out_param, Alpha, expr_info, snp_info;
    uword Ngene;
    int bw;

    parGene(const Mat<unsigned>& X1, const Mat<unsigned>& X2, const cube& expr, const vec& z, const mat& Prjy, 
		    mat& out_param, mat& Alpha, const mat expr_info, const mat snp_info, const uword Ngene, const int bw){
        this -> X1 = X1;
        this -> X2 = X2;
        this -> expr = expr;
		this -> z = z;
        this -> Prjy = Prjy;
        this -> out_param = out_param;
        this -> Alpha = Alpha;
		this -> expr_info = expr_info;
		this -> snp_info = snp_info;
		this -> Ngene = Ngene;
        this -> bw = bw;
    }

    void loop_by_gene(int i);
    void update_by_thread(int thread_id);
    int  next();

};

class parGeneSS{
  public:
    int current_idx=0;

    Mat<unsigned> X1, Xpanel;
    cube expr;
    mat GWAS_SS, Prjy;
    mat out_param, Alpha, expr_info, snp_info;
    uword Ngene;
    double lam;
    int bw;

    parGeneSS(const Mat<unsigned>& X1, const Mat<unsigned>& Xpanel, const cube& expr, const mat& GWAS_SS, const mat& Prjy, 
            mat& out_param, mat& Alpha, const mat expr_info, const mat snp_info, const uword Ngene, const double lam, const int bw){
        this -> X1 = X1;
        this -> Xpanel = Xpanel;
        this -> expr = expr;
        this -> GWAS_SS = GWAS_SS;
        this -> Prjy = Prjy;
        this -> out_param = out_param;
        this -> Alpha = Alpha;
        this -> expr_info = expr_info;
        this -> snp_info = snp_info;
        this -> Ngene = Ngene;
        this -> lam = lam;
        this -> bw = bw;   
    }

    void loop_by_gene(int i);
    void update_by_thread(int thread_id);
    int  next();

};


class parGene_part{
  public:
    int current_idx=0;

    Mat<unsigned> X1, X2;
    cube expr;
    vec z;
    mat Prjy;
    mat out_param, Alpha, chisq, h_y2, expr_info, snp_info;
    uword Ngene;
    int bw;

    parGene_part(const Mat<unsigned>& X1, const Mat<unsigned>& X2, const cube& expr, const vec& z, const mat& Prjy, 
            mat& out_param, mat& Alpha, mat& chisq, mat& h_y2, const mat expr_info, const mat snp_info, const uword Ngene, const int bw){
        this -> X1 = X1;
        this -> X2 = X2;
        this -> expr = expr;
        this -> z = z;
        this -> Prjy = Prjy;
        this -> out_param = out_param;
        this -> Alpha = Alpha;
        this -> chisq = chisq;
        this -> h_y2 = h_y2;
        this -> expr_info = expr_info;
        this -> snp_info = snp_info;
        this -> Ngene = Ngene;
        this -> bw = bw;
    }

    void loop_by_gene(int i);
    void update_by_thread(int thread_id);
    int  next();

};

class parGeneSS_part{
  public:
    int current_idx=0;

    Mat<unsigned> X1, Xpanel;
    cube expr;
    mat GWAS_SS, Prjy;
    mat out_param, Alpha, chisq, h_y2, expr_info, snp_info;
    uword Ngene;
    double lam;
    int bw;

    parGeneSS_part(const Mat<unsigned>& X1, const Mat<unsigned>& Xpanel, const cube& expr, const mat& GWAS_SS, const mat& Prjy, 
            mat& out_param, mat& Alpha, mat& chisq, mat& h_y2, const mat expr_info, const mat snp_info, const uword Ngene, const double lam, const int bw){
        this -> X1 = X1;
        this -> Xpanel = Xpanel;
        this -> expr = expr;
        this -> GWAS_SS = GWAS_SS;
        this -> Prjy = Prjy;
        this -> out_param = out_param;
        this -> Alpha = Alpha;
        this -> chisq = chisq;
        this -> h_y2 = h_y2;
        this -> expr_info = expr_info;
        this -> snp_info = snp_info;
        this -> Ngene = Ngene;
        this -> lam = lam;
        this -> bw = bw;   
    }

    void loop_by_gene(int i);
    void update_by_thread(int thread_id);
    int  next();

};

#endif /* paral_hpp */
