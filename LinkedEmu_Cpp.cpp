// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export("hfunc_linear_Cpp")]]
arma::vec hfunc_linear_Cpp(arma::vec x){ //Linear regression function 
    arma::vec vec_out(2);
    vec_out[0] = 1;
    vec_out[1] = x;
    return(vec_out);
}

// [[Rcpp::export("corrfunc_Gauss_Cpp")]]
double corrfunc_Gauss_Cpp(arma::vec b, arma::vec diff) { //Gaussian correlation function

    double sum1 = 0;
    int nvar = diff.size();
    for (int i=0;i<nvar;i++) {
        double diff_curr = diff[i];
        sum1 += b[i] * pow(diff_curr, 2);
    }
    double res = exp(-sum1);
    return(res);
}

// [[Rcpp::export("corrfunc_Matern32_Cpp")]]
double corrfunc_Matern32_Cpp(arma::vec theta, arma::vec diff) { //Matern correlation function, nu = 3/2

	int nvar = theta.size();
    double res = 1;
    
    for (int i=0;i<nvar;i++) {
    	double theta_curr = theta[i];
        double diff_curr = abs(diff[i]);
        double term1 = (diff_curr * sqrt(3)) / theta_curr;
        double dim_res = (1 + term1) * exp(-term1);
        res *= dim_res;
    }
    return(res);
}

// [[Rcpp::export("corrfunc_Matern52_Cpp")]]
double corrfunc_Matern52_Cpp(arma::vec theta, arma::vec diff) { //Matern correlation function, nu = 5/2

	int nvar = theta.size();
    double res = 1;
    
    for (int i=0;i<nvar;i++) {
    	double theta_curr = theta[i];
        double diff_curr = abs(diff[i]);
        double term1 = (diff_curr * sqrt(5)) / theta_curr;
        double term2 = (5 * pow(diff_curr, 2)) / (3 * pow(theta_curr, 2));
        double dim_res = (1 + term1 + term2) * exp(-term1);
        res *= dim_res;
    }
    return(res);
}

// [[Rcpp::export("corr_pred_Cpp")]]
arma::vec corr_pred_Cpp(arma::vec b, arma::mat xn, arma::vec x, int cov_type){ //Correlation between a point x and the design xn

    int npts = xn.n_rows;
    arma::vec t_vec(npts);
    for (int i=0;i<npts;i++){
        arma::vec xn_row_curr = trans(xn.row(i));
        arma::vec h_vec = x - xn_row_curr;
        if(cov_type==1){
        	t_vec[i] = corrfunc_Gauss_Cpp(b, h_vec);
        }
        else if(cov_type==2){
        	t_vec[i] = corrfunc_Matern32_Cpp(b, h_vec);
        }
        else{
        	t_vec[i] = corrfunc_Matern52_Cpp(b, h_vec);
        }
    }
    return(t_vec);
}

// [[Rcpp::export("nonstdrt")]]
arma::vec nonstdrt(int sampsize, double df, double loc, double scale){ //Sample from univariate non-standardised t distribution
    
    arma::vec norm_vals = arma::randn(sampsize);
    arma::vec chisq_vals = arma::chi2rnd(df, sampsize);
    arma::vec stdt_vals = norm_vals / sqrt(chisq_vals/df);
    arma::vec samp = stdt_vals * scale + loc;
    return(samp);
}

// [[Rcpp::export("GPsamp_multiple")]]
arma::mat GPsamp_multiple(int sampsize, arma::mat predpts, arma::mat xn, arma::vec b, arma::vec yn, arma::mat F_mat, arma::mat corr_mat, int cov_type){

	//Draw a sample from a GP emulator for multiple prediction points

	int npred = predpts.n_rows;
	int ndes = yn.size();
	int dim_reg = F_mat.n_cols;
	arma::mat inv_corr = inv(corr_mat);
	arma::mat inv_mat1 = trans(F_mat) * (inv_corr * F_mat);
	arma::mat mat1 = inv(inv_mat1);
	arma::vec beta_hat = mat1 * (trans(F_mat) * (inv_corr * yn));
	arma::mat mat3 = inv_corr - inv_corr * (F_mat * (mat1 * (trans(F_mat) * inv_corr)));
	double val1 = as_scalar(trans(yn) * (mat3 * yn));
    double nu = ndes - dim_reg;
    double val2 = val1 / nu;
	arma::mat zero_mat = zeros(dim_reg, dim_reg);
	arma::mat mat4_upper = join_rows(zero_mat, trans(F_mat));
	arma::mat mat4_lower = join_rows(F_mat, corr_mat);
	arma::mat mat4 = join_cols(mat4_upper, mat4_lower);
	arma::mat samp(npred, sampsize);
	arma::vec fnew;

	for(int i=0; i < npred; i++){
		arma::vec pred_curr = trans(predpts.row(i));
        if(dim_reg==1){
            fnew = 1;
        }
        else{
            fnew = hfunc_linear_Cpp(pred_curr);
        }
        
        arma::vec rnew = corr_pred_Cpp(b, xn, pred_curr, cov_type);
		arma::vec vec1 = join_cols(fnew, rnew);
		double mu_new = as_scalar(trans(fnew) * beta_hat + trans(rnew) * (inv_corr * (yn - F_mat * beta_hat)));
		double sig2_new = val2 * (1 - as_scalar(trans(vec1) * (inv(mat4) * vec1)));
		double sig_new = sqrt(sig2_new);
		samp.row(i) = trans(nonstdrt(sampsize, nu, mu_new, sig_new));
	}
	
	return(samp);
}

// [[Rcpp::export("condsamps")]]
arma::mat condsamps(arma::vec y1_sample, int y1_loc, arma::mat xn_2, arma::vec b2, arma::vec yn_2, arma::mat F_mat2, arma::mat corr_mat2, arma::vec xnew_2, int cov_type){ 

	//Sample handling for the simulation method

	int sampsize = y1_sample.size();
	int nvar_init = xnew_2.size();
	arma::mat xnew_2a = reshape(xnew_2, 1, nvar_init);
	arma::mat xnew(sampsize, nvar_init+1);
 	if(y1_loc==1){
 		arma::mat xnew_temp = repmat(xnew_2a, sampsize, 1);
 		xnew = join_rows(y1_sample, xnew_temp);
 	}
 	else{
 		arma::mat xnew_temp1 = repmat(xnew_2a.cols(0,y1_loc-2), sampsize, 1);
 		if(nvar_init > y1_loc-1){
 			arma::mat xnew_temp2 = repmat(xnew_2a.cols(y1_loc-1, nvar_init-1), sampsize, 1);
 			arma::mat xnew_c1 = join_rows(xnew_temp1, y1_sample); 
 			xnew = join_rows(xnew_c1, xnew_temp2);
 		}
 		else{
 			xnew = join_rows(xnew_temp1, y1_sample);
 		}
 	}
 	arma::mat y2_samp = GPsamp_multiple(1, xnew, xn_2, b2, yn_2, F_mat2, corr_mat2, cov_type);
 	return(y2_samp);
}

// [[Rcpp::export("PExp_mult_Cpp_full")]]

double PExp_mult_Cpp_full(arma::mat xn, arma::vec a, double beta0, arma::vec theta, double mu1, double sig2_1, arma::vec xnew){

	//Posterior mean of second-model output under linked emulator with inputs other than y1 present in second model

	int m = xn.n_rows;
	int p = xn.n_cols;
	double sum1 = 0;
	for(int i = 0; i < m; i++){
		arma::vec temp1 = trans(xn.row(i));
		arma::vec diff_vec = temp1.subvec(1, p-1) - xnew;
		double corr1 = corrfunc_Gauss_Cpp(theta.subvec(1, p-1), diff_vec);
		double int_curr = exp(- pow(mu1 - xn.at(i, 0), 2) / (2 * sig2_1 + 1 / theta[0])) / sqrt(2 * sig2_1 * theta[0] + 1);
    	sum1 = sum1 + a[i] * corr1 * int_curr;
	}
	return(beta0 + sum1);
}

// [[Rcpp::export("PExp_mult_Cpp_nonewin")]]

double PExp_mult_Cpp_nonewin(arma::mat xn, arma::vec a, double beta0, arma::vec theta, double mu1, double sig2_1){

	//Posterior mean of second-model output under linked emulator when y1 is only input to second model

	int m = xn.n_rows;
	double sum1 = 0;
	for(int i = 0; i < m; i++){
		double int_curr = exp(-pow(mu1 - xn.at(i, 0), 2) / (2 * sig2_1 + 1 / theta[0])) / sqrt(2 * sig2_1 * theta[0] + 1);
    	sum1 = sum1 + a[i] * int_curr;
	}
	return(beta0 + sum1);
}

// [[Rcpp::export("PVar_mult_Cpp_full")]]

double PVar_mult_Cpp_full(arma::mat xn, arma::vec a, arma::vec theta, arma::mat inv_corr, double mu1, double sig2_1, double sig2_2, double nugget, arma::vec xnew){

	//Posterior variance of second-model output under linked emulator with inputs other than y1 present in second model

	int m = xn.n_rows;
	int p = xn.n_cols;
	double sum1 = 0;
	for(int k = 0; k < m; k++){
		for(int l = 0; l < m; l++){
			double term1 = a[k] * a[l] - sig2_2 * inv_corr.at(k, l);
			double sum1a = 0;
			for(int j = 1; j < p; j++){
				double val1 = (pow(xn.at(k, j) - xnew[j-1], 2) + pow(xn.at(l, j) - xnew[j-1], 2)) * theta(j);
    			sum1a = sum1a + val1;
			}
			double corr1 = exp(-sum1a);
			double exp1 = exp(-pow(xn.at(k, 0) - xn.at(l, 0), 2) * theta[0] / 2);
      		double exp2 = exp(-pow((xn.at(k, 0) + xn.at(l, 0))/2 - mu1, 2) / (2 * sig2_1 + 1 / (2 * theta[0])));
			double int_kl = exp1 * exp2 / sqrt(4 * sig2_1 * theta[0] + 1);
      		double res_curr = term1 * corr1 * int_kl;
      		sum1 = sum1 + res_curr;
		}
	}
	double sum2 = 0;
	for(int i = 0; i < m; i++){
		arma::vec temp1 = trans(xn.row(i));
		arma::vec diff_vec = temp1.subvec(1, p-1) - xnew;
		double corr2 = corrfunc_Gauss_Cpp(theta.subvec(1, p-1), diff_vec);
		double int_curr = exp(-pow(mu1 - xn.at(i, 0), 2) / (2 * sig2_1 + 1 / theta[0])) / sqrt(2 * sig2_1 * theta[0] + 1);
    	sum2 = sum2 + a[i] * corr2 * int_curr;
	}
	double var_full = sig2_2 * (1 + nugget) + sum1 - pow(sum2, 2);
	return(var_full);
}

// [[Rcpp::export("PVar_mult_Cpp_nonewin")]]

double PVar_mult_Cpp_nonewin(arma::mat xn, arma::vec a, arma::vec theta, arma::mat inv_corr, double mu1, double sig2_1, double sig2_2, double nugget){

	//Posterior variance of second-model output under linked emulator when y1 is only input to second model

	int m = xn.n_rows;
	double sum1 = 0;
	for(int k = 0; k < m; k++){
		for(int l = 0; l < m; l++){
			double term1 = a[k] * a[l] - sig2_2 * inv_corr.at(k, l);
			double exp1 = exp(-pow(xn.at(k, 0) - xn.at(l, 0), 2) * theta[0] / 2);
      		double exp2 = exp(-pow((xn.at(k, 0) + xn.at(l, 0))/2 - mu1, 2) / (2 * sig2_1 + 1 / (2 * theta[0])));
      		double int_kl = exp1 * exp2 / sqrt(4 * sig2_1 * theta[0] + 1);
      		double res_curr = term1 * int_kl;
      		sum1 = sum1 + res_curr;
		}
	}
	double sum2 = 0;
	for(int i = 0; i < m; i++){
		double int_curr = exp(-pow(mu1 - xn.at(i, 0), 2) / (2 * sig2_1 + 1 / theta[0])) / sqrt(2 * sig2_1 * theta[0] + 1);
    	sum2 = sum2 + a[i] * int_curr;
	}
	double var_full = sig2_2 * (1 + nugget) + sum1 - pow(sum2, 2);
	return(var_full);
}

// [[Rcpp::export("cstar_Cpp")]]
double cstar_Cpp(arma::vec x1, arma::vec x2, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget){ 

    arma::vec xdiff = x1 - x2;
    double cterm = corr_func_Cpp(b, xdiff);
    bool same_x = approx_equal(x1, x2, "absdiff", 1e-7);
    if(same_x==TRUE){
        cterm += nugget;
    }
    arma::vec t1 = corr_pred_Cpp(b, xn, x1, 1);
    arma::vec t2 = corr_pred_Cpp(b, xn, x2, 1);
    arma::vec h1;
    arma::vec h2;
    if(F_mat.n_cols==1){
        h1 = 1;
        h2 = 1;
    }
    else{
        h1 = hfunc_linear_Cpp(x1);
        h2 = hfunc_linear_Cpp(x2);
    }
    double term2 = as_scalar(trans(t1) * (inv_Cmat * t2));
    arma::vec term3a = trans(h1) - trans(t1) * (inv_Cmat * F_mat);
    arma::vec term3b = trans(h2) - trans(t2) * (inv_Cmat * F_mat);
    double term3 = as_scalar(term3a * (W * trans(term3b)));
    double res = cterm - term2 + term3;
    return(res);
}

// [[Rcpp::export("mstar_Cpp")]]
double mstar_Cpp(arma::vec x, arma::vec e, arma::vec b, arma::mat xn, arma::vec beta_hat, int dimh=1){
    arma::vec hterm;
    if(dimh==1){
        hterm = 1;
    }
    else{
        hterm = hfunc_linear_Cpp(x);
    }
    arma::vec tterm = corr_pred_Cpp(b, xn, x, 1);
    double term1 = as_scalar(trans(hterm) * beta_hat);
    double term2 = as_scalar(trans(tterm) * e);
    double res = term1 + term2;
    return(res);
}

//Everything below this is a Monte Carlo integral or multiple integral required for sensitivity analysis, including importance sampling versions for each. Quantities being integrated and the distribution being integrated with respect to are given in the function names.

// [[Rcpp::export("Rint_dG_Cpp_importance")]]
arma::rowvec Rint_dG_Cpp_importance(Function G_func, Function S_func, arma::mat sample, int dimh=1){
    
    int sample_size = sample.n_rows;
    arma::mat int_res(sample_size, dimh);
    arma::vec int_val(dimh);
    arma::vec G_vals_all = Rcpp::as<NumericVector>(G_func(sample));
    arma::vec S_vals_all = Rcpp::as<NumericVector>(S_func(sample));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        double S_val = S_vals_all[i];
        double G_val = G_vals_all[i];
        arma::vec h_val;
        if(dimh==1){
            h_val = 1;
        }
        else{
            h_val = hfunc_linear_Cpp(val_current);
        }
        int_res.row(i) = G_val * h_val / S_val;
    }

    for(int j=0; j < dimh; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Rint_dG_Cpp_exact")]]
arma::rowvec Rint_dG_Cpp_exact(arma::mat sample, int dimh=1){
    
    int sample_size = sample.n_rows;
    arma::mat int_res(sample_size, dimh);
    arma::vec int_val(dimh);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        if(dimh==1){
            int_res.row(i) = 1;
        }
        else{
            int_res.row(i) = hfunc_linear_Cpp(val_current);
        }
    }

    for(int j=0; j < dimh; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Rint_xp_Cpp_importance")]]
arma::rowvec Rint_xp_Cpp_importance(Function G_func, Function S_func, arma::mat sample, arma::vec xp, arma::vec xp_loc, int dimh=1){
    
    int sample_size = sample.n_rows;
    int num_fixed = xp.n_elem;
    arma::mat int_res(sample_size, dimh);
    arma::vec int_val(dimh);
    arma::vec G_vals_all = Rcpp::as<NumericVector>(G_func(sample));
    arma::vec S_vals_all = Rcpp::as<NumericVector>(S_func(sample));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));

        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = xp[j];
            val_current.insert_rows(loc_curr, 1);
            val_current.row(loc_curr) = val_curr;
        }

        double S_val = S_vals_all[i];
        double G_val = G_vals_all[i];
        arma::vec h_val;
        if(dimh==1){
            h_val = 1;
        }
        else{
            h_val = hfunc_linear_Cpp(val_current);
        }
        int_res.row(i) = G_val * h_val / S_val;
    }

    for(int j=0; j < dimh; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Rint_xp_Cpp_exact")]]
arma::rowvec Rint_xp_Cpp_exact(arma::mat sample, arma::vec xp, arma::vec xp_loc, int dimh=1){
    
    int sample_size = sample.n_rows;
    int num_fixed = xp.n_elem;
    arma::mat int_res(sample_size, dimh);
    arma::vec int_val(dimh);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));

        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = xp[j];
            val_current.insert_rows(loc_curr, 1);
            val_current.row(loc_curr) = val_curr;
        }

        if(dimh==1){
            int_res.row(i) = 1;
        }
        else{
            int_res.row(i) = hfunc_linear_Cpp(val_current);
        }
    }

    for(int j=0; j < dimh; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Tint_dG_Cpp_importance")]]
arma::rowvec Tint_dG_Cpp_importance(Function G_func, Function S_func, arma::mat sample, arma::vec b, arma::mat xn){
    
    int sample_size = sample.n_rows;
    int npts = xn.n_rows;
    arma::mat int_res(sample_size, npts);
    arma::vec int_val(npts);
    arma::vec G_vals_all = Rcpp::as<NumericVector>(G_func(sample));
    arma::vec S_vals_all = Rcpp::as<NumericVector>(S_func(sample));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        double S_val = S_vals_all[i];
        double G_val = G_vals_all[i];
        arma::rowvec t_val = trans(corr_pred_Cpp(b, xn, val_current, 1));
        int_res.row(i) = G_val * t_val / S_val;
    }

    for(int j=0; j < npts; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Tint_dG_Cpp_exact")]]
arma::rowvec Tint_dG_Cpp_exact(arma::mat sample, arma::vec b, arma::mat xn){
    
    int sample_size = sample.n_rows;
    int npts = xn.n_rows;
    arma::mat int_res(sample_size, npts);
    arma::vec int_val(npts);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        int_res.row(i) = trans(corr_pred_Cpp(b, xn, val_current, 1));
    }

    for(int j=0; j < npts; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Tint_xp_Cpp_importance")]]
arma::rowvec Tint_xp_Cpp_importance(Function G_func, Function S_func, arma::mat sample, arma::vec b, arma::mat xn, arma::vec xp, arma::vec xp_loc){
    
    int sample_size = sample.n_rows;
    int npts = xn.n_rows;
    int num_fixed = xp.n_elem;
    arma::mat int_res(sample_size, npts);
    arma::vec int_val(npts);
    arma::vec G_vals_all = Rcpp::as<NumericVector>(G_func(sample));
    arma::vec S_vals_all = Rcpp::as<NumericVector>(S_func(sample));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));

        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = xp[j];
            val_current.insert_rows(loc_curr, 1);
            val_current.row(loc_curr) = val_curr;
        }

        double S_val = S_vals_all[i];
        double G_val = G_vals_all[i];
        arma::rowvec t_val = trans(corr_pred_Cpp(b, xn, val_current, 1));
        int_res.row(i) = G_val * t_val / S_val;
    }

    for(int j=0; j < npts; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Tint_xp_Cpp_exact")]]
arma::rowvec Tint_xp_Cpp_exact(arma::mat sample, arma::vec b, arma::mat xn, arma::vec xp, arma::vec xp_loc){
    
    int sample_size = sample.n_rows;
    int npts = xn.n_rows;
    int num_fixed = xp.n_elem;
    arma::mat int_res(sample_size, npts);
    arma::vec int_val(npts);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));

        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = xp[j];
            val_current.insert_rows(loc_curr, 1);
            val_current.row(loc_curr) = val_curr;
        }

        int_res.row(i) = trans(corr_pred_Cpp(b, xn, val_current, 1));
    }

    for(int j=0; j < npts; j++){
        int_val[j] = mean(int_res.col(j));
    }

    arma::rowvec int_final = trans(int_val);
    return(int_final);
}

// [[Rcpp::export("Mstarint_dG_Cpp_importance")]]
double Mstarint_dG_Cpp_importance(Function G_func, Function S_func, arma::mat sample, arma::vec b, arma::mat xn, arma::vec e, arma::vec beta_hat, int dimh=1){
    
    int sample_size = sample.n_rows;
    arma::vec int_vec(sample_size);
    arma::vec G_vals_all = Rcpp::as<NumericVector>(G_func(sample));
    arma::vec S_vals_all = Rcpp::as<NumericVector>(S_func(sample));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        double S_val = S_vals_all[i];
        double G_val = G_vals_all[i];
        double mstar_val = mstar_Cpp(val_current, e, b, xn, beta_hat, dimh);
        int_vec[i] = G_val * pow(mstar_val, 2) / S_val;
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Mstarint_dG_Cpp_exact")]]
double Mstarint_dG_Cpp_exact(arma::mat sample, arma::vec b, arma::mat xn, arma::vec e, arma::vec beta_hat, int dimh=1){
    
    int sample_size = sample.n_rows;
    arma::vec int_vec(sample_size);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        double mstar_val = mstar_Cpp(val_current, e, b, xn, beta_hat, dimh);
        int_vec[i] = pow(mstar_val, 2);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_dG_Cpp_importance")]]
double Cstarint_dG_Cpp_importance(Function G_func, Function S_func, arma::mat sample, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget){
    
    int sample_size = sample.n_rows;
    arma::vec int_vec(sample_size);
    arma::vec G_vals_all = Rcpp::as<NumericVector>(G_func(sample));
    arma::vec S_vals_all = Rcpp::as<NumericVector>(S_func(sample));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        double S_val = S_vals_all[i];
        double G_val = G_vals_all[i];
        double cstar_val = cstar_Cpp(val_current, val_current, inv_Cmat, F_mat, b, xn, W, nugget);
        int_vec[i] = G_val * cstar_val / S_val;
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_dG_Cpp_exact")]]
double Cstarint_dG_Cpp_exact(arma::mat sample, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget){
    
    int sample_size = sample.n_rows;
    arma::vec int_vec(sample_size);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current = trans(sample.row(i));
        int_vec[i] = cstar_Cpp(val_current, val_current, inv_Cmat, F_mat, b, xn, W, nugget);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Mstarint_xp_Cpp_importance")]]
double Mstarint_xp_Cpp_importance(Function Gdp_func, Function Sdp_func, Function Gp_func, Function Sp_func, arma::mat sample1, arma::mat sample2, arma::mat sample3, arma::vec b, arma::mat xn, arma::vec e, arma::vec beta_hat, arma::vec xp_loc, int dimh=1){
    
    int sample_size = sample1.n_rows;
    int num_fixed = sample3.n_cols;
    arma::vec int_vec(sample_size);
    arma::vec G1_vals_all = Rcpp::as<NumericVector>(Gdp_func(sample1));
    arma::vec G2_vals_all = Rcpp::as<NumericVector>(Gdp_func(sample2));
    arma::vec G3_vals_all = Rcpp::as<NumericVector>(Gp_func(sample3));
    arma::vec S1_vals_all = Rcpp::as<NumericVector>(Sdp_func(sample1));
    arma::vec S2_vals_all = Rcpp::as<NumericVector>(Sdp_func(sample2));
    arma::vec S3_vals_all = Rcpp::as<NumericVector>(Sp_func(sample3));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        arma::vec val_current3 = trans(sample3.row(i));
        double S1_val = S1_vals_all[i];
        double S2_val = S2_vals_all[i];
        double S3_val = S3_vals_all[i];
        double G1_val = G1_vals_all[i];
        double G2_val = G2_vals_all[i];
        double G3_val = G3_vals_all[i];
        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = val_current3[j];
            val_current1.insert_rows(loc_curr, 1);
            val_current1.row(loc_curr) = val_curr;
            val_current2.insert_rows(loc_curr, 1);
            val_current2.row(loc_curr) = val_curr;
        }
        
        double mstar_val1 = mstar_Cpp(val_current1, e, b, xn, beta_hat, dimh);
        double mstar_val2 = mstar_Cpp(val_current2, e, b, xn, beta_hat, dimh);
        double prod = mstar_val1 * mstar_val2;
        int_vec[i] = (G1_val * G2_val * G3_val * prod) / (S1_val * S2_val * S3_val);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Mstarint_xp_Cpp_exact")]]
double Mstarint_xp_Cpp_exact(arma::mat sample1, arma::mat sample2, arma::mat sample3, arma::vec b, arma::mat xn, arma::vec e, arma::vec beta_hat, arma::vec xp_loc, int dimh=1){
    
    int sample_size = sample1.n_rows;
    int num_fixed = sample3.n_cols;
    arma::vec int_vec(sample_size);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        arma::vec val_current3 = trans(sample3.row(i));

        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = val_current3[j];
            val_current1.insert_rows(loc_curr, 1);
            val_current1.row(loc_curr) = val_curr;
            val_current2.insert_rows(loc_curr, 1);
            val_current2.row(loc_curr) = val_curr;
        }
        
        double mstar_val1 = mstar_Cpp(val_current1, e, b, xn, beta_hat, dimh);
        double mstar_val2 = mstar_Cpp(val_current2, e, b, xn, beta_hat, dimh);
        int_vec[i] = mstar_val1 * mstar_val2;
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_xp3d_Cpp_importance")]]
double Cstarint_xp3d_Cpp_importance(Function Gdp_func, Function Sdp_func, Function Gp_func, Function Sp_func, arma::mat sample1, arma::mat sample2, arma::mat sample3, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget, arma::vec xp_loc){
    
    int sample_size = sample1.n_rows;
    int num_fixed = sample3.n_cols;
    arma::vec int_vec(sample_size);
    arma::vec G1_vals_all = Rcpp::as<NumericVector>(Gdp_func(sample1));
    arma::vec G2_vals_all = Rcpp::as<NumericVector>(Gdp_func(sample2));
    arma::vec G3_vals_all = Rcpp::as<NumericVector>(Gp_func(sample3));
    arma::vec S1_vals_all = Rcpp::as<NumericVector>(Sdp_func(sample1));
    arma::vec S2_vals_all = Rcpp::as<NumericVector>(Sdp_func(sample2));
    arma::vec S3_vals_all = Rcpp::as<NumericVector>(Sp_func(sample3));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        arma::vec val_current3 = trans(sample3.row(i));
        double S1_val = S1_vals_all[i];
        double S2_val = S2_vals_all[i];
        double S3_val = S3_vals_all[i];
        double G1_val = G1_vals_all[i];
        double G2_val = G2_vals_all[i];
        double G3_val = G3_vals_all[i];
        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = val_current3[j];
            val_current1.insert_rows(loc_curr, 1);
            val_current1.row(loc_curr) = val_curr;
            val_current2.insert_rows(loc_curr, 1);
            val_current2.row(loc_curr) = val_curr;
        }

        double cstar_val = cstar_Cpp(val_current1, val_current2, inv_Cmat, F_mat, b, xn, W, nugget);
        int_vec[i] = (G1_val * G2_val * G3_val * cstar_val) / (S1_val * S2_val * S3_val);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_xp3d_Cpp_exact")]]
double Cstarint_xp3d_Cpp_exact(arma::mat sample1, arma::mat sample2, arma::mat sample3, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget, arma::vec xp_loc){
    
    int sample_size = sample1.n_rows;
    int num_fixed = sample3.n_cols;
    arma::vec int_vec(sample_size);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        arma::vec val_current3 = trans(sample3.row(i));
        for(int j=0; j < num_fixed; j++){
            int loc_curr = xp_loc[j] - 1;
            double val_curr = val_current3[j];
            val_current1.insert_rows(loc_curr, 1);
            val_current1.row(loc_curr) = val_curr;
            val_current2.insert_rows(loc_curr, 1);
            val_current2.row(loc_curr) = val_curr;
        }
        int_vec[i] = cstar_Cpp(val_current1, val_current2, inv_Cmat, F_mat, b, xn, W, nugget);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_xp2d_Cpp_importance")]]
double Cstarint_xp2d_Cpp_importance(Function G_func, Function S_func, arma::mat sample1, arma::mat sample2, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget, double xi, int xi_loc){
    
    int sample_size = sample1.n_rows;
    arma::vec int_vec(sample_size);
    arma::vec G1_vals_all = Rcpp::as<NumericVector>(G_func(sample1));
    arma::vec G2_vals_all = Rcpp::as<NumericVector>(G_func(sample2));
    arma::vec S1_vals_all = Rcpp::as<NumericVector>(S_func(sample1));
    arma::vec S2_vals_all = Rcpp::as<NumericVector>(S_func(sample2));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        double S1_val = S1_vals_all[i];
        double S2_val = S2_vals_all[i];
        double G1_val = G1_vals_all[i];
        double G2_val = G2_vals_all[i];
        val_current1.insert_rows(xi_loc-1, 1);
        val_current1.row(xi_loc-1) = xi;
        val_current2.insert_rows(xi_loc-1, 1);
        val_current2.row(xi_loc-1) = xi;
        double cstar_val = cstar_Cpp(val_current1, val_current2, inv_Cmat, F_mat, b, xn, W, nugget);
        int_vec[i] = (G1_val * G2_val * cstar_val) / (S1_val * S2_val);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_xp2d_Cpp_exact")]]
double Cstarint_xp2d_Cpp_exact(arma::mat sample1, arma::mat sample2, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget, double xi, int xi_loc){
    
    int sample_size = sample1.n_rows;
    arma::vec int_vec(sample_size);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        val_current1.insert_rows(xi_loc-1, 1);
        val_current1.row(xi_loc-1) = xi;
        val_current2.insert_rows(xi_loc-1, 1);
        val_current2.row(xi_loc-1) = xi;
        int_vec[i] = cstar_Cpp(val_current1, val_current2, inv_Cmat, F_mat, b, xn, W, nugget);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_dGdG_Cpp_importance")]]
double Cstarint_dGdG_Cpp_importance(Function G_func, Function S_func, arma::mat sample1, arma::mat sample2, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget){
    
    int sample_size = sample1.n_rows;
    arma::vec int_vec(sample_size);
    arma::vec G1_vals_all = Rcpp::as<NumericVector>(G_func(sample1));
    arma::vec G2_vals_all = Rcpp::as<NumericVector>(G_func(sample2));
    arma::vec S1_vals_all = Rcpp::as<NumericVector>(S_func(sample1));
    arma::vec S2_vals_all = Rcpp::as<NumericVector>(S_func(sample2));
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        double S1_val = S1_vals_all[i];
        double S2_val = S2_vals_all[i];
        double G1_val = G1_vals_all[i];
        double G2_val = G2_vals_all[i];
        double cstar_val = cstar_Cpp(val_current1, val_current2, inv_Cmat, F_mat, b, xn, W, nugget);
        int_vec[i] = (G1_val * G2_val * cstar_val) / (S1_val * S2_val);
    }

    double int_final = mean(int_vec);
    return(int_final);
}

// [[Rcpp::export("Cstarint_dGdG_Cpp_exact")]]
double Cstarint_dGdG_Cpp_exact(arma::mat sample1, arma::mat sample2, arma::mat inv_Cmat, arma::mat F_mat, arma::vec b, arma::mat xn, arma::mat W, double nugget){
    
    int sample_size = sample1.n_rows;
    arma::vec int_vec(sample_size);
    
    for(int i=0; i < sample_size; i++){
        arma::vec val_current1 = trans(sample1.row(i));
        arma::vec val_current2 = trans(sample2.row(i));
        int_vec[i] = cstar_Cpp(val_current1, val_current2, inv_Cmat, F_mat, b, xn, W, nugget);
    }

    double int_final = mean(int_vec);
    return(int_final);
}
