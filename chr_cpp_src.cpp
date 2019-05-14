#include <RcppArmadillo.h> 
#include "chr_cpp_src.hpp"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
	arma::uword getAnc(arma::uword desc, arma::mat edges) {
		bool flag=0;
		int k=0;
		while (!flag) {
			if (edges(k,1)==desc) {
				flag=1;
				desc=edges(k,0);
			}
			k++;
		}
		return desc;
	}
	
// [[Rcpp::export]]
	uword getMRCA(arma::uword tip1, arma::uword tip2, arma::mat edges) {
		arma::uword anci, ancj;	
		// tips should be zero-indexed
		anci=getAnc(tip1+1, edges);
		ancj=getAnc(tip2+1, edges);
		while (anci != ancj) {
			if (anci > ancj) {
				anci=getAnc(anci, edges);
			} else if (ancj > anci) ancj=getAnc(ancj, edges);
		}
// returns zero-indexed anc
		return anci-1;
	}

// [[Rcpp::export]]
	arma::umat mrca_cpp(arma::mat edges) {	
		const unsigned int ntaxa = edges.col(0).min() - 1;
		arma::uword i, j;
		
		arma::umat mrca(ntaxa, ntaxa);
		mrca.zeros();
		
		for (i=0;i<ntaxa;i++) {
			for (j=i;j<ntaxa;j++) {
				if (i!=j) {
					mrca(i,j)=getMRCA(i, j, edges);
				} else mrca(i,j)=i;
			}
		}
		return arma::symmatu(mrca);
	}

// [[Rcpp::export]]
	arma::mat getMRCADates(arma::mat edges, arma::vec dates) {	
		const unsigned int ntaxa = edges.col(0).min() - 1;
		arma::uword i, j;
		
		arma::umat mrca = mrca_cpp(edges);
		arma::mat this_dates(mrca.n_rows, mrca.n_cols);
		for (i=0;i<mrca.n_cols;i++) this_dates.col(i) = dates.elem(mrca.col(i));	
		return this_dates;
//		return dates.elem(mrca);	
//		return mrca.for_each( (arma::mat::elem_type& val) { return dates(val); }) ;
//		return mrca;
	}

// [[Rcpp::export]]
	arma::mat getVCV_BM(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2) {
		// *** assume shiftDates and sigma2 are sorted oldest to youngest
		// *** assume shiftDates and dates (node dates) are rescaled to time since root
		
		const unsigned int ntaxa = edges.col(0).min() - 1;

		arma::mat mrca_dates(ntaxa, ntaxa);
		arma::vec young_dates(ntaxa);			// either the tip dates or the top of the interval or nothing.
		arma::mat old_dates(ntaxa, ntaxa);
		arma::mat VV(ntaxa, ntaxa, arma::fill::zeros);
		arma::vec VV_diag(ntaxa, arma::fill::zeros);

		mrca_dates = getMRCADates(edges, dates);
		for(uword gamma=1; gamma < shiftDates.n_elem; gamma++) { 

			young_dates = mrca_dates.diag();			
			young_dates.elem(find(young_dates > shiftDates(gamma))).fill(shiftDates(gamma));		// tips are younger than the segment interval, young dates are the top of the segment interval
			young_dates.elem(find(young_dates < shiftDates(gamma-1))).fill(shiftDates(gamma-1));	// tips are older than this segment interval,  young dates are the base of the segment interval
					
			old_dates = mrca_dates;
			old_dates.elem(find(old_dates > shiftDates(gamma))).fill(shiftDates(gamma));			// mrca node is younger this segment interval
			old_dates.elem(find(old_dates < shiftDates(gamma-1))).fill(shiftDates(gamma-1));		// mrca node is older than this segment interval

		// shared time before mrca
			VV += sigma2(gamma-1) * (old_dates - shiftDates(gamma-1));		// shared time is the base of the segment to the mrca node
			
		// vcv diagonal
			VV_diag += sigma2(gamma-1) * (young_dates - shiftDates(gamma-1));		// time is the base of the segment to the mrca node to the tip or top of segment interval
		}		

		VV.diag() = VV_diag;
		return VV;
	}

// [[Rcpp::export]]
	arma::mat getVCV_OU(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2, arma::vec alpha) {
		// *** assume shiftDates, sigma2, and alphas are sorted oldest to youngest
		// *** assume shiftDates and dates (node dates) are rescaled to time since root

		const unsigned int ntaxa = edges.col(0).min() - 1;

		arma::mat mrca_dates(ntaxa, ntaxa);
		arma::vec young_dates(ntaxa);							// either the tip dates or the top of the interval or nothing.
		arma::mat old_dates(ntaxa, ntaxa);
		arma::mat this_shared(ntaxa, ntaxa);
		arma::mat sum_shared(ntaxa, ntaxa, arma::fill::zeros);
		arma::mat sum_indep(ntaxa, ntaxa, arma::fill::zeros);
		arma::mat VV(ntaxa, ntaxa, arma::fill::zeros);
		arma::vec VV_diag(ntaxa, arma::fill::zeros);

		mrca_dates = getMRCADates(edges, dates);
		for(uword gamma=1; gamma < shiftDates.n_elem; gamma++) { 

			young_dates = mrca_dates.diag();			
			young_dates.elem(find(young_dates > shiftDates(gamma))).fill(shiftDates(gamma));		// tips are younger than the segment interval, young dates are the top of the segment interval
			young_dates.elem(find(young_dates < shiftDates(gamma-1))).fill(shiftDates(gamma-1));	// tips are older than this segment interval,  young dates are the base of the segment interval
					
			old_dates = mrca_dates;
			old_dates.elem(find(old_dates > shiftDates(gamma))).fill(shiftDates(gamma));			// mrca node is younger this segment interval
			old_dates.elem(find(old_dates < shiftDates(gamma-1))).fill(shiftDates(gamma-1));		// mrca node is older than this segment interval

		// shared time before mrca
			this_shared = exp(-2 * alpha(gamma-1) * shiftDates(gamma-1)) - exp(-2 * alpha(gamma-1) * old_dates);		// shared time is the base of the segment to the mrca node
			sum_shared += sigma2(gamma-1) * (this_shared / (2 * alpha(gamma-1)));
			
		// independent time after mrca
			sum_indep += alpha(gamma-1) * -(old_dates.each_row() - young_dates.t());			// matrix of independent time

		// vcv diagonal
			this_shared = exp(-2 * alpha(gamma-1) * shiftDates(gamma-1)) - exp(-2 * alpha(gamma-1) * young_dates);		// time is the base of the segment to the mrca node to the tip or top of segment interval
			VV_diag += sigma2(gamma-1) * (this_shared / (2 * alpha(gamma-1)));
		}		
		
		VV = exp(-(sum_indep+sum_indep.t())) % sum_shared;
		VV.diag() = VV_diag;
		return VV;
	}

// [[Rcpp::export]]
	arma::mat getW_mat(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec alpha) {	
	// *** assume shiftDates and alphas are sorted oldest to youngest
	// *** assume shiftDates and dates (node dates) are rescaled to time since root
		
		const unsigned int ntaxa = edges.col(0).min() - 1;
		
		arma::vec young_dates(ntaxa);
		arma::mat dur(shiftDates.n_elem-1, ntaxa);		// matrix of durations within each selective regime **** ONLY WORKS FOR TEMPORAL - NOT CLADEWISE *****
		arma::mat W(shiftDates.n_elem-1, ntaxa);		// matrix of weights
	
		for (uword gamma=1; gamma < shiftDates.n_elem; gamma++) {
			young_dates = dates.rows(0, ntaxa-1);
			young_dates.elem(find(young_dates > shiftDates(gamma))).fill(shiftDates(gamma));		// tips are younger than the segment interval, young dates are the top of the segment interval
			young_dates.elem(find(young_dates < shiftDates(gamma-1))).fill(shiftDates(gamma-1));	// tips are older than this segment interval,  young dates are the base of the segment interval
			dur.row(gamma-1) = young_dates.t() - shiftDates(gamma-1);
			W.row(gamma-1) = exp(alpha(gamma-1) * young_dates.t()) - exp(alpha(gamma-1) * shiftDates(gamma-1));
		}

		arma::vec alpha_T = exp(-(alpha.t() * dur).t());		// product of alpha vector at T, summed over all segments (alphas)

		W=cumsum(W,0) % (dur/dur);			// cumulative sum over regimes from root, also eliminates cells for which the branch was not in segment
		W=W.replace(arma::datum::nan, 0);	// replaces NaN (cells for which the branch was not in the segment) with 0
		
		W.each_row() %= alpha_T.t();
		W.insert_rows(0,alpha_T.t());
		W.each_row() /= sum(W, 0);			// each row entry in W is divided by the sum of its row to ensure that the weights for each species sum to 1

		return W.t();
	}
	
// [[Rcpp::export]]
	arma::mat getThetaHat(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2, arma::vec alpha, arma::mat X) {
		arma::mat WW = getW_mat(edges, dates, shiftDates, alpha);
		arma::mat VV = getVCV_OU(edges, dates, shiftDates, sigma2, alpha);
//		return VV.i();
//		return (WW.t() * VV.i() * WW);
		return (WW.t() * VV.i() * WW).i() * WW.t() * VV.i() * X;
	}
	
// [[Rcpp::export]]
	double getPhylogeneticMeans(arma::mat dat, arma::mat V) {
		arma::vec one_vec(dat.n_elem, fill::ones);
		arma::mat x = (one_vec.t() * V.i() * one_vec).i() * (one_vec.t() * V.i() * dat);
		return x(0,0);
	}
	
// [[Rcpp::export]]
	double getLikelihood (arma::mat E, arma::mat V, arma::mat X) {
		const int N = X.n_elem;
		const double pi = std::atan(1.0)*4;
		arma::mat a = X-E;
		arma::mat b = a.t() * V.i() * a;
		double b_ = b(0,0);

//		double c = exp(-0.5*b_);
//		double d = sqrt(pow((2*pi),N) * arma::det(V));
//		return log(c/d);

		return (b_+ ((N * log(2 * pi)) + log(arma::det(V))))/-2;
	}
	
// [[Rcpp::export]]
	double hotLikelihood_BM (arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2, arma::mat data) {
		arma::mat VV = getVCV_BM(edges, dates, shiftDates, sigma2);
		if (arma::det(VV)==0) { return std::numeric_limits<double>::infinity(); 
		} else {
			arma::mat EE(data.n_elem, 1);
			double pmean = getPhylogeneticMeans(data, VV);
			EE.fill(pmean);
			return -getLikelihood(EE, VV, data);
		}
	}
	
// [[Rcpp::export]]
	double hotLikelihood_OU (arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2, arma::vec alpha, arma::vec theta, arma::mat data) {
		arma::mat WW = getW_mat(edges, dates, shiftDates, alpha);
		arma::mat VV = getVCV_OU(edges, dates, shiftDates, sigma2, alpha);
		if (arma::det(VV)==0) { return std::numeric_limits<double>::infinity(); 
		} else return -getLikelihood(WW * theta, VV, data);
	}
	
