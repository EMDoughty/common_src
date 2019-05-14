#ifndef __CHR_CPP__
#define __CHR_CPP__

#include <RcppArmadillo.h> 

	arma::uword getAnc(arma::uword desc, arma::mat edges);
	arma::uword getMRCA(arma::uword tip1, arma::uword tip2, arma::mat edges);
	arma::umat mrca_cpp(arma::mat edges);
	arma::mat getMRCADates(arma::mat edges, arma::vec dates);	
	
	arma::mat getVCV_BM(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2);
	arma::mat getVCV_OU(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2, arma::vec alpha);
	
	arma::mat getW_mat(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec alpha);
	arma::mat getThetaHat(arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2, arma::vec alpha, arma::mat X);
	double getPhylogeneticMeans(arma::mat dat, arma::mat V);
	
	double getLikelihood (arma::mat E, arma::mat V, arma::mat X);
	double hotLikelihood_BM (arma::mat edges, arma::vec dates, arma::vec shiftDates, arma::vec sigma2, arma::mat data);
	double hotLikelihood_OU (arma::mat edges, arma::vec dates, arma::vec model, arma::vec shiftDates, arma:: vec sigma2, arma::vec alpha, arma::vec theta, arma::mat data);

#endif	// __CHR_CPP__
