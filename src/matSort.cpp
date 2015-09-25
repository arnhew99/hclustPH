#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix matSort(IntegerMatrix X, int n, int d) {

	IntegerMatrix Y(n,d);

	for (int i = 0; i < n; i++) {
	
		IntegerVector tmp(d);
		for (int j = 0; j < d; j++) {
			tmp[j] = X(i,j);
		}
		IntegerVector tmp2 = clone(tmp);
		std::sort(tmp2.begin(), tmp2.end());
		for (int j = 0; j < d; j++) {
			Y(i,j) = tmp2[j];
		}
	}

	return(Y);

} 