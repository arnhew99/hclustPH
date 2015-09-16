#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector addSimplices(IntegerVector x, IntegerVector y) {

	int n = x.size();
	int m = y.size();
	
	IntegerVector tmp(n+m);
	
	for (int i = 0; i< n+m; ++i) {
		tmp[i] = 1;
	}
	
	int count = n+m;
	
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			if (x[i] == y[j]) {
				tmp[i] = 0;
				count--;
			} 		
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			if (y[i] == x[j]) {
				tmp[i+n] = 0;
				count--;
			} 		
		}
	}
	
	IntegerVector out(count);
	
	int curout = 0;
	for (int i = 0; i< n+m; ++i) {
	
		if (tmp[i] == 1) {
	
			if (i < n) {
				out[curout] = x[i];
			} else {
				out[curout] = y[i-n];
			}
			curout++;
		}
	}
	
	return(out);
}