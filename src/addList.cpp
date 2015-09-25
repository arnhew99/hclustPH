#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List addList(List simplist_reps, List simplist_adds) {

	int nadds = simplist_reps.size();
	
	//List outlist = clone(simplist_reps);
	
	
	for (int i2 = 0; i2 < nadds; i2++) {
	
		IntegerVector x = simplist_reps[i2];
		IntegerVector y = simplist_adds[i2];
	
		int n = x.size();
		int m = y.size();
		
		IntegerVector tmp(n+m, 1);
		
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
		
		//outlist[i] = out;
		simplist_reps[i2] = out;
	
	}
	
	return(simplist_reps);
}