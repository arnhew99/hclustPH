#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector topsimps(List simplist) {

	int n = simplist.size();
	
	IntegerVector out(n,-1);
	
	for (int i = 0; i < n; i++) {
		IntegerVector tmp = simplist[i];
		int m = tmp.size();
		if (m > 0) {
			out[i] = max(tmp); 
		}
	}
	
	return(out);

}