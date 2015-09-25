#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix findEdges(IntegerMatrix simps, int d, int n) {

	int nrows = d*(d-1)/2;
	IntegerMatrix out(nrows*n,2);
	
	int rowind = 0;
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < d; i++) {
			for (int j = i+1; j < d; j++) {
				out(rowind,0) = simps(k,i);
				out(rowind,1) = simps(k,j);
				rowind++;
			}
		}
	}
	
	return(out);

}