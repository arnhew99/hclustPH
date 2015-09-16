#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector hclustIntermediate(IntegerVector l, IntegerVector r) {

	int n = l.size();
	
	int alter = 1;
	while (alter == 1) {
		alter = 0;
		for (int i = n-1; i >= 0; i--) {
		
			int cur = l[i];
			
			for (int j = i-1; j >= 0; j--) {
				
				if (cur == r[j]) {
					cur = l[j];
					alter = 1;
				}
				
			}
			l[i] = cur;
		
		}
	}
	return(l);
}