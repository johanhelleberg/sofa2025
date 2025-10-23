#include <Rcpp.h>
//[[Rcpp::plugins("cpp14")]]


//[[Rcpp::export]]
double cpp_min_na(Rcpp::NumericVector x){
    bool found = false;
    auto lowest = NA_REAL;
    for (auto &v : x){
        if(!Rcpp::NumericVector::is_na(v)){
            if (!found){
                found = true;
                lowest = v;
            } else if (v < lowest){
                lowest = v;
            }
        }
    }
    return lowest;
}

//[[Rcpp::export]]
Rcpp::NumericVector remove_neg_dups(Rcpp::NumericVector x){
	Rcpp::NumericVector out;
	auto n = x.size();
	//keep track of removed indices
	std::vector<uint8_t> removed(n,0u);
	
	for (auto i = 0u;i<n;++i){
		for (auto j = i;j<n;++j){
			//remove a pair if found
			if (i != j && ((x[j] == x[i] * -1.0) && (removed[j] == 0u))){
				removed[i] = 1u;
				removed[j] = 1u;
				break;
			}
		}
		//if we're not removed, add us to the result
		if (!removed[i]){
			out.push_back(x[i]);
		}
		
	}
	
    return out;
}

//[[Rcpp::export]]
Rcpp::IntegerVector mark_neg_dups(Rcpp::NumericVector x){
	auto n = x.size();
	Rcpp::IntegerVector out(n,0u);
	
	for (auto i = 0u;i<n;++i){
		for (auto j = i;j<n;++j){
			//remove a pair if found
			if (i != j && ((x[j] == x[i] * -1.0) && (out[j] == 0u))){
				out[i] = 1u;
				out[j] = 1u;
				break;
			}
		}
		
	}
	
    return out;
}