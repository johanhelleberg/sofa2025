//2020-11-15
//This file contains fast functions for median, median absolute deviation, and
//a hampel function that returns a vector where element 0 = median, element 1 = mad.

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double cpp_med(Rcpp::NumericVector xx) {
  Rcpp::NumericVector x = Rcpp::clone(xx);
  std::size_t n = x.size() / 2;
  std::nth_element(x.begin(), x.begin() + n, x.end());
  
  if (x.size() % 2) return x[n]; 
  return (x[n] + *std::max_element(x.begin(), x.begin() + n)) / 2.;
}

// [[Rcpp::export]]
NumericVector cpp_quantiles(NumericVector xx, 
                            NumericVector probs = NumericVector::create(0.25,0.75)){
  NumericVector x = clone(xx);
  int n = x.size();
  int pn = probs.size();
  std::size_t maxsort = ceil(1+(n-1)*probs[which_max(probs)]);
  std::size_t minsort = floor(0+(n-1)*probs[which_min(probs)]);
  //return minsort;
  
  std::nth_element(x.begin(), x.begin()+minsort, x.end());
  std::partial_sort(x.begin()+minsort, x.begin()+maxsort, x.end());
  
  NumericVector qs = NumericVector(pn);
  double index;
  int low;
  int high;
  double interp;
  
  for (int i = 0;i<pn;i++){
    index = (n-1)*probs[i];
    low = std::floor(index);
    high = std::ceil(index);
    
    interp = index - low;
    
    qs[i] = x[low]*(1-interp)+x[high]*interp;
  }
  return qs;
}

// [[Rcpp::export]]
NumericVector cpp_tukey(NumericVector x, double k = 1.5){
  NumericVector qs = cpp_quantiles(x);
  NumericVector out = NumericVector(2);
  double IQR = qs[1]-qs[0];
  out[0] = qs[0]-k*IQR;
  out[1] = qs[1]+k*IQR;
  return out;
}

// [[Rcpp::export]]
double cpp_mean(Rcpp::NumericVector x){
  double cumsum = 0;
  int n = x.size();
  for (int i = 0;i<n;i++){
    cumsum+=x[i];
  }
  return cumsum/n;
}

// [[Rcpp::export]]
double cpp_lambden_vd_sum(Rcpp::NumericVector x, 
                         Rcpp::NumericVector id){
  int n = x.size();
  
  if (n == 0){
    return 0;
  }
  
  //initialize the sums
  double na = 0;
  int na_n = 0;
  double ad = 0;
  int ad_n = 0;
  double vp = 0;
  int vp_n = 0;
  
  //sum in a loop
  for (int i = 0; i<n;i++){
    if (id[i] == 226){
      na = na + x[i]; 
      na_n++;
    }
    if (id[i] == 227){
      ad = ad + x[i]; 
      ad_n++;
    }
    if (id[i] == 95){
      vp = vp + x[i]; 
      vp_n++;
    }
    
  }
  if (ad_n == 0) ad_n = 1;
  if (na_n == 0) na_n = 1;
  if (vp_n == 0) vp_n = 1;
  
  return (na/na_n)+(ad/ad_n)+(vp/vp_n)*2.5;
}
// [[Rcpp::export]]
NumericVector cpp_doseAtT(NumericVector x, NumericVector time, double halflife){
  int n = x.size();
  NumericVector out = NumericVector(n);
  double log2 = log(2.0);
  out = x*exp((-1.0 * time / halflife) * log2);
  
  //for (int i = 0;i<n;i++){
  //  out[i] = x[i]*exp((-1.0 * time[i] / halflife) * log(2.0));
  //}
  return(out);
}


// [[Rcpp::export]]
double cpp_mad(Rcpp::NumericVector xx, double k= 1.4826){
  //Get the median and size of the vector
  Rcpp::NumericVector x = Rcpp::clone(xx);
  int n = x.size();
  double median = cpp_med(x);
  //Initialize the array of MADs
  Rcpp::NumericVector absdev(n);
  
  //Rcpp::NumericVector x2 = Rcpp::clone(xx);
  
  for(int i = 0; i < n; ++i) {
    absdev[i] = abs(x[i]-median);
  }
  
  double mad = cpp_med(absdev)*k;
  
  return mad;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_hampel(Rcpp::NumericVector xx, double k = 1.4826){
  Rcpp::NumericVector x = Rcpp::clone(xx);
  int n = x.size();
  double median = cpp_med(x);
  
  //Initialize the array of MADs
  Rcpp::NumericVector absdev(n);
  
  //Rcpp::NumericVector x2 = Rcpp::clone(xx);
  
  for(int i = 0; i < n; ++i) {
    absdev[i] = abs(x[i]-median);
  }
  
  double mad = cpp_med(absdev)*k;
  
  Rcpp::NumericVector ret(2);
  //ret["median"] = median;
  //ret["mad"] = mad;
  ret[0] = median;
  ret[1] = mad;
  return ret;
}


// [[Rcpp::export]]
NumericVector cpp_tukey_slide(NumericVector x, 
                               NumericVector t, 
                               int before = 900,
                               int after = 900,
                               double k = 1.5) {
  int n = x.size();
  NumericVector out = NumericVector(n);
  int start_pos = 0;
  int end_pos = -1;
  std::vector<double> values;
  
  double low, high, interp, index, iqr;
  int m;
  //initialize the vector of probabilities and the quantiles
  std::vector<double> probs = {0.25,0.75};
  std::vector<double> qs(2);
  
  //the main loop
  for (int i = 0;i<n;i++){
    //Section one, generating the sliding, sorted values-vector
    //Increment the end of the window
    while (end_pos<n-1){
      if (t[end_pos+1]<=t[i]+after){
        end_pos++;
        //Insert sorted into the values-vector
        
        values.insert(std::upper_bound(values.begin(),
                                       values.end(),
                                       x[end_pos]),
                                       x[end_pos]);
        
      } else break;
    }
    //increment the beginning of the window
    while(t[start_pos]<t[i]-before){
      //Delete one value equal to x[start_pos]
      values.erase(std::find(values.begin(), values.end(), x[start_pos]));
      start_pos++;
    }
    
    m = values.size();
    //Section two, the outlier algorithm
    //first, calculate the 25th and 75th percentile
    //(that's a painful lambda expression :D :D )
    std::transform(probs.begin(), 
                   probs.end(), 
                   qs.begin(), 
                   [&m, 
                    &index, 
                    &low, 
                    &high, 
                    &interp, 
                    &values](double &p){
                      index = (m-1)*p;
                      low = std::floor(index);
                      high = std::ceil(index);
                      interp = index - low;
                      return values[low]*(1-interp)+values[high]*interp;});
    //calculate the iqr*k
    iqr = (qs[1] - qs[0])*k;
    //Tukey's fences = replace outliers with Q25-IQR*k or Q75+IQR*k
    if (x[i] > qs[1]+iqr) out[i] = qs[1]+iqr; else
      if (x[i] < qs[0] - iqr) out[i] = qs[0]-iqr; else
        out[i] = x[i];
  }
  return out;
}
