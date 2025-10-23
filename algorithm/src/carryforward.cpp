#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::plugins("cpp14")]]

//A way to keep track of last active rate of multiple infusion pumps
NumericVector cpp_roll_sum_of_last(NumericVector X,
                   IntegerVector Times,
                   IntegerVector IDs,
                   int roll = 900){
  int n = X.size();
  int m;
  bool inserted;
  NumericVector out = NumericVector(n);
  //initiate a set of vectors to store the current last known rates
  std::vector<int> id;
  std::vector<int> time;
  std::vector<double> value;
  for (int i = 0;i<n;i++){
    m = id.size();
    inserted = false;
    //If there are no current values just insert
    if (m == 0){
      id.push_back(IDs[i]);
      time.push_back(Times[i]);
      value.push_back(X[i]);
      inserted = true;
    } else {
      //Loop through current elements
      for (int j = 0;j<m;j++){
        
        if (id[j] == IDs[i]){
          //we found an old record with the same ID - update that location
          id[j] = IDs[i];
          time[j] = Times[i];
          value[j] = X[i];
          inserted = true;
        } else if (time[j]+roll<Times[i]){
          //if an element has timed out, remove it from all vectors
          id.erase(id.begin()+j);
          time.erase(time.begin()+j);
          value.erase(value.begin()+j);
          //un-increment loop index
          j--;
          //reduce number of saved elements
          m--;
        }
      }
      //if it is still not inserted, just insert
      if (inserted == false) {
        id.push_back(IDs[i]);
        time.push_back(Times[i]);
        value.push_back(X[i]);
        inserted = true;
      }
    }
    out[i] = std::accumulate(value.begin(), value.end(), 0.0);
  }
  return out;
}

//This function calculates, at each timestamp, the total sum of "Active" rates, i.e. they are not older than "roll" seconds.
//[[Rcpp::export]]
NumericVector cpp_sum_of_last_rates(NumericVector Rate, NumericVector Time, IntegerVector ID, int roll){
  int n = ID.size();
  int m;
  bool inserted;
  
  //These vectors store the last known rates, the pump ID for each rate, and the timestamp of the last rate
  std::vector<double> rates(0);
  std::vector<int> ids(0);
  std::vector<double> times(0);
    
  NumericVector out = NumericVector(n);
  //The loop through all time points
  for (int i = 0;i<n;i++){
    //See how many old rates we have stored
    m = rates.size();
    //Initialize this current total rate at this time point
    inserted = false;
    out[i] = 0;
    //Loop through all old rates
    for (int j = 0;j<m;j++){
      //If we find the same ID in the old rates, replace that entry
      if (ids[j] == ID[i]){
        ids[j] = ID[i];
        rates[j] = Rate[i];
        times[j] = Time[i];
        //add this to the current sum
        out[i] = out[i]+Rate[i];
        inserted = true;
      } else
        //if an entry has timed out, delete it
        if (times[j]+roll<Time[i]){
          ids.erase(ids.begin()+j);
          rates.erase(rates.begin()+j);
          times.erase(times.begin()+j);
          //decrement the number of elements
          m--;
          //decrement loop index of inner loop
          j--;
        } else {
          //if the value is not too old, add in to the sum
          out[i] = out[i]+rates[j];
        }
    }
    //If the ID is new, add it
    if (!inserted){
      ids.push_back(ID[i]);
      rates.push_back(Rate[i]);
      times.push_back(Time[i]);
      out[i] = out[i] + Rate[i];
    }
  }    
  return out;
}

//A way of converting dose rates to NE equivalents according to lambden
double lambden_dose(double x, int id){
  //drugs that are not 226, 227 or 95 have a NE equivalent of 0
  double out = 0.0;
  if (id == 226 || id == 227){
    out = x;
  } else if (id == 95 || id == 1000511) {
    out = x*2.5;
  }
  return out;
}

//Similar to the function above, but calculating Norepi equivalents
//according to lambden
//Note that the variable is called 'rate' and 'rates' but it should actually be
//doses in mcg/kg/h, or units/h to be correct - this needs to be precalculated.
//[[Rcpp::export]]
NumericVector cpp_lambden_sum(NumericVector Rate, NumericVector Time, IntegerVector ID, int roll){
  int n = ID.size();
  int m;
  bool inserted;
  


  //These vectors store the last known rates, the pump ID for each rate, and the timestamp of the last rate
  std::vector<double> rates(0);
  std::vector<int> ids(0);
  std::vector<double> times(0);
  
  NumericVector out = NumericVector(n);
  //The loop through all rates
  for (int i = 0;i<n;i++){
    //See how many old rates we have stored
    m = rates.size();
    //Initialize this current total rate.
    inserted = false;
    out[i] = 0;
    //Loop through all old rates
    for (int j = 0;j<m;j++){
      //If we find the same ID in the old rates, replace that entry
      if (ids[j] == ID[i]){
        ids[j] = ID[i];
        //only calculate this once, for efficiency
        rates[j] = lambden_dose(Rate[i], ID[i]);
        times[j] = Time[i];
        //add this to the current sum
        out[i] = out[i]+rates[j];
        inserted = true;
      } else
        //if an entry has timed out, delete it
        if (times[j]+roll<Time[i]){
          ids.erase(ids.begin()+j);
          rates.erase(rates.begin()+j);
          times.erase(times.begin()+j);
          //decrement the number of elements
          m--;
          //decrement loop index of inner loop
          j--;
        } else {
          //if the value is not too old, add in to the sum
          out[i] = out[i]+rates[j];
        }
    }
    //If the ID is new, add it
    if (!inserted){
      ids.push_back(ID[i]);
      rates.push_back(lambden_dose(Rate[i], ID[i]));
      times.push_back(Time[i]);
      //the newly added, calculated dose...
      out[i] = out[i] + rates.back();
    }
  }    
  return out;
}




//Faster version that loops through x and p
//to keep track of last of x for each p
double cpp_sum_of_last(NumericVector x,
                        NumericVector p){
  NumericVector u = unique(p);
  int m = u.size();
  NumericVector last = NumericVector(m);
  int n = x.size();
  for (int i = 0;i<n;i++){
    for (int j = 0;j<m;j++){
      if (p[i] == u[j]){
        last[j] = x[i];
        break;
      }
    }
  }
  double res = sum(last);
  return res;
}
//Faster version that loops through x and p
//to keep track of last of x for each p
double cpp_lambden_vd_sum_of_last(NumericVector x,
                       NumericVector p){
  NumericVector u = unique(p);
  int m = u.size();
  NumericVector last = NumericVector(m);
  int n = x.size();
  for (int i = 0;i<n;i++){
    for (int j = 0;j<m;j++){
      if (p[i] == u[j]){
        last[j] = x[i];
        break;
      }
    }
  }
  double res = 0;
  for (int i = 0;i<m;i++){
    if (u[i] == 226){
      res = res+last[i];
    }
    if (u[i] == 227){
      res = res+last[i];
    }
    if (u[i] == 95){
      res = res+last[i]*2.5;
    }
  }
  return res;
}


//[[Rcpp::export]]
Rcpp::NumericVector cpp_cf(NumericVector x, NumericVector t, double roll = 0.0) {
  auto n = x.size();
  double last = NA_REAL;
  double last_t = 0.0;
  NumericVector out(n);
  for (auto i = 0u;i<n;++i){
    if (!NumericVector::is_na(x[i])){
      last = x[i];
      last_t = t[i];
      out[i] = x[i];
    } else {
      if (t[i]<=last_t+roll){
        out[i] = last;
      } else {
        out[i] = NA_REAL;
      }
    }
  }
  return out;
}

