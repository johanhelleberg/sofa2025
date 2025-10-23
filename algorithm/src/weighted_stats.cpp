#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <string>
#include <Rcpp.h>
//[[Rcpp::plugins("cpp14")]]


// a simple templated struct to keep a value and weight together
struct VW
{
    double value;
    double weight;
};

//a fast weighted median i
template <typename V>
double weighted_median(const V &v, const V &w, bool midpoint_interpolate = false)
{
    const std::size_t n = v.size();

    // error check the data
    if (n != static_cast<size_t>(w.size())){
        throw std::invalid_argument("Values and Weights must have the same size");
    }

    // return NA if size 0
    if (n == 0){
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<VW> vw;
    vw.reserve(n);

    // sum the weights at this point to avoid another loop for this
    double w_sum = 0;

    // for branchless exception handling
    std::size_t nan_index1 = 0;
    std::size_t neg_index1 = 0;

    // build the structure of value/weight pairs
    for (std::size_t i = 0; i < n; ++i){
        const double vi = static_cast<double>(v[i]);
        const double wi = static_cast<double>(w[i]);
        const bool is_nan = std::isnan(vi) || std::isnan(wi);
        const bool is_neg = wi < 0.0;

        nan_index1 += std::size_t(is_nan & (nan_index1 == 0)) * (i + 1);
        neg_index1 += std::size_t(is_neg & (neg_index1 == 0)) * (i + 1);

        if (!is_nan && !is_neg && wi > 0.0)
        {
            vw.push_back({vi, wi});
            w_sum += wi;
        }
    }

    // deferred exception handling - only the first error is reported
    if (neg_index1){
        const std::size_t i = neg_index1 - 1;
        throw std::invalid_argument("Weights must be non-negative (first at index " + std::to_string(i) + ")");
    }
    if (nan_index1){
        const std::size_t i = nan_index1 - 1;
        throw std::invalid_argument("NaN found in values or weights (first at index " + std::to_string(i) + ")");
    }
    if (vw.empty() || w_sum == 0.0){
        throw std::invalid_argument("At least one weight must be > 0");
    }

    // sort the vector
    std::sort(vw.begin(), vw.end(), [](const VW &a, const VW &b)
              { return a.value < b.value; });

    // set up target for aggregation
    const double w_target = w_sum * 0.5;
    double cum_w_sum = 0;

    // split the logic here depending on path
    if (!midpoint_interpolate){

        // aggregate up to midpoint
        for (const auto &p : vw)
        {
            cum_w_sum += p.weight;
            if (cum_w_sum >= w_target)
            {
                return p.value;
            }
        }

        // fallback - shoudln't happen
        return vw.back().value;
    }

    // midpoint interpolation - slightly more complex as it needs to get unique values in vw
    // setup the unique value vector
    std::vector<VW> u;
    u.reserve(vw.size());
    for (const auto &p : vw)
    {
        if (u.empty() || p.value != u.back().value)
        {
            u.push_back(p);
        }
        else
        {
            u.back().weight += p.weight;
        }
    }
    // numeric tolerance to find 'midpoint'
    const double eps = std::max(1.0, w_sum) * 1e-12;

    // the median loop
    for (std::size_t i = 0u; i < u.size(); ++i)
    {
        cum_w_sum += u[i].weight;
        // strictly crosset the target
        if (cum_w_sum > w_target + eps)
        {
            return u[i].value;
        }

        // exacly at 50%, interpolate between this and the next distinct value (if there is any)
        if (std::abs(cum_w_sum - w_target) <= eps)
        {
            if (i + 1 < u.size())
            {
                return 0.5 * (u[i].value + u[i + 1].value);
            }
            else
            {
                return u[i].value;
            }
        }
    }
    // fallback - shouldn't happen
    return u.back().value;
}


//a fast weighted quantiles
template <typename V>
std::vector<double> weighted_quantiles(const V &v, const V &w, const V &qs, bool midpoint_interpolate = false)
{
    const std::size_t n = v.size();
    const std::size_t nq = qs.size();
    
    // error check the data
    if (n != static_cast<size_t>(w.size())){
        throw std::invalid_argument("Values and Weights must have the same size");
    }

    // return NA if size 0
    if (n == 0 || nq == 0){
        throw std::invalid_argument("Values and Quantiles must contain > 0 elements");
    }

    //unique quantiles
    std::vector<double> uq;
    uq.resize(nq);

    //copy the data
    std::copy(qs.begin(), qs.end(), uq.begin());
    std::sort(uq.begin(), uq.end());
    double last_q = -std::numeric_limits<double>::infinity();
    for (auto q : uq){
        if (q <=0 || q >= 1.0){
            throw std::invalid_argument("Quantiles must fall in range (0-1)");
        }
        if (std::fabs(q-last_q) <= 1e-12){
            throw std::invalid_argument("Quantile values must be unique");
        }
        if (std::isnan(q)){
            throw std::invalid_argument("Quantiles must not include NA/NaN values");
        }
        last_q = q;
    }

    //from now on uq is a sorted vector of unique non-nan values
    std::vector<VW> vw;
    vw.reserve(n);

    // sum the weights at this point to avoid another loop for this
    double w_sum = 0;

    // for branchless exception handling
    std::size_t nan_index1 = 0;
    std::size_t neg_index1 = 0;

    // build the structure of value/weight pairs
    for (std::size_t i = 0; i < n; ++i){
        const double vi = static_cast<double>(v[i]);
        const double wi = static_cast<double>(w[i]);
        const bool is_nan = std::isnan(vi) || std::isnan(wi);
        const bool is_neg = wi < 0.0;

        nan_index1 += std::size_t(is_nan & (nan_index1 == 0)) * (i + 1);
        neg_index1 += std::size_t(is_neg & (neg_index1 == 0)) * (i + 1);

        if (!is_nan && !is_neg && wi > 0.0)
        {
            vw.push_back({vi, wi});
            w_sum += wi;
        }
    }

    // deferred exception handling - only the first error is reported
    if (neg_index1){
        const std::size_t i = neg_index1 - 1;
        throw std::invalid_argument("Weights must be non-negative (first at index " + std::to_string(i) + ")");
    }
    if (nan_index1){
        const std::size_t i = nan_index1 - 1;
        throw std::invalid_argument("NaN found in values or weights (first at index " + std::to_string(i) + ")");
    }
    if (vw.empty() || w_sum == 0.0){
        throw std::invalid_argument("At least one weight must be > 0");
    }

    // sort the vector
    std::sort(vw.begin(), vw.end(), [](const VW &a, const VW &b)
              { return a.value < b.value; });

    // set up target for aggregation
    std::size_t q_idx = 0;
    double w_target = w_sum * uq[q_idx];

    double cum_w_sum = 0;

    //out vector, reserve to number of quantiles
    std::vector<double> out;
    out.reserve(nq);


    // split the logic here depending on path
    if (!midpoint_interpolate){

        // aggregate up to midpoint
        for (const auto &p : vw)
        {
            cum_w_sum += p.weight;
            while (cum_w_sum >= w_target)
            {
                out.push_back(p.value);
                //update target to next weight
                q_idx+=1;
                if (q_idx < nq){
                    w_target = w_sum * uq[q_idx];
                } else break;
                
            }
        }
        
        if (out.size() == nq) {
            return out;
        } else {
            throw std::invalid_argument("This set of data doesn't generate the requested amount of quantiles");
        }
    }

    // midpoint interpolation - slightly more complex as it needs to get unique values in vw
    // setup the unique value vector
    std::vector<VW> u;
    u.reserve(vw.size());
    for (const auto &p : vw)
    {
        if (u.empty() || p.value != u.back().value)
        {
            u.push_back(p);
        }
        else
        {
            u.back().weight += p.weight;
        }
    }
    // numeric tolerance to find 'midpoint'
    const double eps = std::max(1.0, w_sum) * 1e-12;

    // the quantile loop
    for (std::size_t i = 0u; i < u.size(); ++i)
    {
        cum_w_sum += u[i].weight;
        // strictly crossed the target
        while (cum_w_sum > w_target + eps)
        {
            //add the current value
            out.push_back(u[i].value);
            //update target to next weight
            q_idx+=1;
            if (q_idx < nq){
                    w_target = w_sum * uq[q_idx];
            } else break;
        }

        // exactly at 50%, interpolate between this and the next distinct value (if there is any)
        //note that while isn't needed as quantiles are guaranteed to be unique
        if (std::abs(cum_w_sum - w_target) <= eps)
        {
            if (i + 1 < u.size())
            {
                //add the current value
                out.push_back(0.5 * (u[i].value + u[i + 1].value));
                //update target to next weight
                q_idx+=1;
                if (q_idx < nq){
                    w_target = w_sum * uq[q_idx];
                } else break;
            }
            else
            {
                //add the current value
                out.push_back(u[i].value);
                //update target to next weight
                q_idx+=1;
                if (q_idx < nq){
                    w_target = w_sum * uq[q_idx];
                } else break;
            }
        }
    }
    //return the vector

    if (out.size() == nq) {
        return out;
    } else {
        throw std::invalid_argument("This set of data doesn't generate the requested amount of quantiles");
    }
}

template <typename V>
double weighted_mean(const V &v, const V &w){
    std::size_t n = v.size();

    // error check the data
    if (n != static_cast<size_t>(w.size())){
        throw std::invalid_argument("Values and Weights must have the same size");
    }

    // return NA if size 0
    if (n == 0){
        return std::numeric_limits<double>::quiet_NaN();
    }


    double w_sum = 0;
    double sum = 0;

    // for branchless exception handling
    std::size_t nan_index1 = 0;
    std::size_t neg_index1 = 0;

    //run the cumulative sum function
    for (std::size_t i = 0; i < n; ++i){
        const double vi = static_cast<double>(v[i]);
        const double wi = static_cast<double>(w[i]);
        const bool is_nan = std::isnan(vi) || std::isnan(wi);
        const bool is_neg = wi < 0.0;

        nan_index1 += std::size_t(is_nan & (nan_index1 == 0)) * (i + 1);
        neg_index1 += std::size_t(is_neg & (neg_index1 == 0)) * (i + 1);

        w_sum += wi;
        sum += vi*wi;


    }

    // deferred exception handling - only the first error is reported
    if (neg_index1){
        const std::size_t i = neg_index1 - 1;
        throw std::invalid_argument("Weights must be non-negative (first at index " + std::to_string(i) + ")");
    }
    if (nan_index1){
        const std::size_t i = nan_index1 - 1;
        throw std::invalid_argument("NaN found in values or weights (first at index " + std::to_string(i) + ")");
    }
    if (w_sum <= 1e-12){
        throw std::invalid_argument("Sum of weights must be greater than 0");
    }
    //return
    return sum / w_sum;
}



// [[Rcpp::export]]
double cpp_weighted_median(Rcpp::NumericVector values,
                           Rcpp::NumericVector weights,
                           bool interpolate_midpoint = false) {
  try {
    return weighted_median(values, weights, interpolate_midpoint);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_weighted_quantiles(Rcpp::NumericVector values,
    Rcpp::NumericVector weights,
    Rcpp::NumericVector quantiles,
    bool interpolate_midpoint = false) {
  try {
    return  Rcpp::wrap(weighted_quantiles(values, weights, quantiles, interpolate_midpoint));
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}

// [[Rcpp::export]]
double cpp_weighted_mean(Rcpp::NumericVector values,
                           Rcpp::NumericVector weights) {
  try {
    return weighted_mean(values, weights,);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}