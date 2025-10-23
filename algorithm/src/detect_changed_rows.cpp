#include <Rcpp.h>
#include <thread>
#include <vector>
#include <algorithm>

//2025-09-08
//a file for fast detection of identical rows in a data frame of numeric vectors
//multithreaded on groups, by a grouping vector

// Helper to test equality including NA
inline bool check_equality_double_with_na(const double a, const double b) noexcept
{
    //common case first - if both are possible to compare, just do it
    if (a == b ) return true;
    
    //check for NA in each
    const bool a_na = Rcpp::traits::is_na<REALSXP>(a); 
    const bool b_na = Rcpp::traits::is_na<REALSXP>(b);
    return (a_na && b_na);
}

// Compare row index with row index-1 across all numeric columns
inline bool row_identical_to_previous(std::size_t index,
                                      const std::vector<const double *> &cols) noexcept
{
    // check equality per column and return false if any differ
    for (const auto &c : cols){
        if (!check_equality_double_with_na(c[index], c[index - 1u]))
            return false;
    }

    // default - return true
    return true;
}

//column sweep to mark changes
inline void mark_changes_in_group(std::size_t begin, std::size_t end,
                               const std::vector<const double*>& col_ptrs,
                               std::vector<uint8_t>& out) noexcept
{
    //mark the first element in the group as not identical
    if (end <= begin + 1) {
        if (begin < end) out[begin] = 0;
        return;
    }

    //mark the first element in the group as not identical
    out[begin] = 0;   
    
    // rows that could still be false
    std::size_t remaining = (end - begin - 1); 
    //sweep columns
    for (const auto &c : col_ptrs){
        //break the loop if nothing remains
        if (!remaining) break;

        //walk rows, only touch columns that haven't been marked as changed yet
        for (std::size_t j = begin+1u; j < end; ++j) {
            //already marked as non-identical by another column - 
            if (out[j] == 0) continue;

            if (!check_equality_double_with_na(c[j], c[j - 1u])) {
                out[j] = 0;
                // whole group resolved
                if (--remaining == 0) return; 
                
            }
        }
    }
}

// [[Rcpp::export]]
//a fast c++ function to highlight identical rows...
Rcpp::LogicalVector mark_identical_rows(Rcpp::NumericVector groups,
                                          Rcpp::DataFrame df,
                                          int nThreads = 0)
{
    //data size parameters
    const std::size_t n_obs = static_cast<std::size_t>(df.nrows());
    const std::size_t n_group_obs = static_cast<std::size_t>(groups.size());
    const std::size_t n_cols = static_cast<std::size_t>(df.size());

    //assert equal size of df and groups
    if (n_obs == 0 || n_group_obs == 0 || n_group_obs != n_obs){
        Rcpp::stop("'groups' and 'df' must be the same length and > 0");
    }

    //assert number of columns > 0
    if (n_cols == 0){
        Rcpp::stop("`df` must have at least one column");
    }

    //only one row - nothing to compare - return a false vector
    if (n_obs < 2) {
        // nothing to compare; all FALSE
        return Rcpp::LogicalVector(n_obs, false);
    }

    // Coerce all columns to NumericVector once; store raw pointers for thread-safe reads
    std::vector<Rcpp::NumericVector> num_cols;
    num_cols.reserve(n_cols);
    for (std::size_t c = 0; c < n_cols; ++c){
        if (TYPEOF(df[c]) != REALSXP) {
            Rcpp::stop("All columns in `df` must be numeric (REALSXP).");
        }
        Rcpp::NumericVector col = Rcpp::as<Rcpp::NumericVector>(df[c]);
        if (static_cast<std::size_t>(col.size()) != n_obs)
        {
            Rcpp::stop("All columns must have the same length as `df`");
        }
        num_cols.push_back(col);
    }

    //store pointers to columns for thread-safe reads
    std::vector<const double *> col_ptrs;
    col_ptrs.reserve(n_cols);
    for (auto &col : num_cols)
        col_ptrs.push_back(&(col[0]));

    //Group vector pointer
    const double *gptr = &(groups[0]);

    //Precompute contiguous group ranges [begin, end)
    std::vector<std::pair<std::size_t, std::size_t>> ranges;
    ranges.reserve(1024);
    std::size_t start = 0;
    for (std::size_t i = 1; i < n_obs; ++i)
    {
        if (!(check_equality_double_with_na(gptr[i], gptr[start])))
        {
            ranges.emplace_back(start, i);
            start = i;
        }
    }
    //final group emplaced
    ranges.emplace_back(start, n_obs);

    //Determine thread count
    std::size_t user_threads = (nThreads > 0) ? static_cast<std::size_t>(nThreads) : 0;
    std::size_t hw = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    std::size_t nthreads = (user_threads == 0) ? hw : std::min(user_threads, hw);

    //Sanity check to avoid too many threads if nobs is tiny
    nthreads = std::min(nthreads, std::max<std::size_t>(1, n_obs / 1000 + 1));
    //sanity check, maximum of one thread per group
    nthreads = std::min<std::size_t>(nthreads, ranges.size());

    //work array - will be converted to Rcpp::LogicalVector later on
    std::vector<uint8_t> out(n_obs, 1);

    //single-threaded path
    if (nthreads <= 1 || ranges.size() == 1)
    {
        // Single-thread path
        for (const auto &rg : ranges)
        {
            mark_changes_in_group(rg.first, rg.second, col_ptrs, out);


            //const std::size_t begin = rg.first;
            //const std::size_t end = rg.second;
            //if (begin < end)
            //    out[begin] = 0; // first in group
            //for (std::size_t j = begin + 1; j < end; ++j)
            //{
            //    out[j] = row_identical_to_previous(j, col_ptrs) ? 1 : 0;
            //}
        }
    }
    else
    {
        //split the ranges into buckets for the worker threads
        std::vector<std::thread> threads;
        threads.reserve(nthreads);
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> buckets(nthreads);
        //make sure the tasks are about equal in size
        const std::size_t bucket_size_target = n_obs / nthreads;
        std::size_t current_bucket = 0u;
        std::size_t current_bucket_size = 0u;

        //split the work among the threads
        for (std::size_t i = 0u; i < ranges.size(); ++i){
            buckets[current_bucket].push_back(ranges[i]);
            current_bucket_size += ranges[i].second - ranges[i].first;
            //if the current bucket contains enough rows - shift to a new one
            if ((current_bucket_size >= bucket_size_target) && (current_bucket < (nthreads -1))){
                current_bucket += 1u;
                current_bucket_size = 0u;
            }
        }

        //generate the threads
        for (std::size_t t = 0; t < nthreads; ++t)
        {
            //thread logic is identical to the single thread path above
            threads.emplace_back([&, t](){
                const auto& my_ranges = buckets[t];
                for (const auto& rg : my_ranges) {
                    mark_changes_in_group(rg.first, rg.second, col_ptrs, out);
                    //const std::size_t begin = rg.first;
                    //const std::size_t end = rg.second;
                    //if (begin < end)
                    //    out[begin] = 0; // first in group
                    //for (std::size_t j = begin + 1; j < end; ++j){
                    //    out[j] = row_identical_to_previous(j, col_ptrs) ? 1 : 0;
                    //}
                }
            }); 
        }

        for (auto &t : threads)
            t.join();
    }

    // Copy back to LogicalVector
    Rcpp::LogicalVector res(n_obs);
    for (std::size_t i = 0; i < n_obs; ++i)
        res[i] = out[i];

    return res;
}