#include <Rcpp.h>
//[[Rcpp::plugins("cpp11")]]
//[[Rcpp::depends(BH)]]
#include <boost/math/tools/minima.hpp>
//let's GOOO!!!!!!
#include <boost/circular_buffer.hpp>

#include <iostream>

#include <vector>
#include <array>
#include <deque>
#include <queue>
#include <tuple>

#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <future>

#include <string>
#include <chrono>


//2021-04-15
//A dynamically resizing version of the old pool
class threadpool_d
{
public:
	using Task = std::function<void()>;

	explicit threadpool_d(std::size_t max_threads, std::size_t new_thread_threshhold = 8u, std::size_t new_thread_multiplier = 2u, std::size_t min_threads = 1u)
	{
		minThreads = min_threads;
		maxThreads = max_threads;
		threadThreshhold = new_thread_threshhold;
		threadMultiplier = new_thread_multiplier;

		start();
	}
	~threadpool_d()
	{
		stop();
	}

	std::size_t get_size() {
		return threads.size();
	}
    //a debugging method for adding tasks that don't return anything...
	void add_task_lao(Task task)
	{
		{
			std::unique_lock<std::mutex> lock(event_mutex);
			tasks.emplace(std::move(task));
		}
		event_var.notify_one();
	}
    //returns an std::future of the return type of the functor passed onto task...
	template<class T>
	auto add_task(T task)->std::future<decltype(task())>
	{
		//a shared pointer to a packaged_task containing the ()-operator of the task-functor...
        auto wrapper = std::make_shared<std::packaged_task<decltype(task()) ()>>(std::move(task));
		std::size_t waiting_tasks;
		{
			std::unique_lock<std::mutex> lock(event_mutex);
            //get the amount of tasks currently on the pool, to see if we should start new threads
			waiting_tasks = tasks.size();
            //add the functor to the tasks-queue
			tasks.emplace([=] {
				(*wrapper)();
				});

        //This appears to have to be within the same block as the mutex to avoid race conditions that would sometimes
        //crash the pool when the tasks are very small and there are many cores
		    event_var.notify_one();
			
		}
        //the 'dynamic size' is very simplistic and basically starts threads until the size is maximum, 
        //if waiting_tasks > threadThreshhold + threadMultiplier * current amount of active threads
        //i.e. if each thread has more than threadMultiplier tasks yet to be completed
		if ((waiting_tasks > (threadThreshhold +  threadMultiplier * threads.size()) ) & (threads.size() < maxThreads)){
			add_thread();
		}
        //get the future from our shared pointer
		return wrapper->get_future();
	}

private:
	std::vector<std::thread> threads;
	std::condition_variable event_var;
	std::mutex event_mutex;
	bool stop_threads = false;

	std::size_t minThreads;
	std::size_t maxThreads;
	std::size_t threadThreshhold;
	std::size_t threadMultiplier;

	std::queue<Task> tasks;


	void start()
	{
		for (auto i = 0u; i < minThreads; ++i)
		{
			add_thread();
		}
	}

	void add_thread(){
        //this lambda contains the thread logic - 
		threads.emplace_back([=] {
				Task my_task;
				auto tasks_completed = 0u;
				while (true)
				{
					{
						std::unique_lock<std::mutex> lock(event_mutex);
						event_var.wait(lock, [=] {return stop_threads || !tasks.empty(); });

						if (stop_threads && tasks.empty()) break;

						my_task = std::move(tasks.front());
						tasks.pop();
					}
					tasks_completed++;
					//for diagnostic purposes - let's add this clock function
					//auto task_starttime = std::chrono::high_resolution_clock::now();
					my_task();
					//auto task_endtime = std::chrono::high_resolution_clock::now();
					//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(task_endtime - task_starttime);
					//std::cout << "Wo Shi Thread " << std::this_thread::get_id() << ". Wo completed " << tasks_completed << " Ge tasks. Runtime = " << duration.count() << std::endl;

				}
				});

	}

	void stop() noexcept
	{
		{
			std::unique_lock<std::mutex> lock(event_mutex);
			stop_threads = true;
		}
		event_var.notify_all();

		for (auto &t : threads) {
			t.join();
		}
	}
};

//Tukey's fences
class slide_tukey{
private:
	//the vector pointers
    Rcpp::NumericVector *values;
    Rcpp::NumericVector *times;
	
	//pointer to the start of our spans in the vectors
    Rcpp::NumericVector::iterator v_start;
    Rcpp::NumericVector::iterator t_start;
    //size of the spans in the vectors
	std::size_t size;

    //the return vector
    std::vector<double> res;
    //the start of our return vector
	std::vector<double>::iterator r_start;

	//sliding time range
	double before;
	double after;

	//the constant
	double k;


public:
	slide_tukey(Rcpp::NumericVector *v, Rcpp::NumericVector *t, int start, std::size_t group_size, double before_time, double after_time, double K) {
		//initialize our vector pointers
		values = v;
		times = t;
		res.resize(group_size);
		v_start = values->begin()+start;
		t_start = times->begin()+start;
		r_start = res.begin();

		size = group_size;
		before = before_time;
		after = after_time;
		k = K;
	}

	auto operator ()() -> std::vector<double> {
		//The iterators for the sliding time window
		auto window_start = t_start;
		auto window_end = t_start;

        //the end of the times vector
		auto t_end = t_start + size;

        //the current window size
		std::size_t window_size;



		//The vector of values within the time range
		std::vector<double> vals;
        //reserve some memory and push back the first element due to the way the loop is written below
        vals.reserve(64);
        vals.push_back(*v_start);
        
        //the value to be tested when windows are enlarged or shrunk
        double current_val;

		//the values for tukey's algorithm
		std::array<double, 2> probs = { 0.25, 0.75 };
		std::array<double, 2> qs = { 0.0, 0.0 };
		double low, high, interp, index, iqr, test_value;

        //the location in the results vector for the current loop index
		auto vector_index = std::distance(t_start, t_start);

		for (auto it = t_start; it != t_end; ++it) {
			//see if we're not out-of-bounds with the end of the vector
			while (window_end != t_end - 1) {
				//if the time at the end of the window is less than current time + after, add the corresponding value sorted into the values vector
				if (*(window_end + 1) <= *it + after) {
					++window_end;
                    current_val = *(v_start + std::distance(t_start, window_end));
                    vals.insert(std::upper_bound(vals.begin(),
						vals.end(),
						current_val),
						current_val);
					
				}
				else break;

			}

			//see if it's time to increment the beginning of the window
			while (*window_start < *it - before) {
                //if so, deleate one value equal to the one that just went out of the window
                vals.erase(std::find(vals.begin(), vals.end(), *(v_start + std::distance(t_start, window_start))));
				++window_start;
			}

			//calculate the size of the current window
			window_size = vals.size();

			//the value at the current time point
			test_value = *(v_start + std::distance(t_start, it));

			//calculate the quantiles - trivial since the vector is already sorted
			std::transform(probs.begin(),
				probs.end(),
				qs.begin(),
				[&window_size,
				&index,
				&low,
				&high,
				&interp,
				&vals](double& p) {
					index = (window_size - 1) * p;
					low = std::floor(index);
					high = std::ceil(index);
					interp = index - low;
					return vals[low] * (1 - interp) + vals[high] * interp; });


			//calculate the iqr*k
			iqr = (qs[1] - qs[0]) * k;

            //recalculate the vector index
			vector_index = std::distance(t_start, it);
            //Tukey's fences - replace values outside of {q25 - k*IQR, q75 + k*IQR}
			if (test_value > qs[1] + iqr) *(r_start + vector_index) = qs[1] + iqr; else
				if (test_value < qs[0] - iqr) *(r_start + vector_index) = qs[0] - iqr; else
					*(r_start + vector_index) = test_value;

		}
		return res;
	}
};

//this needs to be a global function since it seems difficult to capture a member function as lambda...
//var of a std::vector<double>
double var(std::vector<double>  & vec){
    double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);
    double m =  sum / vec.size();
    double accum = 0.0;
    std::for_each (std::begin(vec), std::end(vec), [&](const double d) {
        accum += (d - m) * (d - m);
    });

    double variance = accum / (vec.size()-1);
    return variance;
}

//Horn's algorithm
class slide_horn{
private:
	//the vector pointers
    Rcpp::NumericVector *values;
    Rcpp::NumericVector *times;
	
	//pointer to the start of our spans in the vectors
    Rcpp::NumericVector::iterator v_start;
    Rcpp::NumericVector::iterator t_start;
    //size of the spans in the vectors
	std::size_t size;

    //the return vector
    std::vector<double> res;
    //the start of our return vector
	std::vector<double>::iterator r_start;

	//sliding time range
	double before;
	double after;

    //the lambda ranges
    double min_lambda;
    double max_lambda;

	//the constant
	double k;



public:
	slide_horn(Rcpp::NumericVector *v, Rcpp::NumericVector *t, int start, std::size_t group_size, double before_time, double after_time, double K, double lambdamin = -2.0, double lambdamax = 2.0) {
		//initialize our vectors
		values = v;
		times = t;
		res.resize(group_size);
		v_start = values->begin()+start;
		t_start = times->begin()+start;
		r_start = res.begin();

        //initialize the constants
		size = group_size;
		before = before_time;
		after = after_time;
		k = K;
 
        min_lambda = lambdamin;
        max_lambda = lambdamax;
	}

	auto operator ()() -> std::vector<double> {
		//The iterators for the sliding time window
		auto window_start = t_start;
		auto window_end = t_start;

        //the end of the times vector
		auto t_end = t_start + size;

        //the current window size
		std::size_t window_size;

		//The vector of values within the time range
		std::vector<double> vals;
        //reserve some memory and push back the first element due to the way the loop is written below
        vals.reserve(64);
        vals.push_back(*v_start);
        
        //the value to be tested when windows are enlarged or shrunk
        double current_val;

		//the values for tukey's algorithm
		std::array<double, 2> probs = { 0.25, 0.75 };
		std::array<double, 2> qs = { 0.0, 0.0 };
		double low, high, interp, index, iqr, test_value;

        //values for the boxcox transformation
        double bestlambda, logsum, variance;
        //the vectors of boxcox-transformed values
        std::vector<double> v;
        std::vector<double> logv;
        //the pair for the return from brent's function
        std::pair<double,double> brentres;
        int bits = 24;



        //the location in the results vector for the current loop index
		auto vector_index = std::distance(t_start, t_start);

		for (auto it = t_start; it != t_end; ++it) {
			//see if we're not out-of-bounds with the end of the vector
			while (window_end != t_end - 1) {
				//if the time at the end of the window is less than current time + after, add the corresponding value sorted into the values vector
				if (*(window_end + 1) <= *it + after) {
					++window_end;
                    current_val = *(v_start + std::distance(t_start, window_end));
                    vals.insert(std::upper_bound(vals.begin(),
						vals.end(),
						current_val),
						current_val);
					
				}
				else break;

			}

			//see if it's time to increment the beginning of the window
			while (*window_start < *it - before) {
                //if so, deleate one value equal to the one that just went out of the window
                vals.erase(std::find(vals.begin(), vals.end(), *(v_start + std::distance(t_start, window_start))));
				++window_start;
			}

			//the value at the current time point
			test_value = *(v_start + std::distance(t_start, it));
            
            //calculate the size of the current window
			window_size = vals.size();
            //resize the vectors of transformed values
            if (v.size() != window_size) v.resize(window_size);
            if (logv.size() != window_size) logv.resize(window_size);

            //log-transform and sum the log values
            //this is needed even if lambda == 0 is never tested in Brent's algorithm, so it is done before
            std::transform(vals.begin(), vals.end(), logv.begin(), [](double &x){return log(x);});
            logsum = std::accumulate(logv.begin(), logv.end(), 0.0);

            //The function call to brent's algorithm
            //this will be used to find the best lambda, for the later transformations.
            //the lambda function is the boxcox likelyhood function as adapted from Python's implementation.

            brentres = boost::math::tools::brent_find_minima([&vals, &v, &logv, &logsum, &variance, &window_size](double lambda) -> double 
                {
                    //box-cox-transformation
                    if (lambda == 0) {
                        variance = var(logv);
                    } else {
                        //box-cox transform the values into v2
                        std::transform(vals.begin(),
                                        vals.end(),
                                        v.begin(), [&lambda](double& x){return pow(x, lambda)/lambda;});
                        variance = var(v);
                    }
                    return -1.0 * ((lambda - 1.0) * logsum  - (window_size/2) * log(variance));
                }, 
                min_lambda, 
                max_lambda, 
                bits);
                                                        
            bestlambda = std::get<0>(brentres);

            //boxcox-transform the values according to the best known lambda and calculate quantiles and IQR
            std::transform(vals.begin(), vals.end(), v.begin(), [&bestlambda](double &x){
            if (bestlambda == 0){
                return log(x);
            } else {
                return (pow(x, bestlambda)-1)/bestlambda;
            }});
            
            //calculate the quantile values in boxcox-transformed space
            std::transform(probs.begin(), 
                        probs.end(), 
                        qs.begin(), 
                        [&window_size, 
                            &index, 
                            &low, 
                            &high, 
                            &interp, 
                            &v](double &p){
                            index = (window_size-1)*p;
                            low = std::floor(index);
                            high = std::ceil(index);
                            interp = index - low;
                            return v[low]*(1-interp)+v[high]*interp;});

            //calculate the iqr*k in box cox transformed space
            iqr = (qs[1] - qs[0])*k;
            //change quantiles values to bounds in box cox transformed space
            qs[0] = qs[0]-iqr;
            qs[1] = qs[1]+iqr;

            //undo the box-cox transformation to get bounds in original space
            std::transform(qs.begin(), qs.end(), qs.begin(),[&bestlambda](double &x){
            if (bestlambda == 0){
                return exp(x);
            } else {
                return exp(log(x*bestlambda+1)/bestlambda);
            }
            });

            //recalculate the vector index
            vector_index = std::distance(t_start, it);
            //Tukey's fences - replace values outside of {q25 - k*IQR, q75 + k*IQR}
            if (test_value > qs[1]) *(r_start + vector_index) = qs[1]; else
                if (test_value < qs[0]) *(r_start + vector_index) = qs[0]; else
                    *(r_start + vector_index) = test_value;
            
            
        }
        return res;
	}
};

//the sliding hampel algorithm
//note, it has been tried with a std::multiset instead of std::vector to keep the sorted values, but
//vectors are much faster even with inserting in the middle.
class slide_hampel{
private:
	//the vector pointers
    Rcpp::NumericVector *values;
    Rcpp::NumericVector *times;
	
	//pointer to the start of our spans in the vectors
    Rcpp::NumericVector::iterator v_start;
    Rcpp::NumericVector::iterator t_start;
    //size of the spans in the vectors
	std::size_t size;

    //the return vector
    std::vector<double> res;
    //the start of our return vector
	std::vector<double>::iterator r_start;

	//sliding time range
	double before;
	double after;

	//the constant
	double k;


public:
	slide_hampel(Rcpp::NumericVector *v, Rcpp::NumericVector *t, int start, std::size_t group_size, double before_time, double after_time, double K) {
		//initialize our vectors
		values = v;
		times = t;
		res.resize(group_size);
		v_start = values->begin()+start;
		t_start = times->begin()+start;
		r_start = res.begin();

        //initialize the constants
		size = group_size;
		before = before_time;
		after = after_time;
		k = K;
	}

	auto operator ()() -> std::vector<double> {

		//The iterators for the sliding time window
		auto window_start = t_start;
		auto window_end = t_start;

		//the end of the times vector
		auto t_end = t_start + size;

		//The vector of values within the time range
		std::vector<double> vals;
        //reserve some memory and push back the first value
        vals.reserve(64);
        vals.push_back(*v_start);

        //the value to compare when increasing window size
        double current_val;

        //initialize the values for the outlier algorithm
        //the vector of absolute deviations
        std::vector<double> dev;
        dev.reserve(64);

        double median, mad;
        //the median element
        std::size_t m;
        //the value at our current loop index
        double test_value;

        //the location in the results vector
		auto vector_index = std::distance(t_start, t_start);

		for (auto it = t_start; it != t_end; ++it) {
			//see if we're not out-of-bounds with the end of the vector
			while (window_end != t_end - 1) {
				//if the time at the end of the window is less than current time + after, add the corresponding value sorted into the values vector
				if (*(window_end + 1) <= *it + after) {
					++window_end;
                    current_val = *(v_start + std::distance(t_start, window_end));
                    vals.insert(std::upper_bound(vals.begin(),
						vals.end(),
						current_val),
						current_val);
					
				}
				else break;

			}

			//see if it's time to increment the beginning of the window
			while (*window_start < *it - before) {
                vals.erase(std::find(vals.begin(), vals.end(), *(v_start + std::distance(t_start, window_start))));
				++window_start;
			}

			//the value at the current time point
			test_value = *(v_start + std::distance(t_start, it));

            //since the values are sorted, median is trivial
            m = vals.size() / 2;
            if (vals.size() % 2) median = vals[m]; 
				else median = (vals[m] + *std::max_element(vals.begin(), vals.begin() + m)) / 2.0;

            //calculate the absolue deviations
            if(dev.size() != vals.size()) dev.resize(vals.size());
            std::transform(vals.begin(), vals.end(), dev.begin(), [&median](double const &v){return std::fabs(v-median);});
            //calculate the mad
            std::nth_element(dev.begin(), dev.begin()+m, dev.end());
            if (dev.size() % 2) mad =  dev[m]; 
				else mad = (dev[m] + *std::max_element(dev.begin(), dev.begin() + m)) / 2.0;
            //scale the MAD by k
            mad = mad *k;

            //find the right position in the results-vector
			vector_index = std::distance(t_start, it);
            //the hampel filter
			if (test_value > median + mad) *(r_start + vector_index) = median; else
				if (test_value < median - mad) *(r_start + vector_index) = median; else
					*(r_start + vector_index) = test_value;

		}
		return res;
	}
};

//sliding algorithm that will quantify how many MADs away from the median the current value is
// - a variant on the hampel filter...
class slide_mads{
private:
	//the vector pointers
    Rcpp::NumericVector *values;
    Rcpp::NumericVector *times;
	
	//pointer to the start of our spans in the vectors
    Rcpp::NumericVector::iterator v_start;
    Rcpp::NumericVector::iterator t_start;
    //size of the spans in the vectors
	std::size_t size;

    //the return vector
    std::vector<double> res;
    //the start of our return vector
	std::vector<double>::iterator r_start;

	//sliding time range
	double before;
	double after;


public:
	slide_mads(Rcpp::NumericVector *v, Rcpp::NumericVector *t, int start, std::size_t group_size, double before_time, double after_time) {
		//initialize our vectors
		values = v;
		times = t;
		res.resize(group_size);
		v_start = values->begin()+start;
		t_start = times->begin()+start;
		r_start = res.begin();

        //initialize the constants
		size = group_size;
		before = before_time;
		after = after_time;
	}

	auto operator ()() -> std::vector<double> {

		//The iterators for the sliding time window
		auto window_start = t_start;
		auto window_end = t_start;

		//the end of the times vector
		auto t_end = t_start + size;

		//The vector of values within the time range
		std::vector<double> vals;
        //reserve some memory and push back the first value
        vals.reserve(64);
        vals.push_back(*v_start);

        //the value to compare when increasing window size
        double current_val;

        //initialize the values for the outlier algorithm
        //the vector of absolute deviations
        std::vector<double> dev;
        dev.reserve(64);

        double median, mad;
        //the median element
        std::size_t m;
        //the value at our current loop index
        double test_value;

        //the location in the results vector
		auto vector_index = std::distance(t_start, t_start);

		for (auto it = t_start; it != t_end; ++it) {
			//see if we're not out-of-bounds with the end of the vector
			while (window_end != t_end - 1) {
				//if the time at the end of the window is less than current time + after, add the corresponding value sorted into the values vector
				if (*(window_end + 1) <= *it + after) {
					++window_end;
                    current_val = *(v_start + std::distance(t_start, window_end));
                    vals.insert(std::upper_bound(vals.begin(),
						vals.end(),
						current_val),
						current_val);
					
				}
				else break;

			}

			//see if it's time to increment the beginning of the window
			while (*window_start < *it - before) {
                vals.erase(std::find(vals.begin(), vals.end(), *(v_start + std::distance(t_start, window_start))));
				++window_start;
			}

			//the value at the current time point
			test_value = *(v_start + std::distance(t_start, it));

            //since the values are sorted, median is trivial
            m = vals.size() / 2;
            if (vals.size() % 2) median = vals[m]; 
				else median = (vals[m] + *std::max_element(vals.begin(), vals.begin() + m)) / 2.0;

            //calculate the absolue deviations
            if(dev.size() != vals.size()) dev.resize(vals.size());
            std::transform(vals.begin(), vals.end(), dev.begin(), [&median](double const &v){return std::fabs(v-median);});
            //calculate the mad
            std::nth_element(dev.begin(), dev.begin()+m, dev.end());
            if (dev.size() % 2) mad =  dev[m]; 
				else mad = (dev[m] + *std::max_element(dev.begin(), dev.begin() + m)) / 2.0;

            //find the right position in the results-vector
			vector_index = std::distance(t_start, it);

			//not absolut values, as the sign of this might be useful at some point
			if (mad == 0){
				mad = 1.0;
			}
			*(r_start + vector_index) = (test_value - median) / mad;

		}
		return res;
	}
};

//2021-04-15
//This is the caller function for the dynamic thread pool, i.e. a pool that automatically resizes based on how many tasks
//are waiting.
//The pool has two controlling arguments - threshhold and multiplier. It will spawn a new thread if there are less than
//max threads, and
//number of tasks waiting > Threshhold + number_of_current_threads * multiplier.
//This keeps the pool from growing too large if the task is very small, but still allows it to leverage all possible threads if
//the task list is large enough.
//[[Rcpp::export]]
Rcpp::NumericVector tp_outliers(Rcpp::IntegerVector groups, 
                                      Rcpp::NumericVector times, 
                                      Rcpp::NumericVector values, 
                                      double before = 0, 
                                      double after = 0, 
                                      double k = 1.5, 
                                      std::string type = "tukey", 
                                      double minlambda = -2.0, 
                                      double maxlambda = 2.0,
									  double vtratio = 1.0,
									  int threadThreshold = 1,
									  int threadMultiplier = 4,
                                      int nThreads = 0,
									  int minThreads = 1,
									  int debug = 0,
                                      bool verbose = false)
{
    //iterators to keep track of the locations in the vectors where each group begins
	auto fstarttime = std::chrono::high_resolution_clock::now();
    //For the grouping algorithm
    int this_group{groups[0]};
    int g_start{0};
    std::size_t group_size;
    std::size_t n = groups.size();
	std::size_t pool_size;

	//some dynamic thread pool calculations
	if (nThreads <= 0){
		nThreads = std::max(1u,std::thread::hardware_concurrency());
	}
	
    //Initialize the vector of futures
    std::vector<std::future<std::vector<double>>> resvec;
    //And the actual output vector
    Rcpp::NumericVector out(n);
    auto r_start = std::begin(out);
    std::size_t offset{0};
    //the scope for the thread pool
    auto starttime = std::chrono::high_resolution_clock::now();
    {
        threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);

        if(verbose){
            Rcpp::Rcout << "Threadpool created, nThreads = " << minThreads << std::endl;
        }

        for (auto i = 0u;i<n;++i){
            //we've found a new group, or the end of the vector
            if ((groups[i] != this_group) | (i == n-1)){
                //Calcuate the size of the group
                if(i == n-1){
                    group_size = i - g_start + 1;
                } else {
                    group_size = i - g_start;
                }

                if (type == "tukey"){
                    resvec.push_back(std::move(pool.add_task(slide_tukey(&values, &times, g_start, group_size, before, after, k))));
                }

				if (type == "hampel"){
                    resvec.push_back(std::move(pool.add_task(slide_hampel(&values, &times, g_start, group_size, before, after, k))));
                }

				if (type == "mads"){
                    resvec.push_back(std::move(pool.add_task(slide_mads(&values, &times, g_start, group_size, before, after))));
                }

                if (type == "horn"){
                    resvec.push_back(std::move(pool.add_task(slide_horn(&values, &times, g_start, group_size, before, after, k, minlambda, maxlambda))));
                }

                this_group = groups[i];
                g_start = i;
                //t_start += group_size;
                //v_start += group_size;
            }
        }

		pool_size = pool.get_size();

    }
    auto endtime = std::chrono::high_resolution_clock::now();
    
    //Duration notification
    if(verbose){
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime);
        Rcpp::Rcout << "Threadpool finished, nThreads = " << pool_size << ", time elapsed = " << duration.count() << " ms." << std::endl;
        starttime = std::chrono::high_resolution_clock::now();
    }
    
    //copy the values back to the results vector
    for (auto &r : resvec){
        auto vec = r.get();
        std::copy_n(vec.begin(), vec.size(), r_start + offset);
        offset+=vec.size();
    }
    
    //Duration notification
    if(verbose){
        endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-fstarttime);
        Rcpp::Rcout << "All tasks finished, total time elapsed = " << duration.count() << " ms." << std::endl;
    }

    return out;
}