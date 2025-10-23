#include <Rcpp.h>
//[[Rcpp::plugins("cpp14")]]

#include <algorithm> //for nth_element
#include <vector>
#include <numeric> //for accumulate
//#include <random>

#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <future>
#include <queue>


//#include <string>
//#include <chrono>


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
		if ((waiting_tasks > (threadThreshhold + threadMultiplier * threads.size())) & (threads.size() < maxThreads)) {
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

	void add_thread() {
		//Rcpp::Rcout << "adding a thread...\n";
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

		for (auto& t : threads) {
			t.join();
		}
	}
};


//just mean
template<typename Tt, typename Tv>
double uniform_na_zero(Tt t_begin, Tt t_end, Tv v_begin, Tv v_end, double time_target, double time_coefficient) {
    if (t_begin == t_end) return 0;
    std::size_t n = std::distance(t_begin, t_end);

    double mean = std::accumulate(v_begin, v_end, 0.0) / n;
    return mean;
}

//mean, with triangular weights
template<typename Tt, typename Tv>
double triangle_na_zero(Tt t_begin, Tt t_end, Tv v_begin, Tv v_end, double time_target, double time_coefficient) {
    if (t_begin == t_end) return 0;
    std::size_t n = std::distance(t_begin, t_end);
    
    //store the sum of the weights
    double weight_sum {0.0};
    std::vector<double> weights(n);
    for (auto i = 0u;i<n;++i){
        weights[i] = 1.0 - std::fabs(*(t_begin+i) - time_target) * time_coefficient;
        weight_sum+=weights[i];
    }

    weight_sum = 1.0/weight_sum;
    double res = std::inner_product(v_begin, v_end, weights.begin(), 0.0)*weight_sum;
    return res;
}

//mean, with epanechnikov weights
template<typename Tt, typename Tv>
double epanechnikov_na_zero(Tt t_begin, Tt t_end, Tv v_begin, Tv v_end, double time_target, double time_coefficient) {
    if (t_begin == t_end) return 0;
    std::size_t n = std::distance(t_begin, t_end);
    
    //store the sum of the weights
    double weight_sum {0.0};
    std::vector<double> weights(n);
    //double dt{0.0};
    for (auto i = 0u;i<n;++i){
        //weights[i] = (3.0*(1.0 - std::pow((*(t_begin+i) - time_target()) * time_coefficient, 2.0)))/4.0;
       //we're gonna scale the sum to 1.0 anyway so skip the 3/4-part
       weights[i] = 1.0 - std::pow((*(t_begin+i) - time_target) * time_coefficient, 2.0);
       weight_sum+=weights[i];
    }

    weight_sum = 1.0/weight_sum;
    double res = std::inner_product(v_begin, v_end, weights.begin(), 0.0)*weight_sum;
    return res;
}

//mean, with tricubic weights
template<typename Tt, typename Tv>
double tricubic_na_zero(Tt t_begin, Tt t_end, Tv v_begin, Tv v_end, double time_target, double time_coefficient) {
    if (t_begin == t_end) return 0;
    std::size_t n = std::distance(t_begin, t_end);
    
    //store the sum of the weights
    double weight_sum {0.0};
    std::vector<double> weights(n);
    //double dt{0.0};
    for (auto i = 0u;i<n;++i){

       // weights[i] = 70.0*std::pow(1.0-std::pow(std::fabs(*(t_begin+i) - time_target()) * time_coefficient,3.0), 3.0)/81.0;
       //we're gonna scale the sum to 1.0 anyway so skip the 70/81-part
       weights[i] = std::pow(1.0-std::pow(std::fabs(*(t_begin+i) - time_target) * time_coefficient,3.0), 3.0);
       weight_sum+=weights[i];
    }

    weight_sum = 1.0/weight_sum;
    double res = std::inner_product(v_begin, v_end, weights.begin(), 0.0)*weight_sum;
    return res;
}

//class for summarizing means using some weighting function
template < typename Td, typename Tg, typename Tv >
class timeweight_summarize {
private:
    //the data vectors
    Tg* groups;
    Td* times;
    Tv* values;

    //the intervals data vectors
    Tg* interval_groups;
    Td* interval_times;

    //the results vector
    Tv* res;

    //our offsets in the data vectors
    std::size_t group_start;
    std::size_t group_end;

    //our offsets in the intervals vector
    std::size_t interval_group_start = 0;
    std::size_t interval_group_end = 0;

    //defining the time intervals
    double interval_time_before;
    double interval_time_after;

    //summary function
    int fun_aggregate = 0;

public:

    timeweight_summarize(Tg* groups, 
    Td* times, 
    Tv* values, 
    Tg* interval_groups, 
    Td* interval_times, 
    Tv* res, 
    std::size_t group_start, 
    std::size_t group_size, 
    std::size_t interval_group_start,
    std::size_t interval_group_end,
    double interval_time_before = 120,
    double interval_time_after = 120,
    int fun_aggregate = 0
    ) {
        //copy the data pointers
        this->groups = groups;
        this->times = times;
        this->values = values;

        this->interval_groups = interval_groups;
        this->interval_times = interval_times;
        this->res = res;
        //and our locations in the vectors
        this->group_start = group_start;
        this->group_end = group_start + group_size;   
        this->interval_group_start = interval_group_start;
        this->interval_group_end = interval_group_end;

        //and the distances to search before and after
        this->interval_time_before = interval_time_before;
        this->interval_time_after = interval_time_after;

        //fetch what sort of aggregation function is wanted
        this->fun_aggregate = fun_aggregate;


    }

    auto operator()() -> void {
        //this required some trial and error.....
        //use a typedef to the iterator to make it less wordy...
        typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::iterator VALUES_ITERATOR_TYPE;
        typedef typename decltype(typename std::remove_pointer < decltype(times) > ::type())::iterator TIMES_ITERATOR_TYPE;
        //function pointer to some sort of aggregation function working on iterators
        double (*f)(VALUES_ITERATOR_TYPE, VALUES_ITERATOR_TYPE, TIMES_ITERATOR_TYPE, TIMES_ITERATOR_TYPE, double, double);

        //store the iterators
        auto vbegin = values->begin();
        auto tbegin = times->begin();
        //also store the time coefficient that will be required to shift delta-t's to the range 0...1
        //notice that the kernel may be truncated, if requested
        double time_coefficient= 1.0/std::max(std::fabs(interval_time_before), std::fabs(interval_time_after));

        switch(fun_aggregate){
            case 0: f = &uniform_na_zero;
                    break;
            case 1: f = &triangle_na_zero;
                    break;
            case 2: f = &epanechnikov_na_zero;
                    break;
            case 3: f = &tricubic_na_zero;
                    break;
            default: f = &uniform_na_zero;
                    break;
        }



        std::size_t aggregate_begin = group_start;
        std::size_t aggregate_end = group_start;
        //std::size_t best_match = group_start;
        //bool matching_found = false;
        //double best_difference = 0;

        //loop through the intervals
        for (auto i = interval_group_start;i<interval_group_end;++i){
            //matching_found = false;
            //move the end index as far forward as possible

            while (((*times)[aggregate_end]<=((*interval_times)[i]+interval_time_after)) && ((aggregate_end<group_end))){
                ++aggregate_end;
            }
            //move the start index as far as possible
            while (((*times)[aggregate_begin]<((*interval_times)[i]-interval_time_before)) && ((aggregate_begin<group_end))){

                ++aggregate_begin;
            }

            //if the indices are different, call the summary function
            if (aggregate_end != aggregate_begin){
                (*res)[i] = (*f)(tbegin+aggregate_begin, tbegin+aggregate_end, vbegin+aggregate_begin, vbegin+aggregate_end, (*interval_times)[i], time_coefficient);
            }
        }
    }

};


//timeweighted summarization with various kernel functions
//[[Rcpp::export]]
Rcpp::NumericVector cpp_timeweight_summarize(Rcpp::NumericVector groups,
    Rcpp::NumericVector times,
    Rcpp::NumericVector values,
    Rcpp::NumericVector int_groups,
    Rcpp::NumericVector int_times,
    Rcpp::String aggregation_function = "uniform",
    double interval_time_before = 300,
    double interval_time_after = 300,
    int nThreads = 0,
    int threadThreshold = 32,
    int threadMultiplier = 64,
    int minThreads = 1){
//setup aggregation function
int fun_aggregate = 0;
if (aggregation_function == "uniform"){
fun_aggregate = 0;
} else if (aggregation_function == "triangle"){
fun_aggregate = 1;
} else if (aggregation_function == "epanechnikov"){
fun_aggregate = 2;
} else if (aggregation_function == "tricubic"){
fun_aggregate = 3;
} else {
Rcpp::stop("'aggregation_function' must be one of 'uniform', 'triangle', 'epanechnikov', or 'tricubic'");
}

if (interval_time_after < (interval_time_before * -1.0)){
Rcpp::stop("'interval_time_after' must be greater or equal to 'interval_time_before'*-1");
}

//some dynamic thread pool calculations
if (nThreads <= 0){
nThreads = std::max(1u,std::thread::hardware_concurrency());
//Rcpp::Rcout << "nThreads: " << nThreads << "\n";
}


std::size_t n = groups.size();
std::size_t n_int = int_groups.size();
//the out vectors, and the vector of futures
//Rcpp::NumericVector res_time(n_int, Rcpp::NumericVector::get_na());
Rcpp::NumericVector res(n_int, Rcpp::NumericVector::get_na());
std::vector<std::future<void>> resvec;

//the location of the groups
std::size_t this_group = groups[0];
std::size_t this_group_start = 0u;
std::size_t this_group_size = 0u;

//the location of the groups within the intervals
std::size_t interval_group_start = 0u;
std::size_t interval_group_end = 0u;

//a flag to check if the current ID has been found within the intervals
bool in_intervals = false;

//the scope for the thread pool
{
threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
//A loop through the first grouping vector
for (auto i = 0u; i < n; ++i) {
//set the size of the group
if ((groups[i] != this_group) | (i == n - 1)) {
if (i == n - 1) {
this_group_size = i - this_group_start + 1;
}
else {
this_group_size = i - this_group_start;
}

//try to find the id in the interval groups
in_intervals = false;

for (auto j = interval_group_end;j<n_int;++j){
if ((!in_intervals) && (int_groups[j] == this_group)){
in_intervals = true;
interval_group_start = j;
interval_group_end = n_int;
}

if (in_intervals && (int_groups[j] != this_group )){
interval_group_end = j;
break;
}
}



//std::cout << "group found! \n";
//std::cout << "this_group_start: " << this_group_start << " this group size: " << this_group_size << "\n";
if (in_intervals){
resvec.push_back(std::move(pool.add_task(


timeweight_summarize < Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
&times,
&values,
&int_groups,
&int_times,
//&int_values, 
&res,
this_group_start,
this_group_size,
interval_group_start,
interval_group_end,
interval_time_before,
interval_time_after,
fun_aggregate
)

)));
}


//resvec.push_back(std::move(pool.add_task(cof_dm<std::vector<double>, std::vector<double>::iterator>(data, coeffs, this_group_start, this_group_size, 4, 1, false))));

//interval_match<std::vector<double>,std::vector<int>>(&groups, &times, &int_groups, &int_begin, &int_end, &res, this_group_start, this_group_size)
this_group = groups[i];
this_group_start = i;
}
}
//Rcpp::Rcout << "threadpool finished, nthreads = " << pool.get_size() << "\n";
}

for (auto &r : resvec){
r.get();
}

return res;
}