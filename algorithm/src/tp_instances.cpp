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

//class for aggregating values over intervals, this time allowing overlapping intervals!
template < typename Td, typename Tg>
class interval_overlap_sum {
private:
    //the data vectors
    Tg* groups;
    Td* times_begin;
    Td* times_end;

    //the intervals data vectors
    Tg* interval_groups;
    Td* interval_begin;
    Td* interval_end;

    //the results vector
    Td* res;

    //our offsets in the data vectors
    std::size_t group_start;
    std::size_t group_end;

    //our offsets in the intervals vector
    std::size_t interval_group_start = 0;
    std::size_t interval_group_end = 0;

    bool include_everything_before;


public:

    interval_overlap_sum(Tg* groups, 
                        Td* times_begin, 
                        Td* times_end, 
                        Tg* interval_groups, 
                        Td* interval_begin, 
                        Td* interval_end, 
                        Td* res, 
                        std::size_t group_start, 
                        std::size_t group_size, 
                        std::size_t interval_group_start, 
                        std::size_t interval_group_end,
                        bool include_everything_before) {
        //copy the data pointers
        this->groups = groups;
        this->times_begin = times_begin;
        this->times_end = times_end;

        this->interval_groups = interval_groups;
        this->interval_begin = interval_begin;
        this->interval_end = interval_end;

        this->res = res;

        this->group_start = group_start;
        this->group_end = group_start + group_size;
        this->interval_group_start = interval_group_start;
        this->interval_group_end = interval_group_end;

        this->include_everything_before = include_everything_before;

    }

    auto operator()() -> void {
        //to keep the number of loops down, and for locf-purposes, store the index of the last element before the start of the current interval
        //std::size_t last_element_before = group_start;
        double overlap_begin;
        double overlap_end;

        //Rcpp::Rcout << "thread = " << std::this_thread::get_id() << ", group = " << (*groups)[group_start] << ", group_start : " << group_start << ", group_end: " << group_end << ", interval_group: " << (*interval_groups)[interval_group_start] << ", interval_group_start: " << interval_group_start << ", interval_group_end: " << interval_group_end << "\n";
        bool new_interval = false;
        //i loops through the intervals
        for (auto i = interval_group_start; i < interval_group_end; ++i){
            //Rcpp::Rcout << "i: " << i << "\n";
            //initialize the sum as 0.0
            (*res)[i] = 0.0;
            new_interval = true;
            //j loops through the groups
            for (auto j = group_start;j<group_end;++j){
                //Rcpp::Rcout << "j: " << j << "\n";
                // we have an overlap
                if ((*times_begin)[j] < (*interval_end)[i] && (*times_end)[j] > (*interval_begin)[i]){
                    
                    //the beginning of the overlap is whichever is greater of interval_begin and times_begin
                    if (include_everything_before){
                        overlap_begin = (*times_begin)[j];
                    } else {
                        if ((*interval_begin)[i] > (*times_begin)[j]){
                            overlap_begin = (*interval_begin)[i];
                        } else {
                            overlap_begin = (*times_begin)[j];
                        }
                    }

                    

                    //the end of the overlap is whichever is smaller of interval_end and times_end
                    if ((*interval_end)[i] < (*times_end)[j]){
                        overlap_end = (*interval_end)[i];
                    } else {
                        overlap_end = (*times_end)[j];
                    }

                    //Rcpp::Rcout << "Overlap between interval " << i << " and group " << j << " at time (" << (*times_begin)[j] << ", " <<(*times_end)[j] << ") and interval ("<<  (*interval_begin)[i] << ", " << (*interval_end)[i] << ")\n";
                    //add the overlap duration to the result

                    //Rcpp::Rcout << "Adding " << overlap_end - overlap_begin << " to res[" << i << "]\n";
                    (*res)[i] += (overlap_end - overlap_begin);
                }

                //if the end time for this times-interval is before the start of the current interval,
                //we won't have to check it again
                if (new_interval && (*times_end)[j] < (*interval_begin)[i]){
                    group_start++;
                } else {
                    new_interval = false;
                }

                //break the inner loop if the start of the times is greater than the end of the intervals
                if ((*times_begin)[j] > (*interval_end)[i]) break;

            }
        }
    }
};


//a specialized and hopefully fast way of summing the overlap time between two sets of intervals
//[[Rcpp::export]]
Rcpp::NumericVector cpp_interval_overlap_sum(Rcpp::NumericVector groups,
                                        Rcpp::NumericVector times_begin,
                                        Rcpp::NumericVector times_end,
                                        Rcpp::NumericVector int_groups,
                                        Rcpp::NumericVector int_begin,
                                        Rcpp::NumericVector int_end,
                                        bool include_everything_before = false,
                                        int nThreads = 0,
                                        int threadThreshold = 32,
                                        int threadMultiplier = 64,
                                        int minThreads = 1){
 
    //some dynamic thread pool calculations
    if (nThreads <= 0){
        nThreads = std::max(1u,std::thread::hardware_concurrency());
        //Rcpp::Rcout << "nThreads: " << nThreads << "\n";
    }
    

    std::size_t n = groups.size();
	std::size_t n_int = int_groups.size();
    //the out vector, and the vector of futures
    Rcpp::NumericVector res(n_int, Rcpp::NumericVector::get_na());
    std::vector<std::future<void>> resvec;

    //the location of the groups
    std::size_t this_group = groups[0];
	std::size_t this_group_start = 0u;
	std::size_t this_group_size;

    //the location of the groups within the intervals
    std::size_t interval_group_start = 0u;
    std::size_t interval_group_end = 0u;

    //a flag to check if the current ID has been found within the intervals
    bool in_intervals = false;

    //the scope for the thread pool
    {
        threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
        //Rcpp::Rcout << "Threadpool started...!\n";
        for (auto i = 0u; i < n; ++i) {
			//set the size of the group
			if ((groups[i] != this_group) | (i == n - 1u)) {
				if (i == n - 1u) {
					this_group_size = i - this_group_start + 1;
				}
				else {
					this_group_size = i - this_group_start;
				}

                //new - check the intervals in the main function and pass locations to the tasks - more efficient for very large number of intervals
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

                    interval_overlap_sum < Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
                        &times_begin,
                        &times_end,
                        &int_groups,
                        &int_begin,
                        &int_end, 
                        &res,
                        this_group_start,
                        this_group_size,
                        interval_group_start,
                        interval_group_end,
                        include_everything_before)

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

//inline function to set *a to the min of a and b
inline void min_inplace(double &a, const double &b){
    if (b < a) a = b;
}
//inline function to set a to the max of a and b
inline void max_inplace(double &a, const double &b){
    if (b > a) a = b;
}
//inline function to set a to a+b
inline void sum_inplace(double &a, const double &b){
    a+=b;
}

//class for aggregating values over intervals, this time allowing overlapping intervals!
template < typename Td, typename Tg>
class interval_overlap_summarize {
private:
    //the data vectors
    Tg* groups;
    Td* times_begin;
    Td* times_end;
    Td* values;

    //the intervals data vectors
    Tg* interval_groups;
    Td* interval_begin;
    Td* interval_end;

    //the results vector
    Td* res;

    //our offsets in the data vectors
    std::size_t group_start;
    std::size_t group_end;

    //our offsets in the intervals vector
    std::size_t interval_group_start = 0;
    std::size_t interval_group_end = 0;

    int fun_aggregate;


public:

    interval_overlap_summarize(Tg* groups, 
                        Td* times_begin, 
                        Td* times_end, 
                        Td* values,
                        Tg* interval_groups, 
                        Td* interval_begin, 
                        Td* interval_end, 
                        Td* res, 
                        std::size_t group_start, 
                        std::size_t group_size, 
                        std::size_t interval_group_start, 
                        std::size_t interval_group_end,
                        int fun_aggregate) {
        //copy the data pointers
        this->groups = groups;
        this->times_begin = times_begin;
        this->times_end = times_end;
        this->values = values;


        this->interval_groups = interval_groups;
        this->interval_begin = interval_begin;
        this->interval_end = interval_end;

        this->res = res;

        this->group_start = group_start;
        this->group_end = group_start + group_size;
        this->interval_group_start = interval_group_start;
        this->interval_group_end = interval_group_end;

        this->fun_aggregate = fun_aggregate;

    }

    auto operator()() -> void {
        //set up the function pointer to the aggregation function

        //this required some trial and error.....
        //use a typedef to the iterator to make it less wordy...
        //typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::iterator ITERATOR_TYPE;
        //typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::value_type VALUES;
        //function pointer to some sort of aggregation function working on iterators
        void (*f)(double&, const double&);

        switch(fun_aggregate){
            case 0: f = &min_inplace;
                    break;
            case 1: f = &max_inplace;
                    break;
            case 2: f = &sum_inplace;
            default: f = &min_inplace;
                break;
        }


        //keep track of when to increase the start position for the inner loop
        bool new_interval;
        //keep track of whether to call the function or not
        bool overlap_found;
        //maybe better to keep this on the stack for performance - right?
        double this_res;
        //i loops through the intervals
        for (auto i = interval_group_start; i < interval_group_end; ++i){
            overlap_found = false;
            new_interval = true;
            //j loops through the groups
            for (auto j = group_start;j<group_end;++j){
                // we have an overlap
                if ((*times_begin)[j] < (*interval_end)[i] && (*times_end)[j] > (*interval_begin)[i]){
                    if (!overlap_found){
                        this_res = (*values)[j];
                    } else {
                        //call aggregation function here!
                        (*f)(this_res, (*values)[j]);
                    }


                }

                //if the end time for this times-interval is before the start of the current interval,
                //we won't have to check it again
                if (new_interval && (*times_end)[j] < (*interval_begin)[i]){
                    group_start++;
                } else {
                    new_interval = false;
                }

                //break the inner loop if the start of the times is greater than the end of the current intervals
                if ((*times_begin)[j] > (*interval_end)[i]) break;

            }

            //yes, this makes a copy but it's only a double...
            if (overlap_found) {
                (*res)[i] = this_res;
            }
        }
    }
};


//a specialized and hopefully fast method for summarizing fast overlap checks
//[[Rcpp::export]]
Rcpp::NumericVector cpp_interval_overlap_summarize(Rcpp::NumericVector groups,
                                        Rcpp::NumericVector times_begin,
                                        Rcpp::NumericVector times_end,
                                        Rcpp::NumericVector values,
                                        Rcpp::NumericVector int_groups,
                                        Rcpp::NumericVector int_begin,
                                        Rcpp::NumericVector int_end,
                                        Rcpp::String aggregation_function,
                                        int nThreads = 0,
                                        int threadThreshold = 32,
                                        int threadMultiplier = 64,
                                        int minThreads = 1){
 
    //some dynamic thread pool calculations
    if (nThreads <= 0){
        nThreads = std::max(1u,std::thread::hardware_concurrency());
        //Rcpp::Rcout << "nThreads: " << nThreads << "\n";
    }

    //setup aggregation function
    int fun_aggregate = 0;
    if (aggregation_function == "min"){
        fun_aggregate = 0;
    } else if (aggregation_function == "max"){
        fun_aggregate = 1;
    } else if (aggregation_function == "sum"){
        fun_aggregate = 2;
    }  else {
        Rcpp::stop("'aggregation_function' must be one of 'min', 'max', 'sum'");
    }
    

    std::size_t n = groups.size();
	std::size_t n_int = int_groups.size();
    //the out vector, and the vector of futures
    Rcpp::NumericVector res(n_int, Rcpp::NumericVector::get_na());
    std::vector<std::future<void>> resvec;

    //the location of the groups
    std::size_t this_group = groups[0];
	std::size_t this_group_start = 0u;
	std::size_t this_group_size;

    //the location of the groups within the intervals
    std::size_t interval_group_start = 0u;
    std::size_t interval_group_end = 0u;

    //a flag to check if the current ID has been found within the intervals
    bool in_intervals = false;

    //the scope for the thread pool
    {
        threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
        //Rcpp::Rcout << "Threadpool started...!\n";
        for (auto i = 0u; i < n; ++i) {
			//set the size of the group
			if ((groups[i] != this_group) | (i == n - 1)) {
				if (i == n - 1) {
					this_group_size = i - this_group_start + 1;
				}
				else {
					this_group_size = i - this_group_start;
				}

                //new - check the intervals in the main function and pass locations to the tasks - more efficient for very large number of intervals
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

                    interval_overlap_summarize < Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
                        &times_begin,
                        &times_end,
                        &values,
                        &int_groups,
                        &int_begin,
                        &int_end, 
                        &res,
                        this_group_start,
                        this_group_size,
                        interval_group_start,
                        interval_group_end,
                        fun_aggregate)

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