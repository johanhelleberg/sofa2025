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

//class for matching intervals, templated on time measurement and grouping variable
template <typename Td, typename Tg>
class interval_match {
private:
	//the data vectors
	Tg* groups;
	Td* times;
	
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
public:
	interval_match(Tg* groups, Td* times, Tg* interval_groups, Td* interval_begin, Td* interval_end, Td* res, std::size_t group_start, std::size_t group_size) {
		//copy the data pointers
		this->groups = groups;
		this->times = times;
		this->interval_groups = interval_groups;
		this->interval_begin = interval_begin;
		this->interval_end = interval_end;
		this->res = res;
		this->group_start = group_start;
		this->group_end = group_start+group_size;
	}

	auto operator()()->void {
		auto this_group = (*groups)[group_start];
		std::size_t intervals_total_size = interval_groups->size();

		//std::cout << "lei ho from group: " << this_group << "\n";
		//std::cout << "group_start: " << group_start << " group_end: " << group_end << "\n";

		//find our id in the intervals
		bool in_intervals = false;
		for (auto i = 0u; i < intervals_total_size; ++i) {
			if ((!in_intervals) && (*interval_groups)[i] == this_group) {
				interval_group_start = i;
				interval_group_end = intervals_total_size;
				in_intervals = true;
			}
			if (in_intervals && (*interval_groups)[i] != this_group) {
				interval_group_end = i;
				break;
			}
		}

		//if we're not in the intervals - just return without trying to match at all
		if (!in_intervals) {
			return;
		}
		//std::size_t nloops = 0;

		std::size_t last_interval_found = interval_group_start;
		//loop through all our points
		for (auto i = group_start; i < group_end; ++i) {
			//loop through the candidate intervals, beginning from the last interval that was matched - i.e. take advantage of the fact that 
			//the intervals and points are sorted on time
			for (auto j = last_interval_found; j < interval_group_end; ++j) {
				//++nloops;
				//semi-open-intervals, i.e. >= start_time but < end_time
				//before the first interval - just stop searching
				if ((*times)[i] < (*interval_begin)[j]) {
					break;
				}
				if ((*times)[i] >= (*interval_begin)[j] && (*times)[i] < (*interval_end)[j]) {
					last_interval_found = j;
					(*res)[i] = (*interval_begin)[j];
					break;
				}
			}

			
			//if (*times)[i] >= 
		}
		//std::cout << "nloops: " << nloops << "\n";
		//std::cout << "interval_group_start: " << interval_group_start << " interval_group_end: " << interval_group_end << "\n";
	}

};

//matching first() and last() in R, for container iterators
template<typename T>
T first_element(T begin, T end){
    return begin;
}

template <typename T>
T last_element(T begin, T end){
    if (begin == end){
        return begin;
    } else {
        return end - 1;
    }
}

//returning the type, not the iterator
template<typename Ti>
double first_na_zero(Ti begin, Ti end){
    if (begin == end){
        return 0;
    } else {
        return *begin;
    }
}

//returning the type, not the iterator
template<typename Ti>
double last_na_zero(Ti begin, Ti end){
    if (begin == end){
        return 0;
    } else {
        return *(end-1);
    }
}

//median, as double
template<typename Ti>
double median_na_zero( Ti begin,  Ti end) {
    if (begin == end) return 0;

    std::size_t n = std::distance(begin, end);
    std::size_t m = n / 2;

    std::vector<double> vec(n);

    std::copy(begin, end, vec.begin());
    std::nth_element(vec.begin(), vec.begin()+m, vec.end());
    if (n % 2) {
        return vec[m];
    }
    else {
        return (vec[m] + *std::max_element(vec.begin(), vec.begin() + m)) * 0.5;
    }
}

//mean, as double
template<typename Ti>
double mean_na_zero( Ti begin,  Ti end) {
    if (begin == end) return 0;
    std::size_t n = std::distance(begin, end);
    double mean = std::accumulate(begin, end, 0.0) / n;
    return mean;
}

//mean, as double
template<typename Ti>
double sum_na_zero( Ti begin,  Ti end) {
    if (begin == end) return 0;
    double sum = std::accumulate(begin, end, 0.0);
    return sum;
}

//variance in two passes
template<typename Ti>
double var_na_zero( Ti begin,  Ti end) {
    if (begin == end) return 0;
    std::size_t n = std::distance(begin, end);
    double mean = std::accumulate(begin, end, 0.0) / n;
    double var = std::accumulate(begin, end, 0.0, [&](const auto r, const auto v){return r+(v - mean) * (v - mean); }) / (n - 1);
    return var;
}

template<typename Ti>
double sd_na_zero( Ti begin,  Ti end){
    return std::sqrt(var_na_zero(begin, end));
}

//min
template<typename Ti>
double min_na_zero( Ti begin,  Ti end) {
    return *std::min_element(begin, end);
}

//max
template<typename Ti>
double max_na_zero( Ti begin,  Ti end) {
    return *std::max_element(begin, end);
}

//class for aggregating values over intervals, this time allowing overlapping intervals!
// template < typename Td, typename Tg, typename Tv >
// class foverlap_aggregate {
// private:
//     //the data vectors
//     Tg* groups;
//     Td* begin_times;
//     Td* end_times;
//     Tv* values;

//     //the intervals data vectors
//     Tg* interval_groups;
//     Td* interval_begin;
//     Td* interval_end;

//     //the results vector
//     Tv* res;
//     Td* res_begin;
//     Td* res_end;
//     Td* res_time_difference;

//     //our offsets in the data vectors
//     std::size_t group_start;
//     std::size_t group_end;

//     //our offsets in the intervals vector
//     std::size_t interval_group_start = 0;
//     std::size_t interval_group_end = 0;

//     //our aggregation function
//     int fun_aggregate = 0;

//     double locf;
//     double nocb;


// public:

//     foverlap_aggregate(Tg* groups, Td* begin_times, Td* end_times, Tv* values, Tg* interval_groups, Td* interval_begin, Td* interval_end, Tv* res, Td* res_begin, Td* res_end, std::size_t group_start, std::size_t group_size, std::size_t interval_group_start, std::size_t interval_group_end, int fun_aggregate, double locf = 0.0, double nocb = 0.0) {
//         //copy the data pointers
//         this->groups = groups;
//         this->begin_times = begin_times;
//         this->end_times = end_times;
//         this->values = values;

//         this->interval_groups = interval_groups;
//         this->interval_begin = interval_begin;
//         this->interval_end = interval_end;

//         this->interval_group_start = interval_group_start;
//         this->interval_group_end = interval_group_end;

//         this->res = res;
//         this->res_begin = res_begin;
//         this->res_end = res_end;
//         //this->res_time_difference = res_time;

//         this->group_start = group_start;
//         this->group_end = group_start + group_size;
//         this->fun_aggregate = fun_aggregate;
//         this->locf = locf;
//         this->nocb = nocb;

//     }

//     auto operator()() -> void {

//         auto this_group = (*groups)[group_start];

//         std::size_t intervals_total_size = interval_groups->size();
//         //std::cout << "lei ho from group: " << this_group << "\n";
//         //std::cout << "group_start: " << group_start << " group_end: " << group_end << "\n";

//         //find our id in the intervals
//         //bool in_intervals = false;

//         // for (auto i = 0u; i < intervals_total_size; ++i) {

//         //     if ((!in_intervals) && (*interval_groups)[i] == this_group) {

//         //         interval_group_start = i;

//         //         interval_group_end = intervals_total_size;

//         //         in_intervals = true;

//         //     }

//         //     if (in_intervals && (*interval_groups)[i] != this_group) {

//         //         interval_group_end = i;

//         //         break;

//         //     }

//         // }

//         // //std::cout << "interval_group_end: " << interval_group_end << "\n";

//         // //if we're not in the intervals - just return without trying to match at all
//         // if (!in_intervals) {

//         //     return;

//         // }

//         //just for diagnostic purposes, add a counter of the number of loops.
//         std::size_t nloops = 0;

//         //to keep the number of loops down, and for locf-purposes, store the index of the last element before the start of the current interval
//         std::size_t last_element_before = group_start;

//         //this required some trial and error.....
//         //use a typedef to the iterator to make it less wordy...
//         typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::iterator ITERATOR_TYPE;
//         //function pointer to some sort of aggregation function working on iterators
//         ITERATOR_TYPE(*f)(ITERATOR_TYPE, ITERATOR_TYPE);

//         typedef typename decltype(typename std::remove_pointer < decltype(times) > ::type())::value_type TIMES;
//         TIMES locf_max = static_cast<TIMES>(locf);
//         TIMES nocb_max = static_cast<TIMES>(nocb);
//         //static_cast<TIMES>
//         //std::cout << "locf_max: " << locf_max << "\n";

//         switch(fun_aggregate){
//             case 0: f = &std::min_element;
//                     break;
//             case 1: f = &std::max_element;
//                     break;
//             case 2: f = &first_element;
//                     break;
//             case 3: f = &last_element;
//                     break;
//             default: f = &std::min_element;
//                     break;
//         }

//         auto aggregate_begin = values->begin();
//         auto aggregate_end = values->begin();

//         //loop through the intervals
//         for (auto i = interval_group_start; i < interval_group_end; ++i) {
//             bool matching_elements_found = false;
//             bool result_set = false;
//             //loop through values
//             for (auto j = last_element_before; j < group_end; ++j) {
//                 ++nloops;
//                 //update the last element before the interval start, i.e. possible locf-value
//                 if ((*times)[j] < (*interval_begin)[i]) {
//                     last_element_before = j;
//                 }
//                 else if (!matching_elements_found && (*times)[j] >= (*interval_begin)[i] && (*times)[j] <= (*interval_end)[i]) {
//                     matching_elements_found = true;
//                     aggregate_begin = values->begin() + j;
                    
//                 }
//                 //break the loop if the current time is larger than the current interval end or we're at the end of the group
//                 if (((*times)[j] > (*interval_end)[i]) || j == group_end - 1) {
//                     //call aggregation function, if we've previously found elements
//                     if (matching_elements_found) {
//                         //the end iterator is either j, j+1, depending on the condition for break
//                         aggregate_end = values->begin() + j;

//                         //make sure we include the last value when it is also within aggregation
//                         if ((j == (group_end - 1)) && ((*times)[j] <= (*interval_end)[i])){
//                             ++aggregate_end;
//                         }

//                         //apply our aggregation function
//                         auto target_position = (*f)(aggregate_begin, aggregate_end);

//                         //save value and time from the correct index
//                         if (target_position != values->end()) {
//                             (*res)[i] = (*target_position);
//                             auto t_index = std::abs(std::distance(values->begin(), target_position));
//                             (*res_time)[i] = (*times)[t_index];
//                             //we are in the interval, so we miss by 0
//                             //(*res_time_difference)[i] = 0;
//                             result_set = true;
//                         }
//                     }

//                     //if we haven't set any result yet, possibly set result using locf
//                     if (!result_set && ((*times)[last_element_before] < (*interval_begin)[i]) && ((*times)[last_element_before] + locf_max >= (*interval_begin)[i])) {
//                         (*res)[i] = (*values)[last_element_before];
//                         (*res_time)[i] = (*times)[last_element_before];
//                         result_set = true;
//                     }

//                     //if that doesn't work, try nocb
//                     if (!result_set && (*times)[j] > (*interval_end)[i] && (*times)[j] - nocb_max <= (*interval_end)[i]) {
//                         (*res)[i] = (*values)[j];
//                         (*res_begin)[i] = (*times)[j];
//                         result_set = true;
//                     }

//                     //std::cout << "interval = " << i << " last_element_before = " << last_element_before << "\n";
//                     //finally, break the value loop
//                     break;
//                 }
//             }
//         }
//     }

// };

//class for aggregating values over intervals, this time allowing overlapping intervals!
template < typename Td, typename Tg, typename Tv >
class interval_aggregate_overlap {
private:
    //the data vectors
    Tg* groups;
    Td* times;
    Tv* values;

    //the intervals data vectors
    Tg* interval_groups;
    Td* interval_begin;
    Td* interval_end;

    //the results vector
    Tv* res;
    Td* res_time;
    Td* res_time_difference;

    //our offsets in the data vectors
    std::size_t group_start;
    std::size_t group_end;

    //our offsets in the intervals vector
    std::size_t interval_group_start = 0;
    std::size_t interval_group_end = 0;

    //our aggregation function
    int fun_aggregate = 0;

    double locf;
    double nocb;


public:

    interval_aggregate_overlap(Tg* groups, Td* times, Tv* values, Tg* interval_groups, Td* interval_begin, Td* interval_end, Tv* res, Td* res_time, std::size_t group_start, std::size_t group_size, std::size_t interval_group_start, std::size_t interval_group_end, int fun_aggregate, double locf = 0.0, double nocb = 0.0) {
        //copy the data pointers
        this->groups = groups;
        this->times = times;
        this->values = values;

        this->interval_groups = interval_groups;
        this->interval_begin = interval_begin;
        this->interval_end = interval_end;

        this->interval_group_start = interval_group_start;
        this->interval_group_end = interval_group_end;

        this->res = res;
        this->res_time = res_time;
        //this->res_time_difference = res_time;

        this->group_start = group_start;
        this->group_end = group_start + group_size;
        this->fun_aggregate = fun_aggregate;
        this->locf = locf;
        this->nocb = nocb;

    }

    auto operator()() -> void {

        auto this_group = (*groups)[group_start];

        std::size_t intervals_total_size = interval_groups->size();
        //std::cout << "lei ho from group: " << this_group << "\n";
        //std::cout << "group_start: " << group_start << " group_end: " << group_end << "\n";

        //find our id in the intervals
        //bool in_intervals = false;

        // for (auto i = 0u; i < intervals_total_size; ++i) {

        //     if ((!in_intervals) && (*interval_groups)[i] == this_group) {

        //         interval_group_start = i;

        //         interval_group_end = intervals_total_size;

        //         in_intervals = true;

        //     }

        //     if (in_intervals && (*interval_groups)[i] != this_group) {

        //         interval_group_end = i;

        //         break;

        //     }

        // }

        // //std::cout << "interval_group_end: " << interval_group_end << "\n";

        // //if we're not in the intervals - just return without trying to match at all
        // if (!in_intervals) {

        //     return;

        // }

        //just for diagnostic purposes, add a counter of the number of loops.
        std::size_t nloops = 0;

        //to keep the number of loops down, and for locf-purposes, store the index of the last element before the start of the current interval
        std::size_t last_element_before = group_start;

        //this required some trial and error.....
        //use a typedef to the iterator to make it less wordy...
        typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::iterator ITERATOR_TYPE;
        //function pointer to some sort of aggregation function working on iterators
        ITERATOR_TYPE(*f)(ITERATOR_TYPE, ITERATOR_TYPE);

        typedef typename decltype(typename std::remove_pointer < decltype(times) > ::type())::value_type TIMES;
        TIMES locf_max = static_cast<TIMES>(locf);
        TIMES nocb_max = static_cast<TIMES>(nocb);
        //static_cast<TIMES>
        //std::cout << "locf_max: " << locf_max << "\n";

        switch(fun_aggregate){
            case 0: f = &std::min_element;
                    break;
            case 1: f = &std::max_element;
                    break;
            case 2: f = &first_element;
                    break;
            case 3: f = &last_element;
                    break;
            default: f = &std::min_element;
                    break;
        }

        auto aggregate_begin = values->begin();
        auto aggregate_end = values->begin();

        //loop through the intervals
        for (auto i = interval_group_start; i < interval_group_end; ++i) {
            bool matching_elements_found = false;
            bool result_set = false;
            //loop through values
            for (auto j = last_element_before; j < group_end; ++j) {
                ++nloops;
                //update the last element before the interval start, i.e. possible locf-value
                if ((*times)[j] < (*interval_begin)[i]) {
                    last_element_before = j;
                }
                else if (!matching_elements_found && (*times)[j] >= (*interval_begin)[i] && (*times)[j] <= (*interval_end)[i]) {
                    matching_elements_found = true;
                    aggregate_begin = values->begin() + j;
                    
                }
                //break the loop if the current time is larger than the current interval end or we're at the end of the group
                if (((*times)[j] > (*interval_end)[i]) || j == group_end - 1) {
                    //call aggregation function, if we've previously found elements
                    if (matching_elements_found) {
                        //the end iterator is either j, j+1, depending on the condition for break
                        aggregate_end = values->begin() + j;

                        //make sure we include the last value when it is also within aggregation
                        if ((j == (group_end - 1)) && ((*times)[j] <= (*interval_end)[i])){
                            ++aggregate_end;
                        }

                        //apply our aggregation function
                        auto target_position = (*f)(aggregate_begin, aggregate_end);

                        //save value and time from the correct index
                        if (target_position != values->end()) {
                            (*res)[i] = (*target_position);
                            auto t_index = std::abs(std::distance(values->begin(), target_position));
                            (*res_time)[i] = (*times)[t_index];
                            //we are in the interval, so we miss by 0
                            //(*res_time_difference)[i] = 0;
                            result_set = true;
                        }
                    }

                    //if we haven't set any result yet, possibly set result using locf
                    if (!result_set && ((*times)[last_element_before] < (*interval_begin)[i]) && ((*times)[last_element_before] + locf_max >= (*interval_begin)[i])) {
                        (*res)[i] = (*values)[last_element_before];
                        (*res_time)[i] = (*times)[last_element_before];
                        result_set = true;
                    }

                    //if that doesn't work, try nocb
                    if (!result_set && (*times)[j] > (*interval_end)[i] && (*times)[j] - nocb_max <= (*interval_end)[i]) {
                        (*res)[i] = (*values)[j];
                        (*res_time)[i] = (*times)[j];
                        result_set = true;
                    }

                    //std::cout << "interval = " << i << " last_element_before = " << last_element_before << "\n";
                    //finally, break the value loop
                    break;
                }
            }
        }
    }

};


//class for aggregating values over intervals, this time allowing overlapping intervals!
template < typename Td, typename Tg, typename Tv >
class interval_summarize {
private:
    //the data vectors
    Tg* groups;
    Td* times;
    Tv* values;

    //the intervals data vectors
    Tg* interval_groups;
    Td* interval_begin;
    Td* interval_end;

    //the results vector
    Tv* res;

    //our offsets in the data vectors
    std::size_t group_start;
    std::size_t group_end;

    //our offsets in the intervals vector
    std::size_t interval_group_start = 0;
    std::size_t interval_group_end = 0;

    //our aggregation function
    int fun_aggregate = 0;

    double locf;
    double nocb;


public:

    interval_summarize(Tg* groups, Td* times, Tv* values, Tg* interval_groups, Td* interval_begin, Td* interval_end, Tv* res, std::size_t group_start, std::size_t group_size, int fun_aggregate, double locf = 0.0, double nocb = 0.0) {
        //copy the data pointers
        this->groups = groups;
        this->times = times;
        this->values = values;

        this->interval_groups = interval_groups;
        this->interval_begin = interval_begin;
        this->interval_end = interval_end;

        this->res = res;

        this->group_start = group_start;
        this->group_end = group_start + group_size;
        this->fun_aggregate = fun_aggregate;
        this->locf = locf;
        this->nocb = nocb;

    }

    auto operator()() -> void {

        auto this_group = (*groups)[group_start];

        std::size_t intervals_total_size = interval_groups->size();
        //std::cout << "lei ho from group: " << this_group << "\n";
        //std::cout << "group_start: " << group_start << " group_end: " << group_end << "\n";

        //find our id in the intervals
        bool in_intervals = false;

        for (auto i = 0u; i < intervals_total_size; ++i) {

            if ((!in_intervals) && (*interval_groups)[i] == this_group) {

                interval_group_start = i;

                interval_group_end = intervals_total_size;

                in_intervals = true;

            }

            if (in_intervals && (*interval_groups)[i] != this_group) {

                interval_group_end = i;

                break;

            }

        }

        //std::cout << "interval_group_end: " << interval_group_end << "\n";

        //if we're not in the intervals - just return without trying to match at all
        if (!in_intervals) {

            return;

        }

        //just for diagnostic purposes, add a counter of the number of loops.
        std::size_t nloops = 0;

        //to keep the number of loops down, and for locf-purposes, store the index of the last element before the start of the current interval
        std::size_t last_element_before = group_start;

        //this required some trial and error.....
        //use a typedef to the iterator to make it less wordy...
        typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::iterator ITERATOR_TYPE;
        //typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::value_type VALUES;
        //function pointer to some sort of aggregation function working on iterators
        double (*f)(ITERATOR_TYPE, ITERATOR_TYPE);

        typedef typename decltype(typename std::remove_pointer < decltype(times) > ::type())::value_type TIMES;
        TIMES locf_max = static_cast<TIMES>(locf);
        TIMES nocb_max = static_cast<TIMES>(nocb);
        //static_cast<TIMES>
        //std::cout << "locf_max: " << locf_max << "\n";

        switch(fun_aggregate){
            case 0: f = &mean_na_zero;
                    break;
            case 1: f = &median_na_zero;
                    break;
            case 2: f = &first_na_zero;
                    break;
            case 3: f = &last_na_zero;
                    break;
            case 4: f = &var_na_zero;
                break;
            case 5: f = &sd_na_zero;
                break;
			case 6: f = &sum_na_zero;
                break;
            case 7: f = &min_na_zero;
                break;
            case 8: f = &max_na_zero;
                break;
            default: f = &last_na_zero;
                    break;
        }

        auto aggregate_begin = values->begin();
        auto aggregate_end = values->begin();

        //loop through the intervals
        for (auto i = interval_group_start; i < interval_group_end; ++i) {
            bool matching_elements_found = false;
            bool result_set = false;
            //loop through values
            for (auto j = last_element_before; j < group_end; ++j) {
                ++nloops;
                //update the last element before the interval start, i.e. possible locf-value
                if ((*times)[j] < (*interval_begin)[i]) {
                    last_element_before = j;
                }
                else if (!matching_elements_found && (*times)[j] >= (*interval_begin)[i] && (*times)[j] <= (*interval_end)[i]) {
                    matching_elements_found = true;
                    aggregate_begin = values->begin() + j;
                    
                }
                //break the loop if the current time is larger than the current interval end or we're at the end of the group
                if (((*times)[j] > (*interval_end)[i]) || j == group_end - 1) {
                    //call aggregation function, if we've previously found elements
                    if (matching_elements_found) {
                        //the end iterator is either j, j+1, depending on the condition for break
                        aggregate_end = values->begin() + j;

                        //make sure we include the last value when it is also within aggregation
                        if ((j == (group_end - 1)) && ((*times)[j] <= (*interval_end)[i])){
                            ++aggregate_end;
                        }
                        //just double check that the iterators don't point to the same element before trying to summarize
                        if (aggregate_begin != aggregate_end){
                            //set the result directly using the summarizing function
                            (*res)[i] = (*f)(aggregate_begin, aggregate_end);
                            result_set = true;
                        }
                    }

                    //if we haven't set any result yet, possibly set result using locf
                    if (!result_set && ((*times)[last_element_before] < (*interval_begin)[i]) && ((*times)[last_element_before] + locf_max >= (*interval_begin)[i])) {
                        (*res)[i] = (*values)[last_element_before];
                    }

                    //if that doesn't work, try nocb
                    if (!result_set && (*times)[j] > (*interval_end)[i] && (*times)[j] - nocb_max <= (*interval_end)[i]) {
                        (*res)[i] = (*values)[j];
                        result_set = true;
                    }
                    break;
                }
            }
        }
    }

};

//[[Rcpp::export]]
Rcpp::DataFrame cpp_interval_aggregate(Rcpp::NumericVector groups,
                                        Rcpp::NumericVector times,
                                        Rcpp::NumericVector values,
                                        Rcpp::NumericVector int_groups,
                                        Rcpp::NumericVector int_begin,
                                        Rcpp::NumericVector int_end,
                                        Rcpp::String aggregation_function = "min",
                                        double locf = 0.0,
                                        double nocb = 0.0,
                                        int nThreads = 0,
                                        int threadThreshold = 32,
                                        int threadMultiplier = 64,
                                        int minThreads = 1){
    //setup aggregation function
    int fun_aggregate = 0;
    if (aggregation_function == "min"){
        fun_aggregate = 0;
    } else if (aggregation_function == "max"){
        fun_aggregate = 1;
    } else if (aggregation_function == "first"){
        fun_aggregate = 2;
    } else if (aggregation_function == "last"){
        fun_aggregate = 3;
    } else {
        Rcpp::stop("'aggregation_function' must be one of 'min', 'max', 'first', or 'last'");
    }

    //some dynamic thread pool calculations
    if (nThreads <= 0){
        nThreads = std::max(1u,std::thread::hardware_concurrency());
        //Rcpp::Rcout << "nThreads: " << nThreads << "\n";
    }
    

    std::size_t n = groups.size();
	std::size_t n_int = int_groups.size();
    //the out vector, and the vector of futures
    Rcpp::NumericVector res(n_int, Rcpp::NumericVector::get_na());
	Rcpp::NumericVector res_time(n_int, Rcpp::NumericVector::get_na());
    //Rcpp::NumericVector res_time_difference(n_int, Rcpp::NumericVector::get_na());
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
                //only add the task if we're in the group
                if (in_intervals){
                    resvec.push_back(std::move(pool.add_task(

                    interval_aggregate_overlap < Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
                        &times,
                        &values,
                        &int_groups,
                        &int_begin,
                        &int_end, 
                        &res,
                        &res_time,
                        //&res_time_difference,
                        this_group_start,
                        this_group_size,
                        interval_group_start,
                        interval_group_end,
                        fun_aggregate,
                        locf,
                        nocb)

                    )));

                }

				//std::cout << "group found! \n";
				//std::cout << "this_group_start: " << this_group_start << " this group size: " << this_group_size << "\n";

				

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
	Rcpp::DataFrame out = Rcpp::DataFrame::create(Rcpp::Named("time") = res_time, Rcpp::Named("value") = res);

    return out;

}


//[[Rcpp::export]]
// Rcpp::DataFrame cpp_interval_overlap_aggregate(Rcpp::NumericVector groups,
//                                         Rcpp::NumericVector begin_times,
//                                         Rcpp::NumericVector end_times,
//                                         Rcpp::NumericVector values,
//                                         Rcpp::NumericVector int_groups,
//                                         Rcpp::NumericVector int_begin,
//                                         Rcpp::NumericVector int_end,
//                                         Rcpp::String aggregation_function = "min",
//                                         double locf = 0.0,
//                                         double nocb = 0.0,
//                                         int nThreads = 0,
//                                         int threadThreshold = 32,
//                                         int threadMultiplier = 64,
//                                         int minThreads = 1){
//     //setup aggregation function
//     int fun_aggregate = 0;
//     if (aggregation_function == "min"){
//         fun_aggregate = 0;
//     } else if (aggregation_function == "max"){
//         fun_aggregate = 1;
//     } else if (aggregation_function == "first"){
//         fun_aggregate = 2;
//     } else if (aggregation_function == "last"){
//         fun_aggregate = 3;
//     } else {
//         Rcpp::stop("'aggregation_function' must be one of 'min', 'max', 'first', or 'last'");
//     }

//     //some dynamic thread pool calculations
//     if (nThreads <= 0){
//         nThreads = std::max(1u,std::thread::hardware_concurrency());
//         //Rcpp::Rcout << "nThreads: " << nThreads << "\n";
//     }
    

//     std::size_t n = groups.size();
// 	std::size_t n_int = int_groups.size();
//     //the out vector, and the vector of futures
//     Rcpp::NumericVector res(n_int, Rcpp::NumericVector::get_na());
// 	Rcpp::NumericVector res_begin(n_int, Rcpp::NumericVector::get_na());
//     Rcpp::NumericVector res_end(n_int, Rcpp::NumericVector::get_na());
//     //Rcpp::NumericVector res_time_difference(n_int, Rcpp::NumericVector::get_na());
//     std::vector<std::future<void>> resvec;

//     //the location of the groups
//     std::size_t this_group = groups[0];
// 	std::size_t this_group_start = 0u;
// 	std::size_t this_group_size;

//     //the location of the groups within the intervals
//     std::size_t interval_group_start = 0u;
//     std::size_t interval_group_end = 0u;

//     //a flag to check if the current ID has been found within the intervals
//     bool in_intervals = false;

//     //the scope for the thread pool
//     {
//         threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
//         //Rcpp::Rcout << "Threadpool started...!\n";
//         for (auto i = 0u; i < n; ++i) {
// 			//set the size of the group
// 			if ((groups[i] != this_group) | (i == n - 1)) {
// 				if (i == n - 1) {
// 					this_group_size = i - this_group_start + 1;
// 				}
// 				else {
// 					this_group_size = i - this_group_start;
// 				}
                
//                 //new - check the intervals in the main function and pass locations to the tasks - more efficient for very large number of intervals
//                 in_intervals = false;
//                 for (auto j = interval_group_end;j<n_int;++j){
//                     if ((!in_intervals) && (int_groups[j] == this_group)){
//                         in_intervals = true;
//                         interval_group_start = j;
//                         interval_group_end = n_int;
//                     }

//                     if (in_intervals && (int_groups[j] != this_group )){
//                         interval_group_end = j;
//                         break;
//                     }
//                 }
//                 //only add the task if we're in the group
//                 if (in_intervals){
//                     resvec.push_back(std::move(pool.add_task(

//                     foverlap_aggregate < Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
//                         &begin_times,
//                         &end_times,
//                         &values,
//                         &int_groups,
//                         &int_begin,
//                         &int_end, 
//                         &res,
//                         &res_begin,
//                         &res_end,
//                         //&res_time_difference,
//                         this_group_start,
//                         this_group_size,
//                         interval_group_start,
//                         interval_group_end,
//                         fun_aggregate,
//                         locf,
//                         nocb)

//                     )));

//                 }

// 				//std::cout << "group found! \n";
// 				//std::cout << "this_group_start: " << this_group_start << " this group size: " << this_group_size << "\n";

				

// 				//resvec.push_back(std::move(pool.add_task(cof_dm<std::vector<double>, std::vector<double>::iterator>(data, coeffs, this_group_start, this_group_size, 4, 1, false))));

// 				//interval_match<std::vector<double>,std::vector<int>>(&groups, &times, &int_groups, &int_begin, &int_end, &res, this_group_start, this_group_size)
// 				this_group = groups[i];
// 				this_group_start = i;
// 			}
// 		}
//         //Rcpp::Rcout << "threadpool finished, nthreads = " << pool.get_size() << "\n";
//     }

//     for (auto &r : resvec){
// 		r.get();
// 	}
// 	Rcpp::DataFrame out = Rcpp::DataFrame::create(Rcpp::Named("begin") = res_begin, Rcpp::Named("end") = res_end, Rcpp::Named("value") = res);

//     return out;

// }

// //just mean
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


//[[Rcpp::export]]
Rcpp::NumericVector cpp_interval_summarize(Rcpp::NumericVector groups,
                                        Rcpp::NumericVector times,
                                        Rcpp::NumericVector values,
                                        Rcpp::NumericVector int_groups,
                                        Rcpp::NumericVector int_begin,
                                        Rcpp::NumericVector int_end,
                                        Rcpp::String aggregation_function = "mean",
                                        double locf = 0.0,
                                        double nocb = 0.0,
                                        int nThreads = 0,
                                        int threadThreshold = 32,
                                        int threadMultiplier = 64,
                                        int minThreads = 1){
    //setup aggregation function
    int fun_aggregate = 0;
    if (aggregation_function == "mean"){
        fun_aggregate = 0;
    } else if (aggregation_function == "median"){
        fun_aggregate = 1;
    } else if (aggregation_function == "first"){
        fun_aggregate = 2;
    } else if (aggregation_function == "last"){
        fun_aggregate = 3;
    } else if (aggregation_function == "var"){
        fun_aggregate = 4;
    } else if (aggregation_function == "sd"){
        fun_aggregate = 5;
	} else if (aggregation_function == "sum"){
        fun_aggregate = 6;
    } else if (aggregation_function == "min"){
        fun_aggregate = 7;
    } else if (aggregation_function == "max"){
        fun_aggregate = 8; 
    } else {
        Rcpp::stop("'aggregation_function' must be one of 'mean', 'median', 'first', 'last', 'var', 'sum' 'min', 'max' or 'sd'");
    }

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

				//std::cout << "group found! \n";
				//std::cout << "this_group_start: " << this_group_start << " this group size: " << this_group_size << "\n";

				resvec.push_back(std::move(pool.add_task(

                    interval_summarize < Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
                        &times,
                        &values,
                        &int_groups,
                        &int_begin,
                        &int_end, 
                        &res,
                        this_group_start,
                        this_group_size,
                        fun_aggregate,
                        locf,
                        nocb)

                )));

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


//class for matching on time or value within a certain time range
template < typename Td, typename Tg, typename Tv >
class interval_match_value {
private:
    //the data vectors
    Tg* groups;
    Td* times;
    Tv* values;

    //the intervals data vectors
    Tg* interval_groups;
    //Td* interval_begin;
    //Td* interval_end;
    Td* interval_times;
    Tv* interval_values;

    //the results vector
    Td* time_res;
    Tv* value_res;

    //our offsets in the data vectors
    std::size_t group_start;
    std::size_t group_end;

    //our offsets in the intervals vector
    std::size_t interval_group_start = 0;
    std::size_t interval_group_end = 0;

    //defining the time intervals
    double interval_time_before;
    double interval_time_after;

    //the penalty coefficients for matching
    double time_penalty_pow;
    double time_penalty_coeff;

    double value_penalty_pow;
    double value_penalty_coeff;



    //double locf;


public:

    interval_match_value(Tg* groups, 
    Td* times, 
    Tv* values, 
    Tg* interval_groups, 
    //Td* interval_begin, 
    //Td* interval_end, 
    Td* interval_times, 
    Tv* interval_values, 
    Td* time_res,
    Tv* value_res, 
    std::size_t group_start, 
    std::size_t group_size, 
    //int fun_aggregate, 
    double interval_time_before = 120,
    double interval_time_after = 120,
    double time_penalty_pow = 0.0,
    double time_penalty_coeff = 0.0,
    double value_penalty_pow = 1.0,
    double value_penalty_coeff = 1.0) {
        //copy the data pointers
        this->groups = groups;
        this->times = times;
        this->values = values;

        this->interval_groups = interval_groups;
        this->interval_times = interval_times;
        this->interval_values = interval_values;
        //this->interval_begin = interval_begin;
        //this->interval_end = interval_end;

        this->time_res = time_res;
        this->value_res = value_res;

        this->group_start = group_start;
        this->group_end = group_start + group_size;
        //this->fun_aggregate = fun_aggregate;

        this->interval_time_before = interval_time_before;
        this->interval_time_after = interval_time_after;
        
        this->time_penalty_pow = time_penalty_pow;
        this->time_penalty_coeff = time_penalty_coeff;

        this->value_penalty_pow = value_penalty_pow;
        this->value_penalty_coeff = value_penalty_coeff;


        //this->locf = locf;

    }

    auto operator()() -> void {

        auto this_group = (*groups)[group_start];

        std::size_t intervals_total_size = interval_groups->size();
        std::size_t values_total_size = values->size();

        //std::cout << "lei ho from group: " << this_group << "\n";
        //std::cout << "group_start: " << group_start << " group_end: " << group_end << "\n";

        //find our id in the intervals
        bool in_intervals = false;

        for (auto i = 0u; i < intervals_total_size; ++i) {

            if ((!in_intervals) && (*interval_groups)[i] == this_group) {

                interval_group_start = i;

                interval_group_end = intervals_total_size;

                in_intervals = true;

            }

            if (in_intervals && (*interval_groups)[i] != this_group) {

                interval_group_end = i;

                break;

            }

        }

        //std::cout << "interval_group_end: " << interval_group_end << "\n";

        //if we're not in the intervals - just return without trying to match at all
        if (!in_intervals) {

            return;

        }

        //just for diagnostic purposes, add a counter of the number of loops.
        std::size_t nloops = 0;

        //to keep the number of loops down, and for locf-purposes, store the index of the last element before the start of the current interval
        std::size_t last_element_before = group_start;

        std::size_t aggregate_begin = 0u;
        std::size_t aggregate_end = 0u;

        //use a typedef to fetch the type of the elements in the times-vector
        //typedef typename decltype(typename std::remove_pointer < decltype(times) > ::type())::value_type TIMES_TYPE;
        //or just hard-code to double to save the headache involved with types from Rcpp-vectors...
        double interval_begin = 0;
        double interval_end = 0;
        //do a static cast of the times before and after
        double interval_before = static_cast<double>(interval_time_before);
        double interval_after = static_cast<double>(interval_time_after);

        bool matching_element_found = false;
        //bool result_set = false;

        //loop through the intervals
        for (auto i = interval_group_start; i < interval_group_end; ++i) {
            matching_element_found = false;
            //result_set = false;

            interval_begin = (*interval_times)[i]-interval_before;
            interval_end = (*interval_times)[i]+interval_after;


            //loop through values
            for (auto j = last_element_before; j < group_end; ++j) {
                ++nloops;
                //update the last element before the interval start, i.e. possible locf-value
                if ((*times)[j] < interval_begin) {
                    last_element_before = j;
                }
                else if (!matching_element_found && (*times)[j] >= interval_begin && (*times)[j] <= interval_end) {
                    matching_element_found = true;
                    aggregate_begin = j;
                    
                }
                //break the loop if the current time is larger than the current interval end or we're at the end of the group
                if (((*times)[j] > interval_end) || j == group_end - 1) {
                    //call aggregation function, if we've previously found elements
                    if (matching_element_found) {
                        //the end point is either j, j+1, depending on the condition for break
                        aggregate_end = j;

                        //make sure we include the last value when it is also within aggregation
                        if ((j == (group_end - 1)) && ((*times)[j] <= interval_end)){
                            ++aggregate_end;
                        }
                        //just double check we actually have at least one element before trying to summarize
                        if (aggregate_begin != aggregate_end){
                            //find closest value
                            double current_difference;
                            double best_difference = std::numeric_limits<double>::max();
                            //double next_best_difference = best_difference;
                            std::size_t best_matching_element = std::numeric_limits<std::size_t>::max();
                            for (auto k = aggregate_begin;k<aggregate_end;++k){
                                //the sum of weighted absolute differences in time and value
                                current_difference = std::pow(std::fabs((*interval_values)[i] - (*values)[k]), value_penalty_pow) * value_penalty_coeff +
                                                     std::pow(std::fabs((*interval_times)[i] - (*times)[k]), time_penalty_pow) * time_penalty_coeff;

                                if (current_difference < best_difference){
                                    best_difference = current_difference;
                                    best_matching_element = k;
                                }
                            }

                            //if we've found a match, set the results
                            if (best_matching_element < values_total_size){
                                (*time_res)[i] = (*times)[best_matching_element];
                                (*value_res)[i] = (*values)[best_matching_element];
                                //result_set = true;
                            }
                            
                            

                            //set the result directly using the summarizing function
                            //(*res)[i] = (*f)(aggregate_begin, aggregate_end);
                            
                        }
                    }

                    //if we haven't set any result yet, possibly set result using locf
                    //if (!result_set && ((*times)[last_element_before] < (*interval_begin)[i]) && ((*times)[last_element_before] + locf_max >= (*interval_begin)[i])) {
                    //    (*res)[i] = (*values)[last_element_before];
                    //}

                    //std::cout << "interval = " << i << " last_element_before = " << last_element_before << "\n";
                    //finally, break the value loop
                    break;
                }
            }
        }
    }

};



//class for matching the closest value on time
//kind-of like a possibly bidirectional rolling join from R's data.table
//but automatically chosing the closest value within the c++ code
template < typename Td, typename Tg, typename Tv >
class interval_match_time {
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

public:

    interval_match_time(Tg* groups, 
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
    double interval_time_after = 120
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
        this->interval_time_before = interval_time_before;
        this->interval_time_after = interval_time_after;
        this->interval_group_start = interval_group_start;
        this->interval_group_end = interval_group_end;


    }

    auto operator()() -> void {

        //auto this_group = (*groups)[group_start];

        //std::size_t intervals_total_size = interval_groups->size();
        //std::size_t values_total_size = values->size();

        //find our id in the intervals
        //bool in_intervals = false;

        // for (auto i = 0u; i < intervals_total_size; ++i) {
        //     if ((!in_intervals) && (*interval_groups)[i] == this_group) {
        //         interval_group_start = i;
        //         interval_group_end = intervals_total_size;
        //         in_intervals = true;
        //     }

        //     if (in_intervals && (*interval_groups)[i] != this_group) {
        //         interval_group_end = i;
        //         break;
        //     }
        // }

        //if we're not in the intervals - just return without trying to match at all
        //if (!in_intervals) {
        //    return;
        //}

        //just for diagnostic purposes, add a counter of the number of loops.
        //std::size_t nloops = 0;

        //to keep the number of loops down, and for locf-purposes, store the index of the last element before the start of the current interval
        //std::size_t last_element_before = group_start;

        //std::size_t best_matching_element = 0;
        std::size_t first_time_after = group_start;
        std::size_t best_match = group_start;
        bool matching_found = false;
        double best_difference = 0;

        //loop through the intervals
        for (auto i = interval_group_start;i<interval_group_end;++i){
            matching_found = false;
            //move the index in times as far as possible,
            //i.e. to the first time point greater or equal to the requested time point
            while (((*times)[first_time_after]<(*interval_times)[i]) && ((first_time_after<group_end))){
                ++first_time_after;
            }

            //now, fill the result with either the value at best_location, or the one before it
            //depending on which one is closest
            
            //step one, assume the one that might be equal to be the closest one
            //note, these &&-operators will be required as the pointer dereferencing could be invalid 
            //if the left part of the && is false (thus, & is impossible)!
            if ((first_time_after != group_end) && (((*times)[first_time_after] - (*interval_times)[i]) <= interval_time_after)){
                matching_found = true;
                best_match = first_time_after;
                best_difference = (*times)[first_time_after] - (*interval_times)[i];
            }
            //step two, check if the one before it might be even closer
            if ((first_time_after > group_start) && (((*interval_times)[i] - (*times)[first_time_after-1]) <= interval_time_before)){
                if (matching_found){
                    if ((*interval_times)[i]- (*times)[first_time_after-1] < best_difference){
                        best_match = first_time_after-1;
                    }
                } else {
                    matching_found = true;
                    best_match = first_time_after-1;
                }
            }

            if (matching_found){
                (*res)[i] = (*values)[best_match];
            }
        }
    }

};


//[[Rcpp::export]]
Rcpp::DataFrame cpp_match_value(Rcpp::NumericVector groups,
                                    Rcpp::NumericVector times,
                                    Rcpp::NumericVector values,
                                    Rcpp::NumericVector int_groups,
                                    Rcpp::NumericVector int_times,
                                    Rcpp::NumericVector int_values,
                                    double interval_time_before = 120,
                                    double interval_time_after = 120,
                                    double time_coeff = 0.0,
                                    double value_coeff = 1.0,
                                    double time_pow = 1.0,
                                    double value_pow = 1.0,
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
    //the out vectors, and the vector of futures
    Rcpp::NumericVector res_time(n_int, Rcpp::NumericVector::get_na());
    Rcpp::NumericVector res_value(n_int, Rcpp::NumericVector::get_na());
    std::vector<std::future<void>> resvec;

    //the location of the groups
    std::size_t this_group = groups[0];
	std::size_t this_group_start = 0u;
	std::size_t this_group_size;

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

				//std::cout << "group found! \n";
				//std::cout << "this_group_start: " << this_group_start << " this group size: " << this_group_size << "\n";

				resvec.push_back(std::move(pool.add_task(

                    interval_match_value < Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
                        &times,
                        &values,
                        &int_groups,
                        &int_times,
                        &int_values, 
                        &res_time,
                        &res_value,
                        this_group_start,
                        this_group_size,
                        interval_time_before,
                        interval_time_after,
                        time_pow,
                        time_coeff,
                        value_pow,
                        value_coeff


                        )

                )));

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

    Rcpp::DataFrame res_frame = Rcpp::DataFrame::create(Rcpp::Named("time") =  res_time, Rcpp::Named("value") = res_value);

    return res_frame;
}

//a simpler, faster version of cpp_match_value that only matches on time and returns only the matched value
//this also does the grouping of intervals in the main calling function, which is probably more efficient
//[[Rcpp::export]]
Rcpp::NumericVector cpp_match_time(Rcpp::NumericVector groups,
                                   Rcpp::NumericVector times,
                                   Rcpp::NumericVector values,
                                   Rcpp::NumericVector int_groups,
                                   Rcpp::NumericVector int_times,
                                   double interval_time_before = 120,
                                   double interval_time_after = 120,
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
    //A loop through the first grouping vectorr
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
              interval_match_time < Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>(&groups,
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
                                                                                                   interval_time_after
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
  
  //Rcpp::DataFrame res_frame = Rcpp::DataFrame::create(Rcpp::Named("time") =  res_time, Rcpp::Named("value") = res_value);
  
  return res;
}