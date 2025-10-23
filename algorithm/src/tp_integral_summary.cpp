// 2024-11-12
// A file for limited integral calculation based on intermittent values
// No guarantees to be the proven optimal mathematical solution
// Clock cycles go a long way when the intellect is lacking

#include <Rcpp.h>
//[[Rcpp::plugins("cpp14")]]

#include <algorithm> //for nth_element
#include <vector>
#include <numeric> //for accumulate
#include <cmath>
#include <iterator>

// for the thread pool
#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <future>
#include <queue>

// 2021-04-15
// A dynamically resizing version of the old pool
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

    std::size_t get_size()
    {
        return threads.size();
    }
    // a debugging method for adding tasks that don't return anything...
    void add_task_lao(Task task)
    {
        {
            std::unique_lock<std::mutex> lock(event_mutex);
            tasks.emplace(std::move(task));
        }
        event_var.notify_one();
    }
    // returns an std::future of the return type of the functor passed onto task...
    template <class T>
    auto add_task(T task) -> std::future<decltype(task())>
    {
        // a shared pointer to a packaged_task containing the ()-operator of the task-functor...
        auto wrapper = std::make_shared<std::packaged_task<decltype(task())()>>(std::move(task));
        std::size_t waiting_tasks;
        {
            std::unique_lock<std::mutex> lock(event_mutex);
            // get the amount of tasks currently on the pool, to see if we should start new threads
            waiting_tasks = tasks.size();
            // add the functor to the tasks-queue
            tasks.emplace([=]
                          { (*wrapper)(); });

            // This appears to have to be within the same block as the mutex to avoid race conditions that would sometimes
            // crash the pool when the tasks are very small and there are many cores
            event_var.notify_one();
        }
        // the 'dynamic size' is very simplistic and basically starts threads until the size is maximum,
        // if waiting_tasks > threadThreshhold + threadMultiplier * current amount of active threads
        // i.e. if each thread has more than threadMultiplier tasks yet to be completed
        if ((waiting_tasks > (threadThreshhold + threadMultiplier * threads.size())) & (threads.size() < maxThreads))
        {
            add_thread();
        }
        // get the future from our shared pointer
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

    void add_thread()
    {
        // Rcpp::Rcout << "adding a thread...\n";
        // this lambda contains the thread logic -
        threads.emplace_back([=]
                             {
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

			} });
    }

    void stop() noexcept
    {
        {
            std::unique_lock<std::mutex> lock(event_mutex);
            stop_threads = true;
        }
        event_var.notify_all();

        for (auto &t : threads)
        {
            t.join();
        }
    }
};

// an absolute solution for straight-up carry forward of last y at any x
// it assumes that index1 is the first index such that (x)[i] < xmin
// and that index2 si the first index such that (x)[i] >= xmax
template <typename Tv>
double cf_time_sum(const Tv &x, const Tv &y, std::size_t index1, std::size_t index2, double xmin, double xmax)
{
    double sum = 0.0;
    double last_y = y[index1];
    double begin_x = xmin;

    // Iterate up to the second-to-last index (if there are multiple elements in range)
    for (auto i = index1 + 1; i < index2 - 1; ++i)
    {
        double end_x = x[i];
        sum += last_y * (end_x - begin_x);
        begin_x = end_x;
        last_y = y[i];
    }

    // Handle the last interval separately
    sum += last_y * (xmax - begin_x);

    return sum;
}

template <typename Tv>
class concentration_at_times{
    private:
    Tv *times;
    Tv *doserates;
    Tv *doserate_times;
    Tv *boluses;
    Tv *bolus_times;
    Tv *res;

    double vd;
    double halflife;
    double max_halflives;
    double max_interval_length;
    
    std::size_t times_begin;
    std::size_t times_end;
    std::size_t doserates_begin;
    std::size_t doserates_end;
    std::size_t boluses_begin;
    std::size_t boluses_end;

    public:
    concentration_at_times(Tv* times,
    Tv* doserates,
    Tv* doserate_times,
    Tv* boluses,
    Tv* bolus_times,
    Tv* res, 
    std::size_t times_begin,
    std::size_t times_end,
    std::size_t doserates_begin,
    std::size_t doserates_end,
    std::size_t boluses_begin,
    std::size_t boluses_end,
    double vd,
    double halflife,
    double max_halflives,
    double max_interval_length){
        //Rcpp::Rcout << "in constructor \n";
        this->times = times;
        this->doserates = doserates;
        this->doserate_times = doserate_times;
        this->boluses = boluses;
        this->bolus_times = bolus_times;
        this->res = res;
        this->times_begin = times_begin;
        this->times_end = times_end;
        this->doserates_begin = doserates_begin;
        this->doserates_end = doserates_end;
        this->boluses_begin = boluses_begin;
        this->boluses_end = boluses_end;

        this->vd = vd;
        this->halflife = halflife;
        this->max_halflives = max_halflives;
        this->max_interval_length = max_interval_length;

        //Rcpp::Rcout << "member variables initialized \n";


    }

    auto operator()()->void{
        double lambda = std::log(2) / halflife;
        double maxtime = max_halflives * halflife;
        std::size_t doserates_n = doserates_end - doserates_begin;
        
        //Rcpp::Rcout << "in operator ()\n";

        //store information about current doserate
        double current_doserate = 0;
        double doserate_start = 0;
        double doserate_end = 0;
        double doserate_duration = 0;
        //store information about smaller steps
        double n_steps = 0;
        double step_size = 0;
        double step_start = 0;
        double step_end = 0;

        //the current dose
        double dose = 0;
        //used instead of (*times)[i] within the loop for legibility
        double t = 0;
        //used instad of (*res)[i] within the loop for legibility
        double concentration = 0;
        for (auto i = times_begin;i<times_end;++i){
            //reset the concentration
            concentration = 0;
            //use t instead of (*times)[i] as it's shorter
            t = (*times)[i];
            if (doserates_n > 1){
                for (auto j = doserates_begin;j < doserates_end-1;++j){
                    current_doserate = (*doserates)[j];
                    doserate_start = (*doserate_times)[j];
                    doserate_end = (*doserate_times)[j+1u];
                    //Rcpp::Rcout << "j: " << j << "\n";
                    //Rcpp::Rcout << "doserate_start : " << doserate_start << " doserate_end: " << doserate_end << "\n";

                    //this doserate could potentially contribute to the concentration at t
                    if ((doserate_start < t) && (t<doserate_end + maxtime)){
                        doserate_end = std::min(t, doserate_end);
                        //Rcpp::Rcout << "hit! doserate_start : " << doserate_start << " doserate_end: " << doserate_end << "\n";
                        doserate_duration = doserate_end - doserate_start;

                        n_steps = std::ceil(doserate_duration / max_interval_length);
                        //Rcpp::Rcout << "in " << n_steps << " steps\n";

                        //if it's shorter than the interval length, count a single bolus dose in the middle of the interval
                        if (doserate_duration <= max_interval_length){
                            dose = current_doserate * doserate_duration;
                            concentration += (dose / vd) * std::exp(-lambda * (t - 0.5 * (doserate_end + doserate_start)));


                        } else {
                            //discretize into equal timesteps shorter than max_interval_length
                            //n_steps = std::ceil(doserate_duration / max_interval_length);
                            step_size = doserate_duration / n_steps;
                            //Rcpp::Rcout << "split into " << n_steps << " steps\n";
                            

                            for (auto k = 0u;k<n_steps;++k){
                                step_start = doserate_start + k * step_size;
                                step_end = doserate_start + (k + 1) * step_size;

                                //this step could potentially contribute to the concentration at t
                                if ((step_start < t) && (t<step_end + maxtime)){
                                    step_end = std::min(t, step_end);
                                    dose = current_doserate * (step_end - step_start);
                                    concentration += (dose / vd) * std::exp(-lambda * (t - 0.5 * (step_end + step_start)));
                                }
                            }
                        }




                    }
                    
                    //if the end of this interval is more than maxtime from the current t, it doesn't need to be checked again
                    if ((t > doserate_end + maxtime)){
                        doserates_begin = j;
                    }

                    //break the loop once the starttimes go beyond t
                    if (doserate_start > t) break;

                }
            }

            //then check contributions from bolus doses
            for (auto j = boluses_begin;j<boluses_end;++j){
                if (((*bolus_times)[j] < t) &&  (t < (*bolus_times)[j]+maxtime)){
                    concentration += ((*boluses)[j] / vd) * std::exp(-lambda * (t - (*bolus_times)[j]));
                }


                //update the loop index if the t is too far off in the future = we don't need to check the next t for this bolus
                if ((*bolus_times)[j]+maxtime<t){
                    boluses_begin = j;
                }

                //break the inner loop if the time is beyond t
                if ((*bolus_times)[j] > t){
                    break;
                }
            }

            (*res)[i] = concentration;
        }

        //loop through every timepoint
        if (false){
            for (auto i = times_begin;i<times_end;++i){
            //reset the concentration
            concentration = 0;
            //use t instead of (*times)[i] as it's shorter
            t = (*times)[i];

            //first check contributions from continous infusions
            if (doserates_n > 1){
                for (auto j = doserates_begin;doserates_end-1;++j){
                    current_doserate = (*doserates)[j];
                    doserate_start = (*doserate_times)[j];
                    doserate_end = (*doserate_times)[j+1u];

                    //this doserate could potentially contribute to the concentration at t
                    if ((doserate_start < t) && (t<doserate_end + maxtime)){
                        doserate_end = std::min(t, doserate_end);
                        doserate_duration = doserate_end - doserate_start;

                        //if it's shorter than the interval length, count a single bolus dose in the middle of the interval
                        if (doserate_duration <= max_interval_length){
                            dose = current_doserate * doserate_duration;
                            concentration += (dose / vd) * std::exp(-lambda * (t - 0.5 * (doserate_end + doserate_start)));


                        } else {
                            //discretize into equal timesteps shorter than max_interval_length
                            std::size_t n_steps = std::ceil(doserate_duration / max_interval_length);
                            step_size = doserate_duration / n_steps;

                            for (auto k = 0u;k<n_steps;++k){
                                step_start = doserate_start + k * step_size;
                                step_end = doserate_start + (k + 1) * step_size;

                                //this step could potentially contribute to the concentration at t
                                if ((step_start < t) && (t<step_end + maxtime)){
                                    step_end = std::min(t, step_end);
                                    dose = current_doserate * (step_end - step_start);
                                    concentration += (dose / vd) * std::exp(-lambda * (t - 0.5 * (step_end + step_start)));
                                }
                            }
                        }




                    }
                    
                    //if the end of this interval is more than maxtime from the current t, it doesn't need to be checked again
                    if ((t > doserate_end + maxtime)){
                        doserates_begin = j;
                    }

                    //break the loop once the starttimes go beyond t
                    if (doserate_start > t) break;

                }
            }

            
            //then check contributions from bolus doses
            for (auto j = boluses_begin;j<boluses_end;++j){
                if (((*bolus_times)[j] < t) &&  (t < (*bolus_times)[j]+maxtime)){
                    concentration += ((*boluses)[j] / vd) * std::exp(-lambda * (t - (*bolus_times)[j]));
                }


                //update the loop index if the t is too far off in the future = we don't need to check the next t for this bolus
                if ((*bolus_times)[j]+maxtime<t){
                    boluses_begin = j;
                }

                //break the inner loop if the time is beyond t
                if ((*bolus_times)[j] > t){
                    break;
                }
            }

            (*res)[i] = concentration;
        }

        }

        
    }

};


struct tv_pair{
    double time;
    double value;

    tv_pair(double time, double value) : time(time), value(value) {}
};

template <typename Tv>
class concentration_at_times_2{
    private:
    Tv *times;
    Tv *doserates;
    Tv *doserate_times;
    Tv *boluses;
    Tv *bolus_times;
    Tv *res;

    double vd;
    double halflife;
    double max_halflives;
    double max_interval_length;
    
    std::size_t times_begin;
    std::size_t times_end;
    std::size_t doserates_begin;
    std::size_t doserates_end;
    std::size_t boluses_begin;
    std::size_t boluses_end;

    public:
    concentration_at_times_2(Tv* times,
    Tv* doserates,
    Tv* doserate_times,
    Tv* boluses,
    Tv* bolus_times,
    Tv* res, 
    std::size_t times_begin,
    std::size_t times_end,
    std::size_t doserates_begin,
    std::size_t doserates_end,
    std::size_t boluses_begin,
    std::size_t boluses_end,
    double vd,
    double halflife,
    double max_halflives,
    double max_interval_length){
        //Rcpp::Rcout << "in constructor \n";
        this->times = times;
        this->doserates = doserates;
        this->doserate_times = doserate_times;
        this->boluses = boluses;
        this->bolus_times = bolus_times;
        this->res = res;
        this->times_begin = times_begin;
        this->times_end = times_end;
        this->doserates_begin = doserates_begin;
        this->doserates_end = doserates_end;
        this->boluses_begin = boluses_begin;
        this->boluses_end = boluses_end;

        this->vd = vd;
        this->halflife = halflife;
        this->max_halflives = max_halflives;
        this->max_interval_length = max_interval_length;

        //Rcpp::Rcout << "member variables initialized \n";


    }

    auto operator()()->void{
        double lambda = std::log(2) / halflife;
        double maxtime = max_halflives * halflife;
        std::size_t doserates_n = doserates_end - doserates_begin;

        //Rcpp::Rcout << "vd: " << vd << " halflife: " << halflife << " max halflives: " << max_halflives << " maxtime: " << maxtime << " lambda: " << lambda << "\n";

        
        //store information about current doserate
        double current_doserate = 0;
        //double current_doserate_time = 0;
        double next_doserate_time = 0;

        //these will be used for the smaller intervals
        double doserate_start = 0;
        double doserate_end = 0;
        double doserate_duration = 0;

        //we're going to split the doserate-times into smaller intervals for the concentration calculations
        std::vector<double> vec_doserate_times;

        //but for the doserates, only track changes...
        std::vector<tv_pair> doserate_updates;

        //store the loop coordinates for of our new doserate vector
        std::size_t vec_doserates_times_begin = 0u;
        std::size_t vec_doserates_times_end = 0u;

        //and for our new doserate-changes-vector
        //std::size_t doserate_updates_begin = 0u;
        //std::size_t doserate_updates_end = 0u;

        //if we have doserates - let's split them into smaller pieces
        if (doserates_n > 1){
            //store information about smaller steps
            std::size_t n_steps = 0;
            double step_size = 0;

            current_doserate = (*doserates)[doserates_begin];
            doserate_updates.emplace_back((*doserate_times)[doserates_begin], current_doserate);
            
            //careful when dereferencing this iterator, Eugene!
            //current_doserate = *(--vec_doserate_updates.end());;

            //the idea is to split the doserates into a std::vector of smaller steps
            for (auto i = doserates_begin;i<doserates_end-1;++i){
                //see if the doserate has been changed so we need to store an update
                if ((*doserates)[i] != current_doserate){
                    current_doserate = (*doserates)[i];       
                    doserate_updates.emplace_back((*doserate_times)[i], current_doserate);
                }

                doserate_start = (*doserate_times)[i];
                doserate_duration = (*doserate_times)[i+1u] - doserate_start;
                    
                n_steps = std::ceil(doserate_duration / max_interval_length);
                step_size = doserate_duration / n_steps;

                //push back intervals of equal size
                for (auto j = 0u;j<n_steps;++j){
                    vec_doserate_times.push_back(doserate_start + step_size * j);
                }
            }
            //add the last one as well
            vec_doserate_times.push_back((*doserate_times)[doserates_end-1]);

            vec_doserates_times_end = vec_doserate_times.size()-1u;

            //store the last update time too if we'll ever need it
            doserate_updates.emplace_back((*doserate_times)[doserates_end-1], (*doserates)[doserates_end-1]);
            //doserate_updates_end = doserate_updates.size() - 1u;
        }
        
        
        // Rcpp::Rcout << "vec_doserate_times: \n";
        // for (auto i = 0u;i<vec_doserate_times.size();++i){
        //    Rcpp::Rcout << vec_doserate_times[i] << "\n";
        // }

        // Rcpp::Rcout << "doserate_updates: \n";
        // for (auto i = 0u;i<doserate_updates.size();++i){
        //    Rcpp::Rcout << doserate_updates[i].time << ", " << doserate_updates[i].value << "\n";
        // }
        


        //the current dose
        double dose = 0;
        //used instead of (*times)[i] within the loop for legibility
        double t = 0;
        //used instad of (*res)[i] within the loop for legibility
        double concentration = 0;
        //track the index of our current doserate within the loop
        std::size_t current_doserate_update_idx = 0;
        std::size_t doserate_idx_begin = 0;

        for (auto i = times_begin;i<times_end;++i){
            //reset the concentration
            concentration = 0;
            //use t instead of (*times)[i] as it's shorter
            t = (*times)[i];


            //Rcpp::Rcout << "doserates_begin: " << doserates_begin << " doserates_end: " << doserates_end << "\n";
            //Rcpp::Rcout << "(*doserates).size(): " << (*doserates).size() << "\n";

            if (doserates_n > 1){
                
                current_doserate_update_idx = doserate_idx_begin;
                current_doserate = doserate_updates[current_doserate_update_idx].value;
                next_doserate_time = doserate_updates[current_doserate_update_idx+1].time;

                //Rcpp::Rcout << "i: " << i << "\n";

                for (auto j = vec_doserates_times_begin;j < vec_doserates_times_end;++j){
                    //current_doserate = (*doserates)[j];
                    doserate_start = vec_doserate_times[j];
                    doserate_end = vec_doserate_times[j+1u];

                    //update the doserate if needed
                    if (doserate_start >= next_doserate_time){
                        current_doserate_update_idx+=1u;
                        current_doserate = doserate_updates[current_doserate_update_idx].value;
                        next_doserate_time = doserate_updates[current_doserate_update_idx+1].time;

                    }  
                    

                    //this doserate could potentially contribute to the concentration at t
                    if ((doserate_start < t) && (t<doserate_end + maxtime)){
                        doserate_end = std::min(t, doserate_end);
                        //Rcpp::Rcout << "hit! doserate_start : " << doserate_start << " doserate_end: " << doserate_end << "\n";
                        doserate_duration = doserate_end - doserate_start;
                        dose = current_doserate * doserate_duration;

                        
                        //Rcpp::Rcout << "doserate_start: " << doserate_start << " current doserate: " << current_doserate << " dose: " << dose << "\n";
                        concentration += (dose / vd) * std::exp(-lambda * (t - 0.5 * (doserate_end + doserate_start)));

                    }
                    
                    //if the end of this interval is more than maxtime from the current t, it doesn't need to be checked again
                    if ((t > doserate_end + maxtime)){
                        vec_doserates_times_begin = j;
                    }

                    //dito for the shorter vector of dose changes
                    if ((t > doserate_updates[doserate_idx_begin+1u].time + maxtime)){
                        doserate_idx_begin+=1;
                    }

                    //break the loop once the starttimes go beyond t
                    if (doserate_start > t) break;

                }
            }

            //then check contributions from bolus doses
            for (auto j = boluses_begin;j<boluses_end;++j){
                if (((*bolus_times)[j] < t) &&  (t < (*bolus_times)[j]+maxtime)){
                    concentration += ((*boluses)[j] / vd) * std::exp(-lambda * (t - (*bolus_times)[j]));
                }


                //update the loop index if the t is too far off in the future = we don't need to check the next t for this bolus
                if ((*bolus_times)[j]+maxtime<t){
                    boluses_begin = j;
                }

                //break the inner loop if the time is beyond t
                if ((*bolus_times)[j] > t){
                    break;
                }
            }

            (*res)[i] = concentration;
        }
        
    }

};


template <typename Tv>
class concentration_at_times_3{
    private:
    Tv *times;
    Tv *doserates;
    Tv *doserate_times;
    Tv *boluses;
    Tv *bolus_times;
    Tv *res;

    double vd;
    double halflife;
    double max_halflives;
    double max_interval_length;
    
    std::size_t times_begin;
    std::size_t times_end;
    std::size_t doserates_begin;
    std::size_t doserates_end;
    std::size_t boluses_begin;
    std::size_t boluses_end;

    public:
    concentration_at_times_3(Tv* times,
    Tv* doserates,
    Tv* doserate_times,
    Tv* boluses,
    Tv* bolus_times,
    Tv* res, 
    std::size_t times_begin,
    std::size_t times_end,
    std::size_t doserates_begin,
    std::size_t doserates_end,
    std::size_t boluses_begin,
    std::size_t boluses_end,
    double vd,
    double halflife,
    double max_halflives,
    double max_interval_length){
        //Rcpp::Rcout << "in constructor \n";
        this->times = times;
        this->doserates = doserates;
        this->doserate_times = doserate_times;
        this->boluses = boluses;
        this->bolus_times = bolus_times;
        this->res = res;
        this->times_begin = times_begin;
        this->times_end = times_end;
        this->doserates_begin = doserates_begin;
        this->doserates_end = doserates_end;
        this->boluses_begin = boluses_begin;
        this->boluses_end = boluses_end;

        this->vd = vd;
        this->halflife = halflife;
        this->max_halflives = max_halflives;
        this->max_interval_length = max_interval_length;

        //Rcpp::Rcout << "member variables initialized \n";


    }

    auto operator()()->void{
        double lambda = std::log(2) / halflife;
        double maxtime = max_halflives * halflife;
        std::size_t doserates_n = doserates_end - doserates_begin;

        //Rcpp::Rcout << "vd: " << vd << " halflife: " << halflife << " max halflives: " << max_halflives << " maxtime: " << maxtime << " lambda: " << lambda << "\n";

        
        //store information about current doserate
        double current_doserate = 0;
        //double current_doserate_time = 0;
        double next_doserate_time = 0;

        //these will be used for the smaller intervals
        double doserate_start = 0;
        double doserate_end = 0;
        double doserate_duration = 0;

        //we're going to split the doserate-times into smaller intervals for the concentration calculations
        //std::vector<double> vec_doserate_times;

        //only track changes of doserate...
        std::vector<tv_pair> doserate_updates;
        //std::size_t doserate_updates_begin = 0u;
        //std::size_t doserate_updates_end;

        //store the loop coordinates for of our new doserate vector
        std::size_t vec_doserates_times_begin = 0u;
        std::size_t vec_doserates_times_end = 0u;

        //and for our new doserate-changes-vector
        std::size_t doserate_updates_begin = 0u;
        std::size_t doserate_updates_end = 0u;

        //if we have doserates - let's split them into smaller pieces
        if (doserates_n > 1){
            //store information about smaller steps
            //std::size_t n_steps = 0;
            //double step_size = 0;

            current_doserate = (*doserates)[doserates_begin];
            doserate_updates.emplace_back((*doserate_times)[doserates_begin], (*doserates)[doserates_begin]);
            
            //careful when dereferencing this iterator, Eugene!
            //current_doserate = *(--vec_doserate_updates.end());;

            //the idea is to split the doserates into a std::vector of smaller steps
            for (auto i = doserates_begin;i<doserates_end;++i){
                //see if the doserate has been changed so we need to store an update
                if ((*doserates)[i] != current_doserate){
                    current_doserate = (*doserates)[i];
                    doserate_updates.emplace_back((*doserate_times)[i], (*doserates)[i]);
                }

                //doserate_start = (*doserate_times)[i];
                //doserate_duration = (*doserate_times)[i+1u] - doserate_start;
                    
                //n_steps = std::ceil(doserate_duration / max_interval_length);
                //step_size = doserate_duration / n_steps;

                //push back intervals of equal size
                //for (auto j = 0u;j<n_steps;++j){
                //    vec_doserate_times.push_back(doserate_start + step_size * j);
                //}
            }
            //add the last one as well
            //vec_doserate_times.push_back((*doserate_times)[doserates_end-1]);

            //vec_doserates_times_end = vec_doserate_times.size()-1u;

            //store the last update time too if we'll ever need it
            //doserate_updates.emplace_back((*doserate_times)[doserates_end-1], (*doserates)[doserates_end-1]);
            //doserate_updates_end = doserate_updates.size() - 1u;
        }
        

        
        
        // Rcpp::Rcout << "vec_doserate_times: \n";
        // for (auto i = 0u;i<vec_doserate_times.size();++i){
        //    Rcpp::Rcout << vec_doserate_times[i] << "\n";
        // }

        // Rcpp::Rcout << "doserate_updates: \n";
        // for (auto i = 0u;i<doserate_updates.size();++i){
        //    Rcpp::Rcout << doserate_updates[i].time << ", " << doserate_updates[i].value << "\n";
        // }
        


        //the current dose
        double dose = 0;
        //used instead of (*times)[i] within the loop for legibility
        double t = 0;
        //used instad of (*res)[i] within the loop for legibility
        double concentration = 0;
        //track the index of our current doserate within the loop
        //std::size_t current_doserate_update_idx = 0;
        //std::size_t doserate_idx_begin = 0;

        std::size_t steps_from_dosechange = 0u;
        std::size_t step_n = 0u;
        double step_start = 0u;
        double step_end = 0u;
        

        for (auto i = times_begin;i<times_end;++i){
            //reset the concentration
            concentration = 0;
            //use t instead of (*times)[i] as it's shorter
            t = (*times)[i];


            //Rcpp::Rcout << "doserates_begin: " << doserates_begin << " doserates_end: " << doserates_end << "\n";
            //Rcpp::Rcout << "(*doserates).size(): " << (*doserates).size() << "\n";

            if (doserates_n > 1){
                for (auto j = doserate_updates_begin;j<doserate_updates_end-1;++j){
                    //break the loop if we've come too far
                    if (doserate_updates[j].time > t) break;

                    //this is part of this doserate-interval that could possibly affect the concentration at t
                    doserate_end_time = std::min(doserate_updates[j+1].time, t);
                    doserate_start_time = std::max(doserate_updates[j].time, t - maxtime - max_interval_length);

                    step_start = doserate_start_time;
                    step_end = doserate_start_time;
                    //the subdivision of the intervals
                    while (step_end < doserate_end_time){
                        step_end = std::min(doserate_start_time + max_interval_length, t);
                        doserate_duration = step_end - step_start;
                        dose = doserate_updates[j].value * doserate_duration;
                        concentration += (dose/vd) * std::exp(-lambda * (t- 0.5 * (step_end + step_start)));
                        step_start = step_end;
                    }

                    //if this index is too far back in time for the t =  (*times)[i], it won't need to be checked for (*times)[i+1], 
                    if (doserate_updates[j+1].time + maxtime < t) {
                        doserate_updates_begin = j;
                    }
                }

            }


            // if (doserates_n > 1){
                
            //     current_doserate_update_idx = doserate_idx_begin;
            //     current_doserate = doserate_updates[current_doserate_update_idx].value;
            //     next_doserate_time = doserate_updates[current_doserate_update_idx+1].time;

            //     for (auto j = 0u;j<)



            //     //Rcpp::Rcout << "i: " << i << "\n";

            //     for (auto j = vec_doserates_times_begin;j < vec_doserates_times_end;++j){
            //         //current_doserate = (*doserates)[j];
            //         doserate_start = vec_doserate_times[j];
            //         doserate_end = vec_doserate_times[j+1u];

            //         //update the doserate if needed
            //         if (doserate_start >= next_doserate_time){
            //             current_doserate_update_idx+=1u;
            //             current_doserate = doserate_updates[current_doserate_update_idx].value;
            //             next_doserate_time = doserate_updates[current_doserate_update_idx+1].time;

            //         }  
                    

            //         //this doserate could potentially contribute to the concentration at t
            //         if ((doserate_start < t) && (t<doserate_end + maxtime)){
            //             doserate_end = std::min(t, doserate_end);
            //             //Rcpp::Rcout << "hit! doserate_start : " << doserate_start << " doserate_end: " << doserate_end << "\n";
            //             doserate_duration = doserate_end - doserate_start;
            //             dose = current_doserate * doserate_duration;

                        
            //             //Rcpp::Rcout << "doserate_start: " << doserate_start << " current doserate: " << current_doserate << " dose: " << dose << "\n";
            //             concentration += (dose / vd) * std::exp(-lambda * (t - 0.5 * (doserate_end + doserate_start)));

            //         }
                    
            //         //if the end of this interval is more than maxtime from the current t, it doesn't need to be checked again
            //         if ((t > doserate_end + maxtime)){
            //             vec_doserates_times_begin = j;
            //         }

            //         //dito for the shorter vector of dose changes
            //         if ((t > doserate_updates[doserate_idx_begin+1u].time + maxtime)){
            //             doserate_idx_begin+=1;
            //         }

            //         //break the loop once the starttimes go beyond t
            //         if (doserate_start > t) break;

            //     }
            // }

            //then check contributions from bolus doses
            for (auto j = boluses_begin;j<boluses_end;++j){
                if (((*bolus_times)[j] < t) &&  (t < (*bolus_times)[j]+maxtime)){
                    concentration += ((*boluses)[j] / vd) * std::exp(-lambda * (t - (*bolus_times)[j]));
                }


                //update the loop index if the t is too far off in the future = we don't need to check the next t for this bolus
                if ((*bolus_times)[j]+maxtime<t){
                    boluses_begin = j;
                }

                //break the inner loop if the time is beyond t
                if ((*bolus_times)[j] > t){
                    break;
                }
            }

            (*res)[i] = concentration;
        }
        
    }

};

//a first shot at templated expressions to allow fast parsing
//of sequences of operations like addition and multiplication

//using curiously recurring template pattern

//base expression class
template <typename derived>
struct texpr{
    const derived& derived() const { return static_cast<const derived&> (*this) ;}

    double eval_x1() const { return derived().eval_x1(); }
    double eval_x2() const { return derived().eval_x2(); }
};

//two value struct
struct tstate : public texpr<tstate> {
    double x1; 
    double x2;
    //constructor from doubles
    tstate(double x1, double x2) : x1(x1), x2(x2) {}
    //constructor from another expr
    template <typename E>
    tstate(const texpr<E> &expr) {
        x1 = expr.derived().eval_x1();
        x2 = expr.derived().eval_x2();
    }

    double eval_x1() const  { return x1; }
    double eval_x2() const  { return x2; }

    template <typename E>
    tstate& operator+=(const texpr<E> & expr){
        x1 += expr.derived().eval_x1();
        x2 += expr.deriver().eval_x2();
        return *this;
    }

    template <typename E>
    tstate& operator-=(const texpr<E> & expr){
        x1 -= expr.derived().eval_x1();
        x2 -= expr.deriver().eval_x2();
        return *this;
    }

    template <typename E>
    tstate& operator*=(const texpr<E> & expr){
        x1 *= expr.derived().eval_x1();
        x2 *= expr.deriver().eval_x2();
        return *this;
    }

    template <typename E>
    tstate& operator/=(const texpr<E> & expr){
        x1 /= expr.derived().eval_x1();
        x2 /= expr.deriver().eval_x2();
        return *this;
    }

    tstate& operator+=(double x){
        x1 += x;
        x2 += x;
        return *this;
    }

    tstate& operator-=(double x){
        x1 -= x;
        x2 -= x;
        return *this;
    }

    tstate& operator*=(double x){
        x1 *= x;
        x2 *= x;
        return *this;
    }

    tstate& operator/=(double x){
        x1 /= x;
        x2 /= x;
        return *this;
    }

}

//one value scalar
struct tscalar : public texpr<tstate> {
    double x;
    //constructor from double
    tscalar(double x) : x(x){}

    double eval_x1() const  { return x; }
    double eval_x2() const  { return x; }
}


//addition expression
template <typename LHS, typename RHS>
struct add_texpr : public texpr<add_texpr<LHS, RHS>> {
    const LHS& lhs;
    const RHS& rhs;

    add_texpr(const LHS& lhs, const RHS& rhs) : lhs(lhs), rhs(rhs) {}

    double eval_x1() const { return lhs.eval_x1() + rhs.eval_x1(); }
    double eval_x2() const { return lhs.eval_x2() + rhs.eval_x2(); }
};

//subtraction expression
template <typename LHS, typename RHS>
struct sub_texpr : public texpr<sub_texpr<LHS, RHS>> {
    const LHS& lhs;
    const RHS& rhs;

    sub_texpr(const LHS& lhs, const RHS& rhs) : lhs(lhs), rhs(rhs) {}

    double eval_x1() const { return lhs.eval_x1() - rhs.eval_x1(); }
    double eval_x2() const { return lhs.eval_x2() - rhs.eval_x2(); }
};

//multiplication expression
template <typename LHS, typename RHS>
struct mult_texpr : public texpr<mult_texpr<LHS, RHS>> {
    const LHS& lhs;
    const RHS& rhs;

    mult_texpr(const LHS& lhs, const RHS& rhs) : lhs(lhs), rhs(rhs) {}

    double eval_x1() const { return lhs.eval_x1() * rhs.eval_x1(); }
    double eval_x2() const { return lhs.eval_x2() * rhs.eval_x2(); }
};

//division expression
template <typename LHS, typename RHS>
struct div_texpr : public texpr<div_texpr<LHS, RHS>> {
    const LHS& lhs;
    const RHS& rhs;

    div_texpr(const LHS& lhs, const RHS& rhs) : lhs(lhs), rhs(rhs) {}

    double eval_x1() const { return lhs.eval_x1() / rhs.eval_x1(); }
    double eval_x2() const { return lhs.eval_x2() / rhs.eval_x2(); }
};

//operator overloading for +
template <typename LHS, typename RHS>
add_texpr<LHS, RHS> operator +(const texpr<LHS> &lhs, const texpr<RHS> &rhs){
    return add_texpr<LHS, RHS>(lhs.derived(), rhs.derived());
}

//operator overloading for -
template <typename LHS, typename RHS>
sub_texpr<LHS, RHS> operator -(const texpr<LHS> &lhs, const texpr<RHS> &rhs){
    return sub_texpr<LHS, RHS>(lhs.derived(), rhs.derived());
}

//operator overloading for *
template <typename LHS, typename RHS>
mult_texpr<LHS, RHS> operator *(const texpr<LHS> &lhs, const texpr<RHS> &rhs){
    return mult_texpr<LHS, RHS>(lhs.derived(), rhs.derived());

//overload for scalar multiplication (double on the left)
template <typename RHS>
mult_texpr<tscalar, RHS> operator*(double scalar, const texpr<RHS>& rhs) {
    return mult_texpr<tscalar, RHS>(tscalar(scalar), rhs.derived());
}

//overload for scalar multiplication (double on the left)
template <typename LHS>
mult_texpr<LHS, tscalar> operator*(const texpr<LHS>& lhs, double scalar) {
    return mult_texpr<LHS, tscalar>(lhs.derived(), tscalar(scalar));
}
}

//operator overloading for /
template <typename LHS, typename RHS>
div_texpr<LHS, RHS> operator /(const texpr<LHS> &lhs, const texpr<RHS> &rhs){
    return div_texpr<LHS, RHS>(lhs.derived(), rhs.derived());
}





struct pk_params{
    double v1;
    double v2;
    double k12;
    double k21;
    double ke1;
    double ke2;
};

//calculate derivatives
tstate calculate_derivatives(const tstate& s, const pk_params& params, double dt, double dose, double doserate) {
    return tstate(
        -(params.k12 / params.v1) * s.x1 + (params.k21 / params.v1) * s.x2 - params.ke1 * (s.x1 / params.v1) + (dose + doserate * dt) / params.v1,
        (params.k12 / params.v2) * s.x1 - (params.k21 / params.v2) * s.x2 - params.ke2 * (s.x2 / params.v2)
    );
}

//Runge-kutta 4 simulation
//use this for speed to skip the division
static constexpr double one_sixth = 1.0 / 6.0;
tstate simulate_rk4(const tstate& s, const pk_params& params, double dt, double dose, double doserate) {
    tstate k1 = calculate_derivatives(s, params, dt, dose, doserate);
    tstate k2 = calculate_derivatives(s + k1 * (0.5 * dt), params, 0.5 * dt, dose, doserate);
    tstate k3 = calculate_derivatives(s + k2 * (0.5 * dt), params, 0.5 * dt, dose, doserate);
    tstate k4 = calculate_derivatives(s + k3 * dt, params, dt, dose, doserate);

    return s + (k1 + k2 * 2 + k3 * 2 + k4) * (dt * one_sixth);
}



// a function that will return pharmacokinetic estimation of drug concentrations at a given timepoint
template <typename Tv>
std::function<double(double)> conc_at_t(
    const Tv &times,
    const Tv &doserates,
    double half_life,
    double Vd,
    double max_halflives)
{
    // calculate the elimination rate constant lambda from half-life
    double lambda = std::log(2) / half_life;
    double max_time_after = half_life * max_halflives;

    // Return a lambda function that computes concentration at any time point t
    return [&times, &doserates, lambda, Vd, max_time_after](double t) -> double
    {
        double concentration = 0.0;
        auto n = times.size() - 1u;
        double start_time = times[0];
        double end_time = times[0];
        double rate = doserates[0];
        double dose = 0;

        //Rcpp::Rcout << "vd: " << Vd <<  " max_time_after: " << max_time_after << " lambda: " << lambda << "\n";


        // Calculate concentration contribution from each infusion interval
        for (auto i = 0u; i < n; ++i)
        {
            rate = doserates[i]; // infusion rate for the interval [x[i], x[i+1]]
            start_time = times[i];
            end_time = times[i + 1];

            // sum those that are lower than the maximum half-lives
            if ((t >= start_time) & (t < end_time + max_time_after))
            {
                // Calculate the effective end time for this interval
                // end_time = (t < end_time) ? t : end_time;
                end_time = std::min(t, end_time);

                // Dose administered in this interval
                dose = rate * (end_time - start_time);

                //Rcpp::Rcout << "start_time: " << start_time << " rate: " << rate << " dose: " << dose << "\n";

                // Add concentration contribution from this interval
                concentration += (dose / Vd) * std::exp(-lambda * (t - 0.5 * (start_time + end_time)));
            }
            // break the loop once the times go beyond t
            if (end_time > t)
                break;
        }

        return concentration;
    };
}

//check for NA values before proceeding instead of crashing the entire thing
bool contains_na(const Rcpp::NumericVector& vec) {
    return Rcpp::is_true(Rcpp::any(Rcpp::is_na(vec)));
}

//[[Rcpp::export]]
Rcpp::NumericVector cpp_conc_at_t(Rcpp::NumericVector times,
                                  Rcpp::NumericVector dosetimes,
                                  Rcpp::NumericVector doserates,
                                  double Vd,
                                  double halflife,
                                  double max_halflives = 8)
{
    if(contains_na(times)){
        Rcpp::stop("times cannot contain NA values");
    }

    if(contains_na(dosetimes)){
        Rcpp::stop("dosetimes cannot contain NA values");
    }

    if(contains_na(doserates)){
        Rcpp::stop("doserates cannot contain NA values");
    }
    
    // Setup the concentration function for this patient
    auto calculate_concentration = conc_at_t<Rcpp::NumericVector>(dosetimes, doserates, halflife, Vd, max_halflives);



    std::size_t n = times.size();
    // the out vector
    Rcpp::NumericVector res(n, Rcpp::NumericVector::get_na());

    for (auto i = 0u; i < n; ++i)
    {
        //Rcpp::Rcout << "i: " << i << "\n";
        res[i] = calculate_concentration(times[i]);
    }

    return res;
}



Rcpp::NumericVector cpp_conc_at_t_discrete_old(Rcpp::NumericVector groups,
    Rcpp::NumericVector times,
    Rcpp::NumericVector vd_groups,
    Rcpp::NumericVector vds,
    Rcpp::NumericVector halflife_groups,
    Rcpp::NumericVector halflives,  
    double max_interval_length = 1800,
    double max_halflives = 8,
    int nThreads = 0,
    int threadThreshold = 1,
    int threadMultiplier = 4,
    int minThreads = 1,
    Rcpp::Nullable<Rcpp::NumericVector> DoseRateGroups = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> DoseRateTimes = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> DoseRates = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> BolusGroups = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> BolusTimes = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> BolusDoses = R_NilValue){
    
    

    //set upp the c++ vectors for the doserates and boluses
    Rcpp::NumericVector doserate_groups;
    Rcpp::NumericVector doserate_times;
    Rcpp::NumericVector doserate_rates;
    Rcpp::NumericVector bolus_groups;
    Rcpp::NumericVector bolus_times;
    Rcpp::NumericVector bolus_doses;

    if (DoseRateGroups.isNull() & BolusGroups.isNull()){
        Rcpp::stop("Either DoseRateGroups or BolusGroups must be provided");
    }

    Rcpp::Rcout << "trying to fetch vectors\n";

    // Only assign if DoseRateGroups is not null
    if (!DoseRateGroups.isNull()) {
        doserate_groups = Rcpp::as<Rcpp::NumericVector>(DoseRateGroups);
        if (contains_na(doserate_groups)) {
            Rcpp::stop("Error: DoseRateGroups contains NA values.\n");
        }
        Rcpp::Rcout << "Fetched DoseRateGroups\n";
    }

    // Only assign if DoseRateTimes is not null
    if (!DoseRateTimes.isNull()) {
        doserate_times = Rcpp::as<Rcpp::NumericVector>(DoseRateTimes);
        if (contains_na(doserate_times)) {
            Rcpp::stop("Error: DoseRateTimes contains NA values.\n");
        }
        Rcpp::Rcout << "Fetched DoseRateTimes\n";
    }

    // Only assign if DoseRateValues is not null
    if (!DoseRates.isNull()) {

        doserate_rates = Rcpp::as<Rcpp::NumericVector>(DoseRates);
        if (contains_na(doserate_rates)) {
            Rcpp::stop("Error: DoseRates contains NA values.\n");
        }
        Rcpp::Rcout << "Fetched DoseRates\n";
    }

    // Similar handling for DoseGroups, DoseTimes, and DoseValues
    if (!BolusGroups.isNull()) {
        bolus_groups = Rcpp::as<Rcpp::NumericVector>(BolusGroups);
        if (contains_na(bolus_groups)) {
            Rcpp::stop("Error: BolusGroups contains NA values.\n");
        }
        Rcpp::Rcout << "Fetched BolusGroups\n";
    }

    if (!BolusTimes.isNull()) {
        bolus_times = Rcpp::as<Rcpp::NumericVector>(BolusTimes);
        if (contains_na(bolus_times)) {
            Rcpp::stop("Error: BolusTimes contains NA values.\n");
        }
        Rcpp::Rcout << "Fetched BolusTimes\n";
    }

    if (!BolusDoses.isNull()) {
        bolus_doses = Rcpp::as<Rcpp::NumericVector>(BolusDoses);
        if (contains_na(bolus_doses)) {
            Rcpp::stop("Error: BolusDoses contains NA values.\n");
        }
        Rcpp::Rcout << "Fetched BolusDoses\n";
    }
    Rcpp::Rcout << "fetched any non-null vectors\n";

    //sanity checks on doserate groups, times and values
    if(contains_na(doserate_groups)){
        Rcpp::stop("DoseRateGroups cannot contain NA values");
    }

    if(contains_na(doserate_times)){
        Rcpp::stop("DoseRateTimes cannot contain NA values");
    }

    if(contains_na(doserate_rates)){
        Rcpp::stop("DoseRates cannot contain NA values");
    }

    if (doserate_groups.size() != doserate_times.size()){
        Rcpp::stop("DoseRateGroups and DoseRateTimes must be of identical size");
    }

    if (doserate_groups.size() != doserate_rates.size()){
        Rcpp::stop("DoseRateGroups and DoseRates must be of identical size");
    }

    //sanity checks on bolus groups, times and values
    if(contains_na(bolus_groups)){
        Rcpp::stop("BolusGroups cannot contain NA values");
    }

    if(contains_na(bolus_times)){
        Rcpp::stop("BolusTimes cannot contain NA values");
    }

    if(contains_na(bolus_doses)){
        Rcpp::stop("BolusDoses cannot contain NA values");
    }

    if (bolus_groups.size() != bolus_times.size()){
        Rcpp::stop("BolusGroups and BolusTimes must be of identical size");
    }

    if (bolus_groups.size() != bolus_doses.size()){
        Rcpp::stop("BolusGroups and BolusDoses must be of identical size");
    }

    if (max_interval_length <= 0){
        Rcpp::stop("max_interval_length must be larger than 0");
    }

    Rcpp::Rcout << "vector sanity checks 2 passed\n";

    //data verified, let's go
    std::size_t n = groups.size();
    Rcpp::NumericVector res(n);
    //keep the results from the threadpool here
    std::vector<std::future<void>> resvec;

    if (nThreads <= 0){
         nThreads = std::max(1u,std::thread::hardware_concurrency());
    }

    //the scope for the pool
    {
        //initialize the pool
        threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);

        Rcpp::Rcout << "thread pool initialized\n";
        bool has_doserates = false;
        bool has_boluses = false;
        bool has_vd = false;
        bool has_halflife = false;

        //the location of the groups
        std::size_t this_group = groups[0];
        std::size_t this_group_start = 0u;
        std::size_t this_group_end = 0u;

        //the location of the groups within the doserates
        std::size_t doserate_group_start = 0u;
        std::size_t doserate_group_end = 0u;

        //the location of the groups within the doserates
        std::size_t bolus_group_start = 0u;
        std::size_t bolus_group_end = 0u;

        //the location in the vd vector
        std::size_t vd_idx = 0u;

        //the location in the halflife vector
        std::size_t halflife_idx = 0u;

        std::size_t n_doserates = doserate_groups.size();
        std::size_t n_boluses = bolus_groups.size();
        std::size_t n_vds = vd_groups.size();
        std::size_t n_halflives = halflife_groups.size();

        //loop through the group vector
        for (auto i = 0u;i<n;++i){
            
            //we've found a new group or we're at the end of the groups vector
            if ((groups[i] != this_group) | (i == n - 1)) {
				if (i == n - 1) {
					this_group_end = i + 1;
				}
				else {
					this_group_end = i;
				}

                
                Rcpp::Rcout << "Found group: " << this_group << " at " << this_group_start << " to " << this_group_end << "\n";

                has_doserates = false;
                has_boluses = false;
                has_vd = false;
                has_halflife = false;

                //find this group in the doserates
                for (auto j = doserate_group_start;j<n_doserates;++j){
                    if ((!has_doserates) && (doserate_groups[j] == this_group)){
                        doserate_group_start = j;
                        has_doserates = true;
                    }

                    if (has_doserates && ((doserate_groups[j] != this_group) | (j == n_doserates - 1))){
                        if (j == n_doserates -1){
                            doserate_group_end = j+1;
                        } else {
                            doserate_group_end = j;
                        }

                        Rcpp::Rcout << "Found doserates: " << doserate_groups[doserate_group_start] << " at " << doserate_group_start << " to " << doserate_group_end << "\n";
                        break;
                    }
                }
                //and in the boluses
                for (auto j = bolus_group_start;j<n_boluses;++j){
                    if ((!has_boluses) && (bolus_groups[j] == this_group)){
                        bolus_group_start = j;
                        has_boluses = true;
                    }

                    if (has_boluses && ((bolus_groups[j] != this_group) | (j == n_boluses - 1))){
                        if (j == n_boluses -1){
                            bolus_group_end = j+1;
                        } else {
                            bolus_group_end = j;
                        }

                        Rcpp::Rcout << "Found boluses: " << bolus_groups[bolus_group_start] << " at " << bolus_group_start << " to " << bolus_group_end << "\n";
                        break;
                    }
                }

                //find the index in the vd
                for (auto j = vd_idx;j<n_vds;++j){
                    if ((!has_vd) && (vd_groups[j] == this_group)){
                        vd_idx = j;
                        has_vd = true;
                        Rcpp::Rcout << "Found vd: " << vds[vd_idx] << " at " << vd_idx << "\n";
                        break;
                    }
                }

                //find the index in the halflives
                for (auto j = halflife_idx;j<n_halflives;++j){
                    if ((!has_halflife) && (halflife_groups[j] == this_group)){
                        halflife_idx = j;
                        has_halflife = true;
                        Rcpp::Rcout << "Found halflife: " << halflives[halflife_idx] << " at " << halflife_idx << "\n";
                        break;
                    }
                }

                //if we have either boluses or doserates, and a vd and halflife : let's calculate
                if ((has_boluses | has_doserates) && has_vd && has_halflife){
                    auto dr_begin = has_doserates ? doserate_group_start : 0u;
                    auto dr_end = has_doserates ? doserate_group_end : 0u;
                    auto b_begin = has_boluses ? bolus_group_start : 0u;
                    auto b_end = has_boluses ? bolus_group_end : 0u;

                    resvec.push_back(std::move(pool.add_task(
                        concentration_at_times<Rcpp::NumericVector> (&times,
                                                &doserate_rates,
                                                &doserate_times,
                                                &bolus_doses,
                                                &bolus_times,
                                                &res, 
                                                this_group_start,
                                                this_group_end,
                                                dr_begin,
                                                dr_end,
                                                b_begin,
                                                b_end,
                                                vds[vd_idx],
                                                halflives[halflife_idx],
                                                max_halflives,
                                                max_interval_length)
                    )));                    
                }
                if (has_doserates){
                    doserate_group_start = doserate_group_end;
                }

                if (has_boluses){
                    bolus_group_start = bolus_group_end;
                }

                this_group = groups[i];
				this_group_start = i;               
            }
        }

    }


    for (auto &r : resvec){
		r.get();
	}

    return res;
}

//[[Rcpp::export]]
Rcpp::NumericVector cpp_conc_at_times_discrete(Rcpp::NumericVector groups,
    Rcpp::NumericVector times,
    Rcpp::NumericVector vd_groups,
    Rcpp::NumericVector vds,
    Rcpp::NumericVector halflife_groups,
    Rcpp::NumericVector halflives,  
    Rcpp::Nullable<Rcpp::NumericVector> DoseRateGroups = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> DoseRateTimes = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> DoseRates = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> BolusGroups = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> BolusTimes = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> BolusDoses = R_NilValue,
    double max_interval_length = 1800,
    double max_halflives = 8,
    int nThreads = 1,
    int threadThreshold = 1,
    int threadMultiplier = 4,
    int minThreads = 1){

        
    //Rcpp::Rcout << "in main function body \n";

    //set upp the c++ vectors for the doserates and boluses
    Rcpp::NumericVector doserate_groups;
    Rcpp::NumericVector doserate_times;
    Rcpp::NumericVector doserate_rates;
    Rcpp::NumericVector bolus_groups;
    Rcpp::NumericVector bolus_times;
    Rcpp::NumericVector bolus_doses;

    if (DoseRateGroups.isNull() & BolusGroups.isNull()){
        Rcpp::stop("Either DoseRateGroups or BolusGroups must be provided");
    }

    //Rcpp::Rcout << "trying to fetch vectors\n";

    // Only assign if DoseRateGroups is not null
    if (!DoseRateGroups.isNull()) {
        doserate_groups = Rcpp::as<Rcpp::NumericVector>(DoseRateGroups);
        if (contains_na(doserate_groups)) {
            Rcpp::stop("Error: DoseRateGroups contains NA values.\n");
        }
        //Rcpp::Rcout << "Fetched DoseRateGroups\n";
    }

    // Only assign if DoseRateTimes is not null
    if (!DoseRateTimes.isNull()) {
        doserate_times = Rcpp::as<Rcpp::NumericVector>(DoseRateTimes);
        if (contains_na(doserate_times)) {
            Rcpp::stop("Error: DoseRateTimes contains NA values.\n");
        }
        //Rcpp::Rcout << "Fetched DoseRateTimes\n";
    }

    // Only assign if DoseRateValues is not null
    if (!DoseRates.isNull()) {

        doserate_rates = Rcpp::as<Rcpp::NumericVector>(DoseRates);
        if (contains_na(doserate_rates)) {
            Rcpp::stop("Error: DoseRates contains NA values.\n");
        }
        //Rcpp::Rcout << "Fetched DoseRates\n";
    }

    // Similar handling for DoseGroups, DoseTimes, and DoseValues
    if (!BolusGroups.isNull()) {
        bolus_groups = Rcpp::as<Rcpp::NumericVector>(BolusGroups);
        if (contains_na(bolus_groups)) {
            Rcpp::stop("Error: BolusGroups contains NA values.\n");
        }
        //Rcpp::Rcout << "Fetched BolusGroups\n";
    }

    if (!BolusTimes.isNull()) {
        bolus_times = Rcpp::as<Rcpp::NumericVector>(BolusTimes);
        if (contains_na(bolus_times)) {
            Rcpp::stop("Error: BolusTimes contains NA values.\n");
        }
        //Rcpp::Rcout << "Fetched BolusTimes\n";
    }

    if (!BolusDoses.isNull()) {
        bolus_doses = Rcpp::as<Rcpp::NumericVector>(BolusDoses);
        if (contains_na(bolus_doses)) {
            Rcpp::stop("Error: BolusDoses contains NA values.\n");
        }
        //Rcpp::Rcout << "Fetched BolusDoses\n";
    }
    //Rcpp::Rcout << "fetched any non-null vectors\n";

    //sanity checks on doserate groups, times and values
    if(contains_na(doserate_groups)){
        Rcpp::stop("DoseRateGroups cannot contain NA values");
    }

    if(contains_na(doserate_times)){
        Rcpp::stop("DoseRateTimes cannot contain NA values");
    }

    if(contains_na(doserate_rates)){
        Rcpp::stop("DoseRates cannot contain NA values");
    }

    if (doserate_groups.size() != doserate_times.size()){
        Rcpp::stop("DoseRateGroups and DoseRateTimes must be of identical size");
    }

    if (doserate_groups.size() != doserate_rates.size()){
        Rcpp::stop("DoseRateGroups and DoseRates must be of identical size");
    }

    //sanity checks on bolus groups, times and values
    if(contains_na(bolus_groups)){
        Rcpp::stop("BolusGroups cannot contain NA values");
    }

    if(contains_na(bolus_times)){
        Rcpp::stop("BolusTimes cannot contain NA values");
    }

    if(contains_na(bolus_doses)){
        Rcpp::stop("BolusDoses cannot contain NA values");
    }

    if (bolus_groups.size() != bolus_times.size()){
        Rcpp::stop("BolusGroups and BolusTimes must be of identical size");
    }

    if (bolus_groups.size() != bolus_doses.size()){
        Rcpp::stop("BolusGroups and BolusDoses must be of identical size");
    }

    
    //sanity checks on groups and times
    if(contains_na(groups)){
        Rcpp::stop("groups cannot contain NA values");
    }

    if(contains_na(times)){
        Rcpp::stop("times cannot contain NA values");
    }

    if (groups.size() != times.size()){
        Rcpp::stop("groups and times must be of identical size");
    }

    //sanity checks on vd groups and times
    if(contains_na(vd_groups)){
        Rcpp::stop("vd_groups cannot contain NA values");
    }

    if(contains_na(vds)){
        Rcpp::stop("vds cannot contain NA values");
    }

    if (vd_groups.size() != vds.size()){
        Rcpp::stop("vd_groups and vds must be of identical size");
    }

    //sanity checks on halflife groups and times
    if(contains_na(halflife_groups)){
        Rcpp::stop("halflife_groups cannot contain NA values");
    }

    if(contains_na(halflives)){
        Rcpp::stop("halflives cannot contain NA values");
    }

    if (halflife_groups.size() != halflives.size()){
        Rcpp::stop("halflife_groups and halflives must be of identical size");
    }

    //Rcpp::Rcout << "vector length sanity check passed\n";

    //Rcpp::Rcout << bolus_groups.size() << "\n";
    //Rcpp::Rcout << doserate_groups.size() << "\n";

    //data verified, let's go
    std::size_t n = groups.size();
    Rcpp::NumericVector res(n);
    //keep the results from the threadpool here
    std::vector<std::future<void>> resvec;

    if (nThreads <= 0){
         nThreads = std::max(1u,std::thread::hardware_concurrency());
    }

    //Rcpp::Rcout << "out and future vectors created\n";
    {
        //Rcpp::Rcout << "starting threadpool...\n";
        threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
        //Rcpp::Rcout << "thread pool initialized\n";

        bool has_doserates = false;
        bool has_boluses = false;
        bool has_vd = false;
        bool has_halflife = false;

        //the location of the groups
        auto this_group = groups[0];
        std::size_t this_group_start = 0u;
        std::size_t this_group_end = 0u;

        //the location of the groups within the doserates
        std::size_t doserate_group_start = 0u;
        std::size_t doserate_group_end = 0u;

        //the location of the groups within the doserates
        std::size_t bolus_group_start = 0u;
        std::size_t bolus_group_end = 0u;

        //the location in the vd vector
        std::size_t vd_idx = 0u;

        //the location in the halflife vector
        std::size_t halflife_idx = 0u;

        std::size_t n_doserates = doserate_groups.size();
        std::size_t n_boluses = bolus_groups.size();
        std::size_t n_vds = vd_groups.size();
        std::size_t n_halflives = halflife_groups.size();

        //Rcpp::Rcout << std::ceil(9/2);

        for (auto i = 0u;i<n;++i){
            //std::cout << "i: " << i << "\n";
            //we've found a new group or we're at the end of the groups vector
            if ((groups[i] != this_group) | (i == n - 1)) {
				if (i == n - 1) {
					this_group_end = i + 1;
				}
				else {
					this_group_end = i;
				}

                
                //Rcpp::Rcout << "Found group: " << this_group << " at " << this_group_start << " to " << this_group_end << "\n";

                has_doserates = false;
                has_boluses = false;
                has_vd = false;
                has_halflife = false;

                //find this group in the doserates
                for (auto j = doserate_group_start;j<n_doserates;++j){
                    if ((!has_doserates) && (doserate_groups[j] == this_group)){
                        doserate_group_start = j;
                        has_doserates = true;
                    }

                    if (has_doserates && ((doserate_groups[j] != this_group) | (j == n_doserates - 1))){
                        if (j == n_doserates -1){
                            doserate_group_end = j+1;
                        } else {
                            doserate_group_end = j;
                        }

                        //Rcpp::Rcout << "Found doserates: " << doserate_groups[doserate_group_start] << " at " << doserate_group_start << " to " << doserate_group_end << "\n";
                        break;
                    }
                }

                //and in the boluses
                for (auto j = bolus_group_start;j<n_boluses;++j){
                    if ((!has_boluses) && (bolus_groups[j] == this_group)){
                        bolus_group_start = j;
                        has_boluses = true;
                    }

                    if (has_boluses && ((bolus_groups[j] != this_group) | (j == n_boluses - 1))){
                        if (j == n_boluses -1){
                            bolus_group_end = j+1;
                        } else {
                            bolus_group_end = j;
                        }

                        //Rcpp::Rcout << "Found boluses: " << bolus_groups[bolus_group_start] << " at " << bolus_group_start << " to " << bolus_group_end << "\n";
                        break;
                    }
                }

                //find the index in the vd
                for (auto j = vd_idx;j<n_vds;++j){
                    if ((!has_vd) && (vd_groups[j] == this_group)){
                        vd_idx = j;
                        has_vd = true;
                        //Rcpp::Rcout << "Found vd: " << vds[vd_idx] << " at " << vd_idx << "\n";
                        break;
                    }
                }

                //find the index in the halflives
                for (auto j = halflife_idx;j<n_halflives;++j){
                    if ((!has_halflife) && (halflife_groups[j] == this_group)){
                        halflife_idx = j;
                        has_halflife = true;
                        //Rcpp::Rcout << "Found halflife: " << halflives[halflife_idx] << " at " << halflife_idx << "\n";
                        break;
                    }
                }
                //if we have either boluses or doserates, and a vd and halflife : let's calculate
                if ((has_boluses | has_doserates) && has_vd && has_halflife){
                    auto dr_begin = has_doserates ? doserate_group_start : 0u;
                    auto dr_end = has_doserates ? doserate_group_end : 0u;
                    auto b_begin = has_boluses ? bolus_group_start : 0u;
                    auto b_end = has_boluses ? bolus_group_end : 0u;
                    //Rcpp::Rcout << "pushing back task with: dr_begin " << dr_begin << " dr_end " << dr_end << " b_begin " << b_begin << " b_end " << b_end << "\n";
                    resvec.push_back(std::move(pool.add_task(
                        concentration_at_times_2<Rcpp::NumericVector> (&times,
                                                &doserate_rates,
                                                &doserate_times,
                                                &bolus_doses,
                                                &bolus_times,
                                                &res, 
                                                this_group_start,
                                                this_group_end,
                                                dr_begin,
                                                dr_end,
                                                b_begin,
                                                b_end,
                                                vds[vd_idx],
                                                halflives[halflife_idx],
                                                max_halflives,
                                                max_interval_length)
                    )));                    
                }

                if (has_doserates){
                    doserate_group_start = doserate_group_end;
                }

                if (has_boluses){
                    bolus_group_start = bolus_group_end;
                }


                this_group = groups[i];
                this_group_start = this_group_end;
                //Rcpp::Rcout << "Group finished\n";
            }
        }
    }
    

    //Rcpp::Rcout << "Threadpool ended, getting results\n";
    for (auto &r : resvec){
		r.get();
	}
    //Rcpp::Rcout << "End of main function\n";

    return res;


}