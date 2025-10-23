
#include <Rcpp.h>
//[[Rcpp::plugins("cpp14")]]

#include <algorithm>
#include <iostream>

#include <vector>
#include <queue>
#include <deque>

#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <future>

#include <string>
#include <chrono>
//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>



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

		for (auto &t : threads) {
			t.join();
		}
	}
};


//2022-02-03
//a distance matrix class based on vectors and modulo operations, with bit shift modulo to speed things up... let's see how it goes!
template<typename Tit>
class dm_vec {
private:
	std::vector<double> coeffs;
	std::vector<Tit> starts;

	uint32_t size, dim;
	uint32_t capacity;
	uint32_t k;
	std::vector<std::vector<double>> mat;

	//uint32_t integers
	uint32_t window_begin = 0;
	uint32_t window_end = 0;

	//some vectors to avoid having to allocate memory all the time
	//the vector for calculating k-th distance
	std::vector<double> k_func_distances;
	//vectors for brute forcing back and forth
	std::vector<double> backfill;
	std::vector<double> frontfill;

	//bit shift hacking!
	//this calculates the i % capacity, if capacity is a power of two
	inline uint32_t im(uint32_t i) {
		return i & (capacity - 1);
	}

	//bit shift hacking 2, find the smallest p larger or equal to n where p is a power of two
	uint32_t nextPowerOf2(uint32_t n)
	{
		uint32_t p = 1;
		if (n && !(n & (n - 1)))
			return n;

		while (p < n)
			p <<= 1;

		return p;
	}

	//Weighted euclidian distance
	//this actually calculates squared distance, to avoid having to call sqrt all the time.
	inline double getdist(uint32_t p1, uint32_t p2) {
		double res = 0.0;
		for (auto i = 0u; i < dim; ++i) {
			res += std::pow((*(starts[i] + p1)) - (*(starts[i] + p2)), 2.0) * coeffs[i];
		}
		//res = std::sqrt(res);
		return res;
	}

	//weighted 1d distance
	//dito - calculate squared distance in 1 dimension
	inline double get_1d_dist(uint32_t p1, uint32_t p2, uint32_t dim) {
		double res = std::pow((*(starts[dim] + p1)) - (*(starts[dim] + p2)), 2.0) * coeffs[dim];
		//res = std::sqrt(res);
		return res;
	}

	//get the size of our current window
	inline uint32_t matrix_size() {
		return window_end - window_begin;
	}

	//'remove' the first element of the matrix, i.e. just change the index of the start point
	inline void distmatrix_pop_front() {
		window_begin++;
	}
	//'remove' the last element of the matrix, i.e. just change the index of the end point
	inline void distmatrix_pop_back() {
		window_end--;
	}

	//increase the size of the matrix in the forward direction, i.e. adds the element at window_end
	void distmatrix_push_back() {
		//remove the first element of the matrix if we're at full capacity
		while (matrix_size() >= capacity) {
			distmatrix_pop_front();
		}

		//update the matrix
		//ir is the row in the matrix for our new point
		uint32_t ir = im(window_end);
		//ii is the column
		uint32_t ii = 0u;
		for (auto i = window_begin; i < window_end; ++i) {
			ii = im(i);
			//calculate the distance between the new point (window_end) and the current old point (i)
			//and add it to the row for our new point (ir = im(window_end))
			mat[ir][ii] = getdist(window_end, i);
			//make the matrix square!
			mat[ii][ir] = mat[ir][ii];
		}
		//add our self distance, i.e. 0
		mat[ir][ir] = 0;
		//update window_end
		window_end++;
	}
	//increase the size of the matrix in the forward direction, i.e. adds the element at window_begin -1
	void distmatrix_push_front() {
		//remove the last element of the matrix if we're at full capacity
		while (matrix_size() >= capacity) {
			distmatrix_pop_back();
		}

		//update the matrix
		//ir is the current row in the matrix...
		uint32_t ir = im(window_begin - 1);
		//ii is the column
		uint32_t ii = 0;
		for (auto i = window_begin; i < window_end; ++i) {
			ii = im(i);
			//calculate the distance between the new point (window_begin -1) and the current old point (i)
			//and add it to the row for our new point (iw = im(window_begin -1))
			mat[ir][ii] = getdist(window_begin - 1, i);
			//make the matrix square!
			mat[ii][ir] = mat[ir][ii];
		}
		//add our self distance, i.e. 0
		mat[ir][ir] = 0;
		//update window_end
		window_begin--;
	}

	//copy a row from the matrix into a vector, in the correct order
	template<typename Ta>
	inline void matrix_to_container(Ta& target, uint32_t index) {
		auto mat_size = matrix_size();
		target.resize(mat_size);
		auto ir = im(index);
		for (auto i = 0u; i < mat_size; ++i) {
			target[i] = mat[ir][im(i + window_begin)];
		}
	}

	//better resize target first!!
	template<typename Ta>
	inline void matrix_to_container2(Ta& target, uint32_t index) {
		//get the matrix size
		auto mat_size = matrix_size();
		//copy our row
		auto ir = im(index);
		for (auto i = 0u; i < mat_size; ++i) {
			target[i] = mat[ir][im(i + window_begin)];
		}
	}

	//get the kth element from a row in the matrix
	template<typename Ta>
	inline double get_k_dist_matrix(Ta& target, uint32_t k, uint32_t index) {
		auto mat_size = matrix_size();
		target.resize(mat_size);
		//copy the matrix row to start of container
		matrix_to_container2(target, index);
		std::nth_element(target.begin(), target.begin() + k, target.end());
		return target[k];
	}

	//get the kth element from a row in the matrix and a container
	template<typename Ta, typename Tb>
	inline double get_k_dist_matrix(Ta& target, uint32_t k, uint32_t index, const Tb& b) {
		auto mat_size = matrix_size();
		target.resize(mat_size + b.size());
		//copy the matrix row
		matrix_to_container2(target, index);
		//copy b
		std::copy(b.begin(), b.end(), target.begin() + mat_size);
		std::nth_element(target.begin(), target.begin() + k, target.end());
		return target[k];
	}

	//get the kth element from a row in the matrix and two containers
	template<typename Ta, typename Tb, typename Tc>
	inline double get_k_dist_matrix(Ta& target, uint32_t k, uint32_t index, const Tb& b, const Tc& c) {
		auto mat_size = matrix_size();
		target.resize(mat_size + b.size() + c.size());
		//copy the matrix row
		matrix_to_container2(target, index);
		//copy b
		std::copy(b.begin(), b.end(), target.begin() + mat_size);
		//copy c
		std::copy(c.begin(), c.end(), target.begin() + mat_size + b.size());
		std::nth_element(target.begin(), target.begin() + k, target.end());
		return target[k];
	}

	//breaking this code out to a separate function to avoid duplicating it in the various distance calls
	inline void shrink_matrix_to_fit(uint32_t index, uint32_t k, double k_dist) {
		while ((matrix_size() > k + 1) && window_begin < index && (get_1d_dist(window_begin, index, 0) > k_dist)) {
			distmatrix_pop_front();
		}
	}

	//breaking away the part where the matrix is updated and the k_distance is calculated
	//note, it will change the backfill and frontfill vectors
	inline double get_k_dist(uint32_t index, uint32_t k) {
		//enlarge the matrix forwards if it's too small
		while ((matrix_size() < k + 1) && (window_end < size)) {
			distmatrix_push_back();
		}
		//enlarge it backwards if it's still too small
		while ((matrix_size() < k + 1) && (window_begin > 0)) {
			distmatrix_push_front();
		}

		//if we've somehow managed to not insert our own row into the matrix, force it forwards until it's there
		while (window_end <= index) {
			distmatrix_push_back();
		}
		//dito, in the reverse direction
		while (window_begin > index && window_begin > 0) {
			distmatrix_push_front();
		}

		//a guess of k_distance based on the current distance matrix
		double k_dist = get_k_dist_matrix(k_func_distances, k, index);

		//try to enlarge the matrix forward, if we're not at the end of our view of the vector
		while (window_end < size) {
			//if the v0-distance is <= k-dist, and not at capacity -> enlarge forward
			if ((get_1d_dist(index, window_end, 0) <= k_dist) && (matrix_size() < capacity)) {
				distmatrix_push_back();
				//recalculate the k distance 
				k_dist = get_k_dist_matrix(k_func_distances, k, index);
			}
			else break;
		}

		//do a brute force seach backwards if there might be points earlier that are closer than our current kth-distance
		//this is instead of enlarging the matrix backwards just for a single outlier...
		backfill.clear();
		//now with uint32_t integers!
		uint32_t backfill_offset = 1;
		while (window_begin >= backfill_offset) {
			if (get_1d_dist(index, window_begin - backfill_offset, 0) <= k_dist) {
				backfill.push_back(getdist(index, window_begin - backfill_offset));
				k_dist = get_k_dist_matrix(k_func_distances, k, index, backfill);
				backfill_offset++;
			}
			else break;
		}

		//dito in the forward direction if the matrix is at capacity.
		frontfill.clear();
		uint32_t frontfill_offset = 0;
		while (window_end + frontfill_offset < size) {
			if (get_1d_dist(index, window_end + frontfill_offset, 0) <= k_dist) {
				frontfill.push_back(getdist(index, window_end + frontfill_offset));
				frontfill_offset++;
				k_dist = get_k_dist_matrix(k_func_distances, k, index, backfill, frontfill);
			}
			else break;
		}

		return k_dist;

	}

	//debugging, printing the full matrix
	void print_matrix() {
		std::cout << "Full data matrix:\n";

		for (auto& v : mat) {
			for (auto e : v) {
				std::cout << std::round(e) << " ";
			}
			std::cout << "\n";
		}
	}

	void print_window() {
		for (auto i = window_begin; i < window_end; ++i) {
			auto ir = im(i);
			for (auto j = window_begin; j < window_end; ++j) {
				std::cout << mat[ir][im(j)] << " ";
			}
			std::cout << "\n";
		}
	}

public:
	dm_vec(std::vector<Tit>& starts, std::vector<double>& coeffs, uint32_t size, uint32_t k) {
		this->starts = starts;
		//copy the coeffs
		this->coeffs = coeffs;

		//copy the constants
		this->size = size;
		this->k = k;

		//very computer science-y!
		//+4 is a empiric value found from benchmarking - it gets expensive when the capacity is too close to k, it would seem...
		//this->capacity = nextPowerOf2(k + 4u);
		this->capacity = nextPowerOf2(k + 16u);
		dim = coeffs.size();
		window_begin = 0;
		window_end = 0;

		for (auto i = 0u; i < this->capacity; ++i) {
			std::vector<double> vec(this->capacity);
			mat.push_back(vec);
		}
	}

	//get the k nearest neighbors, with optional sorting of results on either distance or index, and may or may not include ties.
	std::vector<std::pair<uint32_t, double>> knn_distances(uint32_t index, uint32_t k, int sort_results_by = 0, bool include_ties = false, int include_self = 0) {
		//get our k-distance, and our neighborhood will be stored in our corresponding matrix row, plus maybe some points in the 
		//backfill or frontfill vectors - if they are outside the window shown by the matrix.
		double k_dist = get_k_dist(index, k);
		std::vector<std::pair<uint32_t, double>> res;
		res.reserve(k + 4);
		uint32_t ir = im(index);

		//include ourselves if that is requested
		if (include_self == 1) {
			res.emplace_back(index, 0.0);
		}

		//this will keep track of if we're already sorted the out-vector when we removed any ties...
		int sorted_by = 0;
		//only for the first group, to keep this sane
		//add the points from the matrix that are part of our neighborhood
		for (auto i = 0u; i < matrix_size(); ++i) {
			auto ii = im(i + window_begin);
			if (mat[ir][ii] <= k_dist && (window_begin + i != index)) {
				res.emplace_back(std::make_pair<uint32_t, double>(window_begin + i, std::sqrt(mat[ir][ii])));
			}
		}


		//and the point from our backfill that are part of our knn
		for (auto i = 0u; i < backfill.size(); ++i) {
			if (backfill[i] <= k_dist) {
				res.emplace_back(std::make_pair<uint32_t, double>(window_begin - 1 - i, std::sqrt(backfill[i])));
			}
		}

		//and the points from the frontfill
		for (auto i = 0u; i < frontfill.size(); ++i) {
			if (frontfill[i] <= k_dist) {
				res.emplace_back(std::make_pair<uint32_t, double>(window_end + i, std::sqrt(frontfill[i])));
			}
		}


		if (!include_ties) {
			//sort the results on distance
			std::sort(res.begin(), res.end(), [](const auto& p1, const auto& p2) {return p1.second < p2.second; });
			//make sure we only keep the k first ones, plus one more if we include ourselves
			res.resize(k + include_self);
			//so, we've sorted by p.second = distance = 2;
			sorted_by = 2;
		}

		//check if we need to sort the results
		if (sort_results_by > 0 && sort_results_by != sorted_by) {
			switch (sort_results_by) {
			case 1: //sort by index, to make it look more like a row in a sparse distance matrix
				std::sort(res.begin(), res.end(), [](const auto& p1, const auto& p2) {return p1.first < p2.first; });
				break;
			case 2: //sort by distance, probably more useful but we'll see if that is the case in future function calls...
				std::sort(res.begin(), res.end(), [](const auto& p1, const auto& p2) {return p1.second < p2.second; });
				break;
			}
		}

		//shrink our matrix
		shrink_matrix_to_fit(index, k, k_dist);

		return res;
	}

	//the k distance neighborhood query function
	std::pair<double, std::vector<uint32_t>> k_dist_neighborhood(uint32_t index, uint32_t k) {
		//get our k-distance, and our neighborhood will be stored in our corresponding matrix row, plus maybe some points in the 
		//backfill or frontfill vectors - if they are outside the window shown by the matrix.
		double k_dist = get_k_dist(index, k);
		//the vector of the indexes of our closest neighborhood
		std::vector<uint32_t> locations;
		//reserve some extra memory for this vector
		locations.reserve(k + 2);
		//precalculate im(index), i.e. our current row
		uint32_t ir = im(index);

		//add the points in the matrix that are part of our k-distance-neighborhood
		for (auto i = 0u; i < matrix_size(); ++i) {
			auto ii = im(i + window_begin);
			//<= to allow k-neighborhood larger than k if there are ties
			if (mat[ir][ii] <= k_dist && (window_begin + i) != index) {
				locations.push_back(window_begin + i);
			}
		}

		//add the points from the backfill that are part of the neigborhood
		for (auto i = 0u; i < backfill.size(); ++i) {
			//<= to allow k-neighborhood larger than k if there are ties
			if (backfill[i] <= k_dist) {
				locations.push_back(window_begin - 1 - i);
			}
		}

		//and the parts of the frontfill
		for (auto i = 0u; i < frontfill.size(); ++i) {
			//<= to allow k-neighborhood larger than k if there are ties
			if (frontfill[i] <= k_dist) {
				locations.push_back(window_end + i);
			}
		}

		//finally, try to shrink the matrix
		shrink_matrix_to_fit(index, k, k_dist);

		//calculate our k-distance
		k_dist = std::sqrt(k_dist);
		//return k-distance and the indexes of our neighborhood
		return std::make_pair(k_dist, std::move(locations));
	}
};


//2022-04-20
//A class that performs linear regression in one variable on subsets of the data
template<typename Td>
class llr {
private:
	//pointers to the data vectors
	Td* times;
	Td* values;
	Td* res;

	//how to select data window
	uint32_t type = 0;
	//what to output
	uint32_t output = 0;

    //which beta to output
    uint32_t output_beta = 0;

	//parameter for knn
	uint32_t k = 0;
	//weights for distance calculations in knn
	std::vector<double> coeffs;

	uint32_t group_size;
	uint32_t group_start;

	//parameter for time windows
	uint32_t before = 0;
	uint32_t after = 0;

    uint32_t degree = 1;

    //a convenience function to copy values and times form the values-vector into an Eigen matrix
    template <typename Tit>
    void knn_to_matrix(Tit &times_begin, Tit &values_begin, Eigen::MatrixXd &X, Eigen::VectorXd &y, std::vector<std::pair<uint32_t, double>> &knn, uint32_t &knn_size, uint32_t &degree){
        for (auto j = 0u; j < knn_size; ++j) {
			X(j, 1) = *(times_begin + knn[j].first);
			y(j) = *(values_begin + knn[j].first);
		}
        //a secondary, nested loop to fit higher degree polynom forms of x
        for (auto d = 1u; d < degree;++d){
            for (auto j = 0u;j<knn_size;++j){
                X(j, d+1) = std::pow(X(j, 1), d+1);
            }
        }
    }

public:
	llr(Td* times, Td* values, Td* res, std::vector<double>& c, uint32_t group_start, uint32_t group_size, uint32_t type, uint32_t output, uint32_t output_beta, uint32_t degree, uint32_t k, uint32_t before, uint32_t after) {
		//the main data vectors
		this->times = times;
		this->values = values;
		this->res = res;

        this->degree = degree;

		//the parameters for this group
		this->group_start = group_start;
		this->group_size = group_size;

		//the parameters for type and output
		this->type = type;
		this->output = output;
        //which beta coefficient to output, if requested
        this->output_beta = std::max(0u, std::min(degree, output_beta));

		//the parameters for time window - for both prefilter and the regression...!
		this->before = before;
		this->after = after;

		//the parameters for knn
		//note that it is impossible to find more than size - 1 neighbors...
		this->k = std::max(0u, std::min(group_size - 1, k));
		//copy the coeffs
		coeffs.resize(c.size());
		std::copy(c.begin(), c.end(), coeffs.begin());
	}

	auto operator()()->void {
		//set up iterators to pass to the knn-method...
		auto values_begin = values->begin() + group_start;
		auto times_begin = times->begin() + group_start;

		//knn based neighborhood
		if (type == 0) {
			//use a typedef to the iterator to make it less wordy...
			typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::iterator ITERATOR_TYPE;
			std::vector<ITERATOR_TYPE> starts;
			starts.push_back(times_begin);
			starts.push_back(values_begin);

			//the distance matrix
			dm_vec<ITERATOR_TYPE> distmatrix(starts, coeffs, group_size, k);

			//reserving memory for our vector k distance neighborhood
			std::vector<std::pair<uint32_t, double>> knn;

			//storing the size of the knn-neighborhood...
			uint32_t knn_size = k + 1;

			//these two can be initialized at once since they will be needed for all output
			//X = the predictor matrix
			Eigen::MatrixXd X(knn_size, degree+1); // knn_size x (1+degree) matrix
			for (auto i = 0u; i < knn_size; ++i) {
				//set our intercept term to one for all elements...
				X(i, 0) = 1.0;
			}
			//y = target vector
			Eigen::VectorXd y(knn_size);

            //for efficiency, this is written to be similar to a unique function per requested output, with unique scopes and variables and so on
            //it does mean that the knn-code is repeated a bit, though...
			switch (output) {
				case 0: //standardized residuals
				{
					Eigen::VectorXd residuals(knn_size);
					Eigen::MatrixXd H(knn_size, knn_size);
					double RSE = 0;
                    double DOF = knn_size - degree - 1;

					for (auto i = 0u; i < group_size; ++i) {
						//parameters : don't sort, don't include ties (which means it will be sorted anyway), and include ourselves (at the 0 index!)
						knn = std::move(distmatrix.knn_distances(i, k, 0, false, 1));
						//set matrix and vector elements to our knn
                        knn_to_matrix(times_begin, values_begin, X, y, knn, knn_size, degree);
						//calculate residuals
						residuals = y - X * (X.transpose() * X).inverse() * X.transpose() * y;
						//calcaulte hat matrix
						H = X * (X.transpose() * X).inverse() * X.transpose();
						//calculate Residual Standard Error
						RSE = std::sqrt(residuals.squaredNorm() / DOF);
						(*res)[group_start + i] = residuals(0) / (RSE * std::sqrt(1.0 - H(0, 0)));
					}
				}
				break;

				case 1: //yhat
				{
					Eigen::VectorXd yhat(knn_size);
					for (auto i = 0u; i < group_size; ++i) {
						//parameters : don't sort, don't include ties (which means it will be sorted anyway), and include ourselves (at the 0 index!)
						knn = distmatrix.knn_distances(i, k, 0, false, 1);
						//set matrix and vector elements to our knn
						knn_to_matrix(times_begin, values_begin, X, y, knn, knn_size, degree);
						//calculate residuals
						yhat = X * (X.transpose() * X).inverse() * X.transpose() * y;
						//return yhat for ourselves
						(*res)[group_start + i] = yhat(0);
					}
				}
				break;

				case 2: //R2
				{
					Eigen::MatrixXd H(knn_size, knn_size);
					Eigen::MatrixXd I = Eigen::MatrixXd::Identity(knn_size, knn_size);
					Eigen::VectorXd ivec = Eigen::VectorXd::Ones(knn_size);
					Eigen::MatrixXd M = ivec * (ivec.transpose() * ivec).inverse() * ivec.transpose();

					for (auto i = 0u; i < group_size; ++i) {
						//parameters : don't sort, don't include ties (which means it will be sorted anyway), and include ourselves (at the 0 index!)
						knn = distmatrix.knn_distances(i, k, 0, false, 1);
						//set matrix and vector elements to our knn
						knn_to_matrix(times_begin, values_begin, X, y, knn, knn_size, degree);

						//calcaulte hat matrix
						H = X * (X.transpose() * X).inverse() * X.transpose();

						//R2 in matrix notation
						(*res)[group_start + i] = 1.0 - ((y.transpose() * (I - H) * y) / (y.transpose() * (I - M) * y))(0,0);
					}
				}
				break;

				case 3: //beta-coefficient
				{
					Eigen::VectorXd betahat(degree+1);

					for (auto i = 0u; i < group_size; ++i) {
						//parameters : don't sort, don't include ties (which means it will be sorted anyway), and include ourselves (at the 0 index!)
						knn = distmatrix.knn_distances(i, k, 0, false, 1);
						//set matrix and vector elements to our knn
						knn_to_matrix(times_begin, values_begin, X, y, knn, knn_size, degree);

						//calcaulte OLS solution
						betahat = (X.transpose() * X).inverse() * X.transpose() * y;

						//beta 0
						(*res)[group_start + i] = betahat(output_beta);
					}
				}
				break;
			}
		} else if (type == 4){
			//use a typedef to the iterator to make it less wordy...
			typedef typename decltype(typename std::remove_pointer < decltype(values) > ::type())::iterator ITERATOR_TYPE;
			std::vector<ITERATOR_TYPE> starts;
			starts.push_back(times_begin);
			starts.push_back(values_begin);

			//the distance matrix
			dm_vec<ITERATOR_TYPE> distmatrix(starts, coeffs, group_size, k);

			//reserving memory for our vector k distance neighborhood
			std::vector<std::pair<uint32_t, double>> knn;

			//storing the size of the knn-neighborhood...
			uint32_t knn_size = k;

			//these two can be initialized at once since they will be needed for all output
			//X = the predictor matrix
			Eigen::MatrixXd X(knn_size, degree+1); // knn_size x (1+degree) matrix
			for (auto i = 0u; i < knn_size; ++i) {
				//set our intercept term to one for all elements...
				X(i, 0) = 1.0;
			}
			//y = target vector
			Eigen::VectorXd y(knn_size);

			//self-matrix
			Eigen::MatrixXd Xs(1, degree+1);
			Xs(0,0) = 1.0;

            //for efficiency, this is written to be similar to a unique function per requested output, with unique scopes and variables and so on
            //it does mean that the knn-code is repeated a bit, though...
			switch (output) {
				case 0: //standardized residuals, here defined as residual(target) / std.dev (residuals)
				{
					Eigen::VectorXd residuals(knn_size);
					double stddev;
					double error;

					for (auto i = 0u; i < group_size; ++i) {
						//set self matrix
						for (auto j =0u;j<degree;++j){
							Xs(0,j+1) = std::pow(*(times_begin+i), j+1.0);
						}
						//parameters : don't sort, don't include ties (which means it will be sorted anyway), and don't include ourselves (at the 0 index!)
						knn = std::move(distmatrix.knn_distances(i, k, 0, false, 0));
						//set matrix and vector elements to our knn
                        knn_to_matrix(times_begin, values_begin, X, y, knn, knn_size, degree);

						//calculate residuals
						residuals = y - X * (X.transpose() * X).inverse() * X.transpose() * y;
						stddev = sqrt((residuals.array() - residuals.mean()).square().sum() / (residuals.size() - 1));
						error = *(values_begin + i) - (Xs * (X.transpose() * X).inverse() * X.transpose() * y)(0);

						//calculate Residual Standard Error
						(*res)[group_start + i] = std::abs(error / stddev);
					}
				}
				break;
				case 1: //yhat
				{
					for (auto i = 0u; i < group_size; ++i) {
						//set self matrix
						for (auto j =0u;j<degree;++j){
							Xs(0,j+1) = std::pow(*(times_begin+i), j+1.0);
						}
						//parameters : don't sort, don't include ties (which means it will be sorted anyway), and don't include ourselves (at the 0 index!)
						knn = std::move(distmatrix.knn_distances(i, k, 0, false, 0));
						//set matrix and vector elements to our knn
                        knn_to_matrix(times_begin, values_begin, X, y, knn, knn_size, degree);

						//get yhat from the OLS solution
						(*res)[group_start + i] = (Xs * (X.transpose() * X).inverse() * X.transpose() * y)(0);
					}
				}
				break;
			}
		}
	}
};


//[[Rcpp::export]]
Rcpp::NumericVector tp_local_regression_outlier(Rcpp::IntegerVector groups,
                                Rcpp::NumericVector times,
                                Rcpp::NumericVector values,
								Rcpp::Nullable<Rcpp::NumericVector> Coeffs = R_NilValue,
                                std::string type = "timewindow",
								std::string output = "z_score",
                                int degree = 1,
                                int k = 5,
								int before = 900,
								int after = 900,
                                int nThreads = 0,
                                int threadThreshold = 1,
                                int threadMultiplier = 4,
                                int minThreads = 1,
								int dist_type = 0,
                                int which_beta = 0,
                                bool verbose = false)
{
	uint32_t n = groups.size();
	//the output vector
	Rcpp::NumericVector out(n);
	//get the pointers to the vectors
	Rcpp::IntegerVector *g = &groups;
	Rcpp::NumericVector *t = &times;
	Rcpp::NumericVector *v = &values;
	Rcpp::NumericVector *r = &out;
    
    uint32_t nthread = nThreads;


    //some dynamic thread pool calculations
  if (nThreads <= 0){
    nThreads = std::max(1u,std::thread::hardware_concurrency());
    if (verbose) Rcpp::Rcout << "nThreads: " << nThreads << "\n";
  }

	int type_int = 0;
	int output_int = 0;
    int output_beta = which_beta;
	
	std::vector<double> coeffs;

	//set the type correctly and throw error if illegal type
	if (type == "knn"){
		//set the distance coeffs for the knn-method
		if (Coeffs.isNotNull()){
			Rcpp::NumericVector c = Rcpp::NumericVector::create();
			c = Coeffs;
			coeffs.resize(c.size());
			std::copy(c.begin(), c.end(), coeffs.begin());
		}
		while(coeffs.size()<2){
			coeffs.push_back(1.0);
		}
		
		coeffs.resize(2);
		type_int = 0;
	}  else if (type == "knn_nonself"){
		//set the distance coeffs for the knn-method
		if (Coeffs.isNotNull()){
			Rcpp::NumericVector c = Rcpp::NumericVector::create();
			c = Coeffs;
			coeffs.resize(c.size());
			std::copy(c.begin(), c.end(), coeffs.begin());
		}
		while(coeffs.size()<2){
			coeffs.push_back(1.0);
		}
		
		coeffs.resize(2);
		type_int = 4;
	} else if (type == "timewindow"){
		type_int = 1;
	} else if (type == "window"){
		type_int = 2;
	} else {
		Rcpp::stop("'type' must be one of 'knn', 'timewindow', or 'window'");
	}

	//set output parameter and throw error if illegal output type
	if (output == "z_score"){
		output_int = 0;
	} else if (output == "y_hat"){
		output_int = 1;
	} else if (output == "r2"){
		output_int = 2;
	} else if (output == "b0"){
		output_int = 3;
        output_beta = 0;       
	} else if (output == "b1"){
		output_int = 3;
        output_beta = 1;
	} else if (output == "beta"){
		output_int = 3;
	} else {
		Rcpp::stop("'output' must be one of 'z_score', 'y_hat', 'r2', 'b0', 'b1' or 'beta'");
	}


	//set up the grouping parameters
	auto this_group = groups[0];
	uint32_t this_group_start = 0;
	uint32_t this_group_size = 0;
	//the vector of futures
	std::vector<std::future<void>> resvec;

	//the scope of the thread pool
	{
		threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
        if (verbose){
            Rcpp::Rcout << "pool created\n";
        }

		for (uint32_t i = 0;i<n;++i){
			if ((groups[i] != this_group) | (i == n - 1)) {

				//set the correct size of the currently found group
				if (i == n - 1) {
					this_group_size = i - this_group_start + 1;
				} else {
					this_group_size = i - this_group_start;
				}

				//resvec.push_back(std::move(pool.add_task()));
                resvec.push_back(std::move(pool.add_task(llr<Rcpp::NumericVector>(&times, &values, &out, coeffs, this_group_start, this_group_size, type_int, output_int, output_beta, degree, k, before, after))));

				//store this as the last found group
				this_group = groups[i];
				this_group_start = i;

			}
		}
        if (verbose){
            Rcpp::Rcout << "TP finished, nthreads = " << pool.get_size() << "\n";
        }
	}
	
	//get the futures
	for (auto &f : resvec){
		f.get();
	}

	return out;
  
}