//2022-12-29
//A file that implements multithreaded versions of:
//LOF
//COF
//kNN-avg
//LOOP
//Throughout the file, templated on either 32 or 64 bit unsigned integers for indexing the arrays.
//Now with a revised version of COF that is faster for larger Ks, reducing the amount of redundant distance calculations.


#include <Rcpp.h>
//[[Rcpp::plugins("cpp14")]]

#include <algorithm>
#include <iostream>

#include <vector>
#include <queue>
#include <array>

#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <future>

#include <string>
#include <chrono>
#include <cmath> //for erf


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


//2022-12-20
//a struct-of-arrays-style object to return from knn-calculations, including a continous vector of all calculated pairwise distances calculated for this index.
template<typename Tint>
struct knn_row {
	//how far from the data start does this distrow begin
	Tint offset;
	//used for fast 'between'-test
	Tint diff;
	//a sorted vector of the k nearest neighbors to this point
	std::vector<Tint> knn;
	//the row of distances
	std::vector<double> dist;
};

//2022-12-21
//A class that generates k-nearest neighbors sequentially by sliding a submatrix through the diagonal of the full distance matrix and calculating
//the pairwise distances within the submatrix.
//3 return types are possible:
//1. A pair of K-distance, and a vector of all points <= the k-distance
//2. A (optionally sorted) vector of pairs of nearest neighbors and distance to that neighbor
//3. A struct containing the k-nearest neighbors, and the row from the distance matrix of the point in question as a continous vector, with an offset from the start iterator
template<typename Tit, typename Tint>
class dm_vec {
private:
	//std::vector<double> coeffs;
	//2022-12-18
	//keep the coefficients on the stack, should be faster - and who could ever want more than 8 dimensions?
	//std::array<double,8> coeffs;
	std::vector<double> coeffs;
	std::vector<Tit> starts;
	Tint size, dim;
	Tint capacity;
	Tint capacity_bits;
	Tint k;
	std::vector<std::vector<double>> mat;

	//Tint integers
	Tint window_begin = 0;
	Tint window_end = 0;

	//some vectors to avoid having to allocate memory all the time
	//the vector for calculating k-th distance
	std::vector<double> k_func_distances;
	//vectors for brute forcing back and forth
	std::vector<double> backfill;
	std::vector<double> frontfill;

	//bit shift hacking!
	//this calculates the i % capacity, if capacity is a power of two
	//inline Tint im(Tint i) {
	//	return i & (capacity - 1);
	//}

	//2022-12-18 - the capacity-1 expression and see if it helps speed a bit
	inline Tint im(const Tint i) {
		return i & capacity_bits;
	}

	//bit shift hacking 2, find the smallest p larger or equal to n where p is a power of two
	Tint nextPowerOf2(const Tint n)
	{
		Tint p = 1;
		if (n && !(n & (n - 1)))
			return n;

		while (p < n)
			p <<= 1;

		return p;
	}

	//Weighted euclidian distance
	//this actually calculates squared distance, to avoid having to call sqrt all the time.
	inline double getdist(const Tint p1, const Tint p2) {
		double res = 0.0;
		for (Tint i = 0u; i < dim; ++i) {
			res += std::pow((*(starts[i] + p1)) - (*(starts[i] + p2)), 2.0) * coeffs[i];
		}
		//res = std::sqrt(res);
		return res;
	}

	//weighted 1d distance
	//dito - calculate squared distance in 1 dimension
	inline double get_1d_dist(const Tint p1, const Tint p2, const Tint idx) {
		double res = std::pow((*(starts[idx] + p1)) - (*(starts[idx] + p2)), 2.0) * coeffs[idx];
		//res = std::sqrt(res);
		return res;
	}

	//get the size of our current window
	inline Tint matrix_size() {
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
		Tint ir = im(window_end);
		//ii is the column
		Tint ii = 0;
		for (Tint i = window_begin; i < window_end; ++i) {
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
	//increase the size of the matrix in the backward direction, i.e. adds the element at window_begin -1
	void distmatrix_push_front() {
		//remove the last element of the matrix if we're at full capacity
		while (matrix_size() >= capacity) {
			distmatrix_pop_back();
		}

		//update the matrix
		//ir is the current row in the matrix...
		Tint ir = im(window_begin - 1);
		//ii is the column
		Tint ii = 0;
		for (Tint i = window_begin; i < window_end; ++i) {
			ii = im(i);
			//calculate the distance between the new point (window_begin -1) and the current old point (i)
			//and add it to the row for our new point (iw = im(window_begin -1))
			mat[ir][ii] = getdist(window_begin - 1u, i);
			//make the matrix square!
			mat[ii][ir] = mat[ir][ii];
		}
		//add our self distance, i.e. 0
		mat[ir][ir] = 0;
		//update window_begin
		window_begin--;
	}

	//copy a row from the matrix into a vector, in the correct order
	template<typename Ta>
	inline void matrix_to_container(Ta& target, const Tint index) {
		Tint mat_size = matrix_size();
		target.resize(mat_size);
		Tint ir = im(index);
		for (Tint i = 0u; i < mat_size; ++i) {
			target[i] = mat[ir][im(i + window_begin)];
		}
	}

	//better resize target first!!
	template<typename Ta>
	inline void matrix_to_container2(Ta& target, const Tint index) {
		//get the matrix size
		Tint mat_size = matrix_size();
		//copy our row
		Tint ir = im(index);
		for (Tint i = 0u; i < mat_size; ++i) {
			target[i] = mat[ir][im(i + window_begin)];
		}
	}

	//get the kth element from a row in the matrix
	template<typename Ta>
	inline double get_k_dist_matrix(Ta& target, const Tint k, const Tint index) {
		Tint mat_size = matrix_size();
		target.resize(mat_size);
		//copy the matrix row to start of container
		matrix_to_container2(target, index);
		std::nth_element(target.begin(), target.begin() + k, target.end());
		return target[k];
	}

	template<typename Ta>
	inline double update_k_dist_matrix(Ta& target, const Tint k, const double newvalue) {
		target.push_back(newvalue);
		std::nth_element(target.begin(), target.begin() + k, target.end());
		return target[k];
	}

	//get the kth element from a row in the matrix and a container
	template<typename Ta, typename Tb>
	inline double get_k_dist_matrix(Ta& target, const Tint k, const Tint index, const Tb& b) {
		Tint mat_size = matrix_size();
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
	inline double get_k_dist_matrix(Ta& target, const Tint k, const Tint index, const Tb& b, const Tc& c) {
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
	inline void shrink_matrix_to_fit(Tint const index, const Tint k, const double k_dist) {
		while ((matrix_size() > k + 1u) && window_begin < index && (get_1d_dist(window_begin, index, 0) > k_dist)) {
			distmatrix_pop_front();
		}
	}

	//breaking away the part where the matrix is updated and the k_distance is calculated
	//note, it will change the backfill and frontfill vectors
	inline double get_k_dist(const Tint index, const Tint k) {
		//enlarge the matrix forwards if it's too small
		while ((matrix_size() < k + 1u) && (window_end < size)) {
			distmatrix_push_back();
		}
		//enlarge it backwards if it's still too small
		while ((matrix_size() < k + 1u) && (window_begin > 0u)) {
			distmatrix_push_front();
		}

		//if we've somehow managed to not insert our own row into the matrix, force it forwards until it's there
		while (window_end <= index) {
			distmatrix_push_back();
		}
		//dito, in the reverse direction
		while (window_begin > index && window_begin > 0u) {
			distmatrix_push_front();
		}



		//a guess of k_distance based on the current distance matrix
		double k_dist = get_k_dist_matrix(k_func_distances, k, index);
		//our row in the matrix, which will be useful in the future
		Tint ir = im(index);
		//this might be useful to preallocate as it will be used a lot
		double newdist;

		//try to enlarge the matrix forward, if we're not at the end of our view of the vector
		while (window_end < size) {
			//if the v0-distance is <= k-dist, and not at capacity -> enlarge forward
			if ((get_1d_dist(index, window_end, 0) <= k_dist) && (matrix_size() < capacity)) {
				distmatrix_push_back();
				//recalculate the k distance 
				k_dist = update_k_dist_matrix(k_func_distances, k, mat[ir][im(window_end - 1u)]);
			}
			else break;
		}

		//do a brute force seach backwards if there might be points earlier that are closer than our current kth-distance
		//this is instead of enlarging the matrix backwards just for a single outlier...
		backfill.clear();
		//now with Tint integers!
		Tint backfill_offset = 1;
		while (window_begin >= backfill_offset) {
			if (get_1d_dist(index, window_begin - backfill_offset, 0u) <= k_dist) {
				newdist = getdist(index, window_begin - backfill_offset);
				backfill.push_back(newdist);
				k_dist = update_k_dist_matrix(k_func_distances, k, newdist);
				backfill_offset++;
			}
			else break;
		}

		//dito in the forward direction if the matrix is at capacity.
		frontfill.clear();
		Tint frontfill_offset = 0;
		while (window_end + frontfill_offset < size) {
			if (get_1d_dist(index, window_end + frontfill_offset, 0u) <= k_dist) {
				newdist = getdist(index, window_end + frontfill_offset);
				frontfill.push_back(newdist);
				k_dist = update_k_dist_matrix(k_func_distances, k, newdist);
				frontfill_offset++;

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
		for (Tint i = window_begin; i < window_end; ++i) {
			Tint ir = im(i);
			for (Tint j = window_begin; j < window_end; ++j) {
				std::cout << mat[ir][im(j)] << " ";
			}
			std::cout << "\n";
		}
	}

public:
	dm_vec(const std::vector<Tit>& starts, const std::vector<double>& coeffs, const Tint size, const Tint k, const Tint matrix_buffer) {
		this->starts = starts;
		dim = starts.size();

		this->coeffs = coeffs;
		//copy the coeffs
		// for (auto i = 0u;i<dim;++i){
		//     this->coeffs[i] = coeffs[i];
		// }

		//copy the constants
		this->size = size;
		this->k = k;

		//very computer science-y!
		//+4 is a empiric value found from benchmarking - it gets expensive when the capacity is too close to k, it would seem...
		//this->capacity = nextPowerOf2(k + 4u);
		capacity = nextPowerOf2(k + std::max(1, static_cast<int>(matrix_buffer)));
		//precalculate this to keep it on the stack in all future im-calls
		capacity_bits = capacity - 1u;
		window_begin = 0;
		window_end = 0;

		//reserve some space for the distance calculations
		k_func_distances.reserve(2 * capacity);

		//populate the matrix
		mat.reserve(capacity);
		for (Tint i = 0u; i < capacity; ++i) {
			mat.emplace_back(std::vector<double>(capacity));
			//std::vector<double> vec(this->capacity);
			//mat.push_back(vec);
		}
	}

	//knn_row<Tint> get_knn_with_row(Tint index, Tint k){
	//this returns a knn_row object, i.e. a sorted vector of k-nearest neighbors (by distance), a vector of doubles corresponding to the row in
	//the distance matrix for the index, and an offset and the size
	knn_row<Tint> get_knn_with_row(const Tint index, const Tint k) {
		//this won't be used but the function call is used for its side effects (i.e. updating matrix and backfill/frontfill vectors)
		double k_dist{ get_k_dist(index, k) };

		//precalculate this as it will be used a couple of times
		Tint knn_row_size{ matrix_size() + backfill.size() + frontfill.size() };

		//initialize our result object
		knn_row<Tint> res{ window_begin - backfill.size(),
			knn_row_size - 1u,
			std::vector<Tint>(0),
			std::vector<double>(0) };
		//reserve space in the vectors
		res.knn.reserve(knn_row_size - 1u);
		res.dist.reserve(knn_row_size);

		//store this to keep calling back/frontfill.size() too much
		Tint vsize{ backfill.size() };

		//first add the backfill distances in reverse order
		//another old-school c++ syntax hack :)
		for (Tint i = vsize; i--;) {
			//res.knn.emplace_back(window_begin - vsize + i);
			res.dist.emplace_back(std::sqrt(backfill[i]));
		}
		//second, add all points from the matrix except for self
		Tint ir{ im(index) };

		for (Tint i = window_begin; i < window_end; ++i) {
			//don't add ourselves to our knn as it will be a hassle to remove us later on...
			if (i != index) {
				res.knn.emplace_back(i);
			}
			res.dist.emplace_back(std::sqrt(mat[ir][im(i)]));
		}

		//third, add the points from the backfill in non-reverse order...
		for (Tint i = 0u; i < vsize; ++i) {
			res.knn.emplace_back(window_begin - 1u - i);
		}

		//add the frontfill points
		vsize = frontfill.size();
		for (Tint i = 0u; i < vsize; ++i) {
			res.knn.emplace_back(window_end + i);
			res.dist.emplace_back(std::sqrt(frontfill[i]));
		}
		//sort the index by distance
		//this nth-element reduces the work done when the size of the matrix row >> k
		std::nth_element(res.knn.begin(), res.knn.begin() + k, res.knn.end(), [&](const auto& a, const auto& b) {return res.dist[a - res.offset] < res.dist[b - res.offset]; });
		std::sort(res.knn.begin(), res.knn.begin() + k, [&](const auto& a, const auto& b) {return res.dist[a - res.offset] < res.dist[b - res.offset]; });
		//std::sort(res.knn.begin(), res.knn.end(), [&](const auto& a, const auto& b) {return res.dist[a - res.offset] < res.dist[b - res.offset]; });
		//shrink knn to the correct number...
		res.knn.resize(k);
		//shrink our matrix
		shrink_matrix_to_fit(index, k, k_dist);

		return res;
	}

	//get the k nearest neighbors, with optional sorting of results on either distance or index, and may or may not include ties.
	std::vector<std::pair<Tint, double>> knn_distances(const Tint index, const Tint k, int sort_results_by = 0, bool include_ties = false, int include_self = 0) {
		//get our k-distance, and our neighborhood will be stored in our corresponding matrix row, plus maybe some points in the 
		//backfill or frontfill vectors - if they are outside the window shown by the matrix.
		double k_dist = get_k_dist(index, k);
		std::vector<std::pair<Tint, double>> res;
		res.reserve(k + 4);
		Tint ir = im(index);

		//include ourselves if that is requested
		if (include_self == 1) {
			res.emplace_back(index, 0.0);
		}

		//this will keep track of if we're already sorted the out-vector when we removed any ties...
		int sorted_by = 0;
		//only for the first group, to keep this sane
		//add the points from the matrix that are part of our neighborhood
		for (Tint i = 0u; i < matrix_size(); ++i) {
			Tint ii = im(i + window_begin);
			if (mat[ir][ii] <= k_dist && (window_begin + i != index)) {
				//res.emplace_back(std::make_pair<Tint, double>(window_begin + i, std::sqrt(mat[ir][ii])));
				res.emplace_back(window_begin + i, std::sqrt(mat[ir][ii]));
			}
		}


		//and the point from our backfill that are part of our knn
		for (Tint i = 0u; i < backfill.size(); ++i) {
			if (backfill[i] <= k_dist) {
				//res.emplace_back(std::make_pair<Tint, double>(window_begin - 1 - i, std::sqrt(backfill[i])));
				res.emplace_back(window_begin - 1 - i, std::sqrt(backfill[i]));
			}
		}

		//and the points from the frontfill
		for (Tint i = 0u; i < frontfill.size(); ++i) {
			if (frontfill[i] <= k_dist) {
				//res.emplace_back(std::make_pair<Tint, double>(window_end + i, std::sqrt(frontfill[i])));
				res.emplace_back(window_end + i, std::sqrt(frontfill[i]));
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
	std::pair<double, std::vector<Tint>> k_dist_neighborhood(const Tint index, const Tint k) {
		//get our k-distance, and our neighborhood will be stored in our corresponding matrix row, plus maybe some points in the 
		//backfill or frontfill vectors - if they are outside the window shown by the matrix.
		double k_dist = get_k_dist(index, k);
		//the vector of the indexes of our closest neighborhood
		std::vector<Tint> locations;
		//reserve some extra memory for this vector
		locations.reserve(k + 2);
		//precalculate im(index), i.e. our current row
		Tint ir = im(index);

		//add the points in the matrix that are part of our k-distance-neighborhood
		for (Tint i = 0u; i < matrix_size(); ++i) {
			Tint ii = im(i + window_begin);
			//<= to allow k-neighborhood larger than k if there are ties
			if (mat[ir][ii] <= k_dist && (window_begin + i) != index) {
				locations.push_back(window_begin + i);
			}
		}

		//add the points from the backfill that are part of the neigborhood
		for (Tint i = 0u; i < backfill.size(); ++i) {
			//<= to allow k-neighborhood larger than k if there are ties
			if (backfill[i] <= k_dist) {
				locations.push_back(window_begin - 1 - i);
			}
		}

		//and the parts of the frontfill
		for (Tint i = 0u; i < frontfill.size(); ++i) {
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

//2022-12-21
//A rewrite of the COF-class, this time with the aim to reduce the amount of distance recalculations by returning a custom sparse distance matrix from the
//knn-query
template <typename T, typename Tit, typename Tint, typename Tres>
class cof_new {
private:
	std::vector<Tit> starts;
	//why use std::vector when an array could possibly be a bit faster here...
	//std::array<double,8> coeffs;
	std::vector<double> coeffs;
	Tres* out;
	Tint start, size, dim, k, matrix_buffer;
	//Tint foundfirst = 0;
	//Tint foundsecond = 0;
	//Tint recounts = 0;

	//since we're working with actual distances, not squared distances...
	inline double getdist(const Tint p1, const Tint p2) {
		double res = 0.0;
		for (auto i = 0u; i < dim; ++i) {
			res += std::pow((*(starts[i] + p1)) - (*(starts[i] + p2)), 2.0) * coeffs[i];
		}
		res = std::sqrt(res);
		return res;
	}

	//get the lowest value from the lower triangle of a row in the dist matrix, i.e. the shortest linking distance for that element.
	template <typename Ta, typename Tb>
	double min_linking_distance(const Tint index, const Tint row_number, const Ta& knnres, Tb& dist_row) {
		//since c++ uses 0 based indexes, the row numbers will go from 0... k-1 whereas the number of points we're looking at is 1...k
		//dist_row.resize(row_number + 1u);
		//the id of point one, i.e. the row_numberth nearest neighbor to the point at index
		Tint id = knnres[index].knn[row_number];
		//the first column is equal to the first row in the full matrix == the sorted vector of knns for the index, so the distance is guaranteed to be known
		dist_row[0] = knnres[index].dist[id - knnres[index].offset];

		//the loop to populate the rest of the row
		Tint targetid = 0u;
		for (Tint i = 0u; i < row_number; ++i) {
			//the id of point2, i.e. the ith nearest neighbor to the point at index
			targetid = knnres[index].knn[i];

			//first check if point2 is within point1s known distances, if so just fetch it
			//then, check if point1 is within point2s known distances, if so just fetch it
			//if it's a new pair, the distance has to be recalculated, of course.
			//if (targetid >= knnres[id].offset && targetid < knnres[id].end_of_row) {
			
			//this is a quick way of checking if an unsigned integer lies between two values
			if (targetid - knnres[id].offset <= knnres[id].diff) {
				dist_row[i + 1] = knnres[id].dist[targetid - knnres[id].offset];
			}
			//it would seem from intuition that it is faster to see if we've precalculated these distances before
			//but from benchmarking, it is marginally faster to just recalculate all distances not known from knnres[id].dist.
			//probably due to cache-friendlyness, I guess...
			// else if (id - knnres[targetid].offset <= knnres[targetid].diff) {
			// 	dist_row[i + 1] = knnres[targetid].dist[id - knnres[targetid].offset];
			// }
			else {
				dist_row[i + 1] = getdist(id, targetid);
			}
		}
		//the lowest element from relevant part of the vector, i.e. up to row_number + 1
		return *(std::min_element(dist_row.begin(), dist_row.begin() + row_number + 1u));
	}

public:
	cof_new(std::vector<T*>& Data, Tres* out, std::vector<double>& coeffs, const Tint start, const Tint size, const Tint k = 3, const Tint matrix_buffer = 16) {
		this->k = k;
		this->size = size;
		this->start = start;
		this->matrix_buffer = matrix_buffer;
		//copy the coeffs
		dim = Data.size();

		this->coeffs = coeffs;

		this->out = out;

		for (auto& v : Data) {
			starts.push_back(v->begin() + start);
		}

		//make sure k isn't larger than size - 1, otherwise it's impossible to find k neighbors
		if (k >= size) {
			this->k = std::max(static_cast<int>(size) - 1, 1);
		}
	}

	auto operator()()->void {
		//vector for the knn-neighbors
		std::vector<knn_row<Tint>> knnres;
		knnres.reserve(size);

		{
			dm_vec<Tit, Tint> distmatrix(starts, coeffs, size, k, matrix_buffer);
			//for now, just print the knn and see how that goes...
			for (Tint i = 0u; i < size; ++i) {
				knnres.emplace_back(distmatrix.get_knn_with_row(i, k));
			}
		}

		//reuse this vector to avoid having to reallocate memory all the time...
		std::vector<double> dist_row(k + 1u);

		//this will contain the sum of the weighted chaining distances, per object
		std::vector<double> sum_chaining_distances(size);
		

		//precalculate this entity which will be used in each loop
		double kk1 = k * (k + 1.0);
		//calculate the sum of weighted chaining distances per object
		for (Tint i = 0u; i < size; ++i) {
			//the linking distances function, i.e. the smallest
			sum_chaining_distances[i] = 0.0;
			for (Tint j = 0u; j < k; ++j) {
				//the weighted linking distance
				sum_chaining_distances[i] += ((2.0 * (k - j)) / kk1) * min_linking_distance(i, j, knnres, dist_row);
			}
		}

		//preallocate these entities which will be used many times.
		double avg_chain_dist_obs;
		double avg_chain_dist_nn;

		for (Tint i = 0u; i < size; ++i) {
			avg_chain_dist_obs = sum_chaining_distances[i] * k;
			avg_chain_dist_nn = 0.0;
			for (Tint j = 0u; j < k; ++j) {
				avg_chain_dist_nn += sum_chaining_distances[knnres[i].knn[j]];
			}
			//cof = avg_chain_dist_obs / avg_chain_dist_nn
			//std::cout << "i: " << i << " cof: " << avg_chain_dist_obs / avg_chain_dist_nn << "\n";

			(*out)[start + i] = avg_chain_dist_obs / avg_chain_dist_nn;
		}
	}
};



//2022-02-05
//COF implemented as c++ class, a trial...
template <typename T, typename Tit, typename Tint, typename Tres>
class cof {
private:
	std::vector<Tit> starts;
    //why use std::vector when an array could possibly be a bit faster here...
    //std::array<double,8> coeffs;
    std::vector<double> coeffs;
    Tres* out;
	Tint start, size, dim, k, matrix_buffer;

	//since we're working with actual distances, not squared distances...
	inline double getdist(Tint p1, Tint p2){
		double res = 0.0;
		for (auto i = 0u; i < dim; ++i) {
			res += std::pow((*(starts[i] + p1)) - (*(starts[i] + p2)), 2.0) * coeffs[i];
		}
		res = std::sqrt(res);
		return res;
	}

	//get the lowest value from the lower triangle of a row in the dist matrix, i.e. the shortest linking distance for that element.
	template <typename Ta, typename Tb>
	double min_linking_distance(Tint index, Tint row_number,  const Ta& knn,  Tb &dist_row) {
		//since c++ uses 0 based indexes, the row numbers will go from 0... k-1 whereas the number of points we're looking at is 1...k
		dist_row.resize(row_number + 1u);
		//the first column is equal to the first row in the full matrix == the sorted vector of knns for the index, so the distance is guaranteed to be known
		dist_row[0] = knn[index][row_number].second;
        
		//the loop to populate the rest of the row
		for (auto i = 0u; i < row_number; ++i) {
			//dist_row[i + 1] = find_dist(knn[index][row_number].first, knn[index][i].first, knn);

			//from benchmarking, it's actually faster to hard recalculate these distances than fetching them from the precalculated distances!
			dist_row[i + 1] = getdist(knn[index][row_number].first, knn[index][i].first);
		}
		return *(std::min_element(dist_row.begin(), dist_row.end()));
	}

public:
	cof(std::vector<T*>& Data, Tres *out, std::vector<double>& coeffs, Tint start, Tint size, Tint k = 3, Tint matrix_buffer = 16) {
		this->k = k;
		this->size = size;
		this->start = start;
		this->matrix_buffer = matrix_buffer;
        //copy the coeffs
        dim = Data.size();

        this->coeffs = coeffs;
        //for (auto i = 0u;i<dim;++i){
        //    coeffs[i] = coeffs_knn[i];
        //}

        this->out = out;

		for (auto& v : Data) {
			starts.push_back(v->begin() + start);
		}

		//make sure k isn't larger than size - 1, otherwise it's impossible to find k neighbors
		if (k >= size) {
			this->k = std::max(static_cast<int>(size) - 1, 1);
		}
	}

	auto operator()()->void {
		std::vector<std::vector<std::pair<Tint, double>>> knn;
		//std::vector<double> cof(size);
		//get the knn-neighbors
		{
			dm_vec<Tit, Tint> distmatrix(starts, coeffs, size, k, k+matrix_buffer);
			//for now, just print the knn and see how that goes...
			for (Tint i = 0u; i < size; ++i) {
				auto knn_result = distmatrix.knn_distances(i, k, 2);
				//res[i] = knn_result.size();
				knn.emplace_back(std::move(knn_result));
			}
		}
		
		//second, try to recreate the distance matrix...
		//reuse these vectors to avoid having to reallocate memory all the time...
		//std::vector<double> link_distances(k);
		//std::vector<double> avg_chaining_distances(k);		
		std::vector<double> dist_row(k);

		//this will contain the sum of the weighted chaining distances, per object
		std::vector<double> sum_chaining_distances(size);

        //calculate the sum of weighted chaining distances per object
		for (Tint i = 0u; i < size; ++i) {
			//the linking distances function, i.e. the smallest
            sum_chaining_distances[i] = 0.0;
			for (Tint j = 0u; j < k; ++j) {
			    //link_distances[j] = min_linking_distance(i, j, knn, dist_row);
                //avg_chaining_distances[j] = ((2.0 * (k - j)) / (k * (k + 1.0))) * link_distances[j];
                //sum_chaining_distances[i]+=avg_chaining_distances[j];
				sum_chaining_distances[i]+=((2.0 * (k - j)) / (k * (k + 1.0))) *min_linking_distance(i, j, knn, dist_row);
			}
		}

		

		double avg_chain_dist_obs;
		double avg_chain_dist_nn;

		for (Tint i = 0u; i < size; ++i) {
			avg_chain_dist_obs = sum_chaining_distances[i] * k;
			avg_chain_dist_nn = 0.0;
			for (Tint j = 0u; j < k; ++j) {
				avg_chain_dist_nn += sum_chaining_distances[knn[i][j].first];
			}
            //cof = avg_chain_dist_obs / avg_chain_dist_nn
			(*out)[start+i] = avg_chain_dist_obs / avg_chain_dist_nn;
		}
	}
};

//2022-02-01
//this version splits the distance matrix object into a separate class, to enable (in the future) multiple kinds of distance objects to be used
template <typename T, typename Tit, typename Tint, typename Tres>
class lof {
private:
	std::vector<Tit> starts;
    //std::array<double,8> coeffs;
    std::vector<double> coeffs;
    Tres* out;
	Tint start, size, dim, k, matrix_buffer;
	
	//Weighted euclidian distance
	inline double getdist(Tint p1, Tint p2) {
		double res = 0.0;
		for (Tint i = 0u; i < dim; ++i) {
			res += std::pow((*(starts[i] + p1)) - (*(starts[i] + p2)), 2.0) * coeffs[i];
		}
		res = std::sqrt(res);
		return res;
	}

public:
	lof(std::vector<T*>& Data, Tres *out, std::vector<double>& coeffs, Tint start, Tint size, Tint k = 3, Tint matrix_buffer = 16) {
		this->k = k;
		this->size = size;
		this->start = start;
		this->matrix_buffer = matrix_buffer;
        this->out = out;
        this->coeffs = coeffs;

        //make sure we have the right amount of coefficients
		dim = Data.size();
        //for (auto i = 0u;i<dim;++i){
        //    coeffs[i] = coeffs_knn[i];
        //}

		for (auto& v : Data) {
			starts.push_back(v->begin() + start);
		}

		//make sure k isn't larger than size - 1, otherwise it's impossible to find k neighbors
		if (k >= size) {
			this->k = std::max(static_cast<int>(size) - 1, 1);
		}
	}

	auto operator ()()->void {
		//reserving memory for our vectors
		std::vector<std::pair<double, std::vector<Tint>>> knn;
		std::vector<double> lrd;

		knn.reserve(size);
		lrd.resize(size);


        {
            dm_vec<Tit, Tint> distmatrix(starts, coeffs, size, k, matrix_buffer);
			//get the knn-neighborhood at each point
			for (auto i = 0u; i < size; ++i) {
				knn.emplace_back(distmatrix.k_dist_neighborhood(i, k));
			}
        }
        //return;

		//step 2, get the LRD at each location
		for (auto i = 0u; i < size; ++i) {
			double reachability_distance_sum = 0.0;
			auto knn_size = knn[i].second.size();
			for (auto j : knn[i].second) {
				reachability_distance_sum += std::max(knn[j].first, getdist(i, j));
			}
			lrd[i] = knn_size / reachability_distance_sum;
		}

		//step 3, calculate LOF
		for (auto i = 0u; i < size; ++i) {
			double lrd_sum = 0.0;
			auto knn_size = knn[i].second.size();
			for (auto j : knn[i].second) {
				lrd_sum += lrd[j];
			}
			(*out)[start+i] = lrd_sum / (lrd[i] * knn_size);
		}
	}
};


//2022-05-15
//kNN - average distance to k nearest neighbors
template <typename T, typename Tit, typename Tint, typename Tres>
class knn_avg {
private:
	std::vector<Tit> starts;
	//std::vector<double> coeffs_knn;
    std::vector<double> coeffs;
    //std::array<double,8> coeffs;
    Tres* out;
	Tint start, size, dim, k, matrix_buffer;
	
	//Weighted euclidian distance
	inline double getdist(Tint p1, Tint p2) {
		double res = 0.0;
		for (Tint i = 0u; i < dim; ++i) {
			res += std::pow((*(starts[i] + p1)) - (*(starts[i] + p2)), 2.0) * coeffs[i];
		}
		res = std::sqrt(res);
		return res;
	}

public:
	knn_avg(std::vector<T*>& Data, Tres *out, std::vector<double>& coeffs, Tint start, Tint size, Tint k = 3, Tint matrix_buffer = 16) {
		this->k = k;
		this->size = size;
		this->start = start;
        this->out = out;
        this->matrix_buffer = matrix_buffer;
        this->coeffs = coeffs;
        dim = Data.size();
        //for (auto i = 0u;i<dim;++i){
        //    coeffs[i] = coeffs_knn[i];
        //}

		for (auto& v : Data) {
			starts.push_back(v->begin() + start);
		}

		//make sure k isn't larger than size - 1, otherwise it's impossible to find k neighbors
		if (k >= size) {
			this->k = std::max(static_cast<int>(size) - 1, 1);
		}
	}

	auto operator()()->void {
		//std::vector<double> cof(size);

		//std::vector<double> res(size);
		//get the knn-neighbors
		{
			dm_vec<Tit, Tint> distmatrix(starts, coeffs, size, k, matrix_buffer);
			//for now, just print the knn and see how that goes...
			Tint m;
			double sum;
			for (auto i = 0u; i < size; ++i) {
				auto knn_result = distmatrix.knn_distances(i, k, 2);
				m = knn_result.size();
				sum = 0.0;
                for (auto &p : knn_result){
                        sum+=p.second;
                }
				(*out)[start+i]  = sum / m;
			}
		}
	}
};


//2022-05-15
//local outlier probability - LOOP.
template <typename T, typename Tit, typename Tint, typename Tres>
class loop {
private:
	std::vector<Tit> starts;
    //std::array<double,8> coeffs;
    std::vector<double> coeffs;
    Tres* out;
	Tint start, size, dim, k, matrix_buffer;
    double lambda;
	
	//Weighted euclidian distance
	inline double getdist(Tint p1, Tint p2) {
		double res = 0.0;
		for (Tint i = 0u; i < dim; ++i) {
			res += std::pow((*(starts[i] + p1)) - (*(starts[i] + p2)), 2.0) * coeffs[i];
		}
		res = std::sqrt(res);
		return res;
	}

public:
	loop(std::vector<T*>& Data, Tres *out, std::vector<double>& coeffs, Tint start, Tint size, Tint k = 3, double lambda = 2.0, Tint matrix_buffer = 16) {
		this->k = k;
		this->size = size;
		this->start = start;
        this->out = out;
        this->matrix_buffer = matrix_buffer;

        dim = Data.size();
        this->coeffs = coeffs;
        //for (auto i = 0u;i<dim;++i){
        //    coeffs[i] = coeffs_knn[i];
        //}

        this->lambda = lambda;

		for (auto& v : Data) {
			starts.push_back(v->begin() + start);
		}

		//make sure k isn't larger than size - 1, otherwise it's impossible to find k neighbors
		if (k >= size) {
			this->k = std::max(static_cast<int>(size) - 1, 1);
		}
	}

	auto operator()()->void {
		std::vector<std::vector<std::pair<Tint, double>>> knn;
		//std::vector<double> loop(size);
		std::vector<double> distSD(size);


		//declare on stack for quick mean calculations
		double sum = 0.0;
		std::size_t m = 0;
		//get the knn-neighbors
		{
			dm_vec<Tit, Tint> distmatrix(starts, coeffs, size, k, matrix_buffer);
			//for now, just print the knn and see how that goes...
			for (auto i = 0u; i < size; ++i) {
				auto knn_result = distmatrix.knn_distances(i, k, 2);
				//res[i] = knn_result.size();
				//get the SD of the distances while we're at it
				m = knn_result.size();
				sum = 0.0;
                for (auto &p : knn_result){
					    sum+=p.second*p.second;
				}

				distSD[i] = std::sqrt(sum / m);

				//store the result
				knn.emplace_back(std::move(knn_result));
			}
		}

		//return distSD;
		//get the outlier-ishness of our point
		std::vector<double> plof(size);
		double nplof = 0.0;

		for (auto i = 0u;i<size;++i){
			m = knn[i].size();
			//calculate the mean of the sd of the distances of our neighbors
			sum = 0.0;
			for (auto &p : knn[i]){
				sum+=distSD[p.first];
			}
			plof[i] = (distSD[i] / (sum/m))-1.0;
			nplof += plof[i]*plof[i];
		}
		//return plof;
		nplof = lambda*std::sqrt(nplof/size) * std::sqrt(2.0);
		for (auto i = 0u;i<size;++i){
			(*out)[start+i] = std::max(std::erf(plof[i]/nplof), 0.0);
		}
	}
};

//2022-12-19
//A new version that is somewhat cleaned up compared to before, including less redundant code.
//[[Rcpp::export]]
Rcpp::NumericVector density_outliers_im(Rcpp::IntegerVector groups,
                                Rcpp::DataFrame df,
                                Rcpp::Nullable<Rcpp::NumericVector> coeffs_knn = R_NilValue,
                                std::string type = "knn_avg",
                                int k = 5,
                                int nThreads = 0,
                                int threadThreshold = 1,
                                int threadMultiplier = 4,
                                int minThreads = 1,
								double lambda = 2.0, 
								int matrix_buffer = 16,
                                bool verbose = false,
                                bool use_64bit_index = false)
{
    std::size_t n = groups.size();
    Rcpp::NumericVector out(n);
    std::vector<std::future<void>> resvec;
    //we're gonna pass data as a vector of iterators
    //this vector of vectors represents our data as numericvectors
    std::vector<Rcpp::NumericVector> data_num;
     //vector of pointers to the vectors
     std::vector<Rcpp::NumericVector*> data_ptr;
     std::vector<double> c_knn;
     //some dynamic thread pool calculations
     if (nThreads <= 0){
         nThreads = std::max(1u,std::thread::hardware_concurrency());
         if (verbose) Rcpp::Rcout << "nThreads: " << nThreads << "\n";
    }
    //this will be by reference for numeric, but copy if it's integer (since it's cast to numeric)
    for (auto &v:df){
        //14 numeric, 13 integer
        if (TYPEOF(v) == 14 || TYPEOF(v) == 13){
            data_num.push_back(v);
            if (verbose) Rcpp::Rcout << "pushing back vector at: " << &v << "\n";
        }
    }
    
    //make sure we only pass pointers to the function
    for (auto &v:data_num){
        data_ptr.push_back(&v);
        if (verbose) Rcpp::Rcout << "pointer addresses: " << &v << "\n";
    }
    //copy the coeffs
    if (coeffs_knn.isNotNull()){
        //fetch the vector
        Rcpp::NumericVector c = Rcpp::NumericVector::create();
        c = coeffs_knn;
        //copy into std::vector
        c_knn.resize(c.size());
        std::copy(c.begin(), c.end(), c_knn.begin());
    }

    //make sure we have coeffs for each variable and substitute 1.0 if missing
    while (c_knn.size() < data_ptr.size()){
        c_knn.push_back(1.0);
        if (verbose){
            Rcpp::Rcout << "creating new kNN-coefficient: 1.0\n";
        }
    }

    //also shrink the coeffient vectors if they are too large
    c_knn.resize(data_ptr.size());


	if (verbose){
        
		Rcpp::Rcout << "Coefficients (kNN): ";
        for (auto cc : c_knn){
            Rcpp::Rcout << cc <<" ";
        }
        Rcpp::Rcout << "\n";
    }


    
    if (!use_64bit_index){
        //n = number of rows in total
        //grouping variables
        int this_group = groups[0];
        uint32_t this_group_start = 0;
        uint32_t this_group_size = 0;
        //the scope of the pool...
        {
            //threadpool_d pool(nThreads);
            threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
            //threadpool_d pool(1,1,1,1);
            if (verbose){
                Rcpp::Rcout << "pool created\n";
            }
            
            //the grouping loop
            for (uint32_t i = 0u; i < n; ++i) {
                //set the size of the group
                if ((groups[i] != this_group) | (i == n - 1)) {
                    if (i == n - 1) {
                        this_group_size = i - this_group_start + 1;
                    }
                    else {
                        this_group_size = i - this_group_start;
                    }
                    
                    //copy the vector of iterators
                    if (verbose){
                        Rcpp::Rcout << "found a group at " << this_group_start << " of size " << this_group_size << "\n";
                    }

                    if (type == "lof" || type == "lof_vec"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task lof(uint32_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(lof<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint32_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));

                    } else if (type == "knn" || type == "knn_avg"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task knn_avg(uint32_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(knn_avg<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint32_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));
                        
                    } else if (type == "loop"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task loop(uint32_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(loop<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint32_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, lambda, matrix_buffer))));
                    } else if (type == "cof_old"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task cof(uint32_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(cof<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint32_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));
                    }else if (type == "cof"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task cof_new(uint32_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(cof_new<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint32_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));
                    }
					this_group = groups[i];
        			this_group_start = i;
                }
            }
        }
    } else {
        //n = number of rows in total
        //grouping variables
        int this_group = groups[0];
        uint64_t this_group_start = 0;
        uint64_t this_group_size = 0;
        //the scope of the pool...
        {
            //threadpool_d pool(nThreads);
            threadpool_d pool(nThreads, threadThreshold, threadMultiplier, minThreads);
            //threadpool_d pool(1,1,1,1);
            if (verbose){
                Rcpp::Rcout << "pool created\n";
            }
            
            //the grouping loop
            for (uint64_t i = 0u; i < n; ++i) {
                //set the size of the group
                if ((groups[i] != this_group) | (i == n - 1)) {
                    if (i == n - 1) {
                        this_group_size = i - this_group_start + 1;
                    }
                    else {
                        this_group_size = i - this_group_start;
                    }
                    
                    //copy the vector of iterators
                    if (verbose){
                        Rcpp::Rcout << "found a group at " << this_group_start << " of size " << this_group_size << "\n";
                    }

                    if (type == "lof" || type == "lof_vec"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task lof(uint64_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(lof<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint64_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));

                    } else if (type == "knn"|| type == "knn_avg"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task knn_avg(uint64_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(knn_avg<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint64_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));
                        
                    } else if (type == "loop"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task loop(uint64_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(loop<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint64_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, lambda, matrix_buffer))));
                    } else if (type == "cof_old"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task cof(uint64_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(cof<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint64_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));
                    }else if (type == "cof"){
                        if (verbose){
                            Rcpp::Rcout << "pushing back task cof_new(uint64_t)...\n";
                        }
                        resvec.push_back(std::move(pool.add_task(cof_new<Rcpp::NumericVector, Rcpp::NumericVector::iterator, uint64_t, Rcpp::NumericVector>(data_ptr, &out, c_knn, this_group_start, this_group_size, k, matrix_buffer))));
                    }
					this_group = groups[i];
        			this_group_start = i;
                }
            }
        }
    }

	//get the futures, i.e. wait for the threads to be finished
    for (auto &r : resvec){
        r.get();
    }

    return out;
}

