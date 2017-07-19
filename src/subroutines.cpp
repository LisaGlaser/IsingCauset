#include "subroutines.h"

//Generate random total order for vector U
void randomTotalOrder(std::vector<unsigned int> &U, const int N)
{
	#if UNIT_TEST
	assert (N > 0);
	#endif

	U.resize(N);

	std::iota(U.begin(), U.end(), 0);
	std::random_shuffle(U.begin(), U.end());
	/*printf("random total order: (");
	for (int i = 0; i < N; i++)
		printf(" %d,", U[i]);
	printf(" )\n");*/
}

//Measure Causal Set Action
//Algorithm has been optimized using minimal bitwise operations
//Requires the existence of the whole adjacency matrix
//This will calculate all cardinality intervals by construction
bool measureAction_v3(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon)
{
	#if DEBUG
	assert (cardinalities != NULL);
	assert (adj.size() > 0);
	assert (stdim >= 2 && stdim <= 4);
	assert (N > 0);
	#endif

	int n = N + N % 2;
	uint64_t npairs = static_cast<uint64_t>(n) * (n - 1) / 2;
	uint64_t clone_length = adj[0].getNumBlocks();

	memset(cardinalities, 0, sizeof(uint64_t) * N * omp_get_max_threads());
	action = 0.0;

	//The first element will be N
	cardinalities[0] = N;

	unsigned int nthreads = omp_get_max_threads();
	#ifdef AVX2_ENABLED
	nthreads >>= 1;
	#endif

	//Compare all pairs of elements
	#ifdef _OPENMP
	#pragma omp parallel for schedule (dynamic, 64) if (npairs > 10000) num_threads (nthreads)
	#endif
	for (uint64_t k = 0; k < npairs; k++) {
		unsigned int tid = omp_get_thread_num();
		//Choose a pair
		uint64_t i = k / (n - 1);
		uint64_t j = k % (n - 1) + 1;
		uint64_t do_map = i >= j ? 1ULL : 0ULL;
		i += do_map * ((((n >> 1) - i) << 1) - 1);
		j += do_map * (((n >> 1) - j) << 1);

		//Ignore pairs which are not connected
		if (static_cast<int>(j) == N) continue;	//Arises due to index padding
		if (!adj[i].read(j)) continue;	//If elements i and j are not related, continue

		//Save the cardinality
		uint64_t length = j - i + 1;
		adj[i].clone(workspace[tid], 0ULL, clone_length);
		workspace[tid].partial_intersection(adj[j], i, length);
		cardinalities[tid*N+workspace[tid].partial_count(i, length)+1]++;
	}

	//Reduction for OpenMP
	for (int i = 1; i < omp_get_max_threads(); i++)
		for (int j = 0; j < N; j++)
			cardinalities[j] += cardinalities[i*N+j];

	action = calcAction(cardinalities, stdim, epsilon);
	if (action != action) return false;

	return true;
}

bool measureAction_v2(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon)
{
	#if DEBUG
	assert (cardinalities != NULL);
	assert (adj.size() > 0);
	assert (stdim >= 2 && stdim <= 4);
	assert (N > 0);
	#endif

	memset(cardinalities, 0, sizeof(uint64_t) * N * omp_get_max_threads());
	action = 0.0;

	//The first element will be N
	cardinalities[0] = N;

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (!adj[i].read(j)) continue;
			int num = 0;
			for (int k = i + 1; k < j; k++)
				if (adj[i].read(k) && adj[j].read(k)) num++;
			cardinalities[num+1]++;
		}
	}

	return true;
}

bool measureAction_v1(uint64_t * const cardinalities, double &action, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const unsigned int stdim, const int N, const double epsilon)
{
	#if DEBUG
	assert (cardinalities != NULL);
	assert (stdim >= 2 && stdim <= 4);
	assert (N > 0);
	#endif


	memset(cardinalities, 0, sizeof(uint64_t) * N * omp_get_max_threads());
	action = 0.0;

	//The first element will be N
	cardinalities[0] = N;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			int num = 0;
			if (U[i] < U[j] && V[i] < V[j]) {
				for (int k = 0; k < N; k++)
					if (U[i] < U[k] && U[k] < U[j] &&
					    V[i] < V[k] && V[k] < V[j])
						num++;
				cardinalities[num+1]++;
			}
		}
	}

	return true;
}

void printmatrix(Bitvector &mat, const int N)
{

	for(int i=0; i<N;i++)
	{
			mat[i].printBitset();
	}
}

void updateRelations(Bitvector &new_adj, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const int N)
{
	#if UNIT_TEST
	assert (new_adj.size() > 0);
	assert (U.size() > 0);
	assert (V.size() > 0);
	assert (N > 0);
	#endif

	/*printf("\nU:");
	for (int i = 0; i < N; i++)
		printf(" %d", U[i]);
	printf("\n");
	printf("V:");
	for (int i = 0; i < N; i++)
		printf(" %d", V[i]);
	printf("\n");*/

	for (int i = 0; i < N; i++)
		new_adj[i].reset();

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if ((U[i] < U[j] && V[i] < V[j]) || (U[j] < U[i] && V[j] < V[i])) {
				//printf("Adding link [%d-%d]\n", i, j);
				//new_adj[i].set(j);
				//new_adj[j].set(i);
				//printf("Adding link [%d-%d]\n", U[i], U[j]);
				new_adj[U[i]].set(U[j]);
				new_adj[U[j]].set(U[i]);
			}
			//} else
			//	printf("NOT adding link [%d-%d]\n", i,j);
		}
	}
}

void updateLinks(Bitvector &new_adj,Bitvector &new_link,  const int N)
{
	#if UNIT_TEST
	assert (new_adj.size() > 0);
	assert (new_link.size() > 0);
	assert (N > 0);
	#endif

	for (int i = 0; i < N; i++)
		new_link[i].reset();

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{

      if (!new_adj[i].read(j)) continue;

        if (!new_adj[i].partial_vecprod(new_adj[j], i, j-i+1))
				{

            new_link[i].set(j);

            new_link[j].set(i);
				}
		}
	}


}
