#include "subroutines.h"

// found this here: http://www.cplusplus.com/forum/general/33606/
// Our special 'read a line as a vector' object
template <typename T>
struct linevector_t: public std::vector <T> { };

// Cast function (mutable object)
template <typename T>
linevector_t <T> &
linevector( std::vector <T> & v )
  {
  return static_cast <linevector_t <T> &> ( v );
  }
// Cast functions (const object)

template <typename T>
const linevector_t <T> &

linevector( const std::vector <T> & v )
  {
  return static_cast <const linevector_t <T> &> ( v );
  }

// Input operation to read a single line into a vector <T>
template <typename T>
std::istream& operator >> ( std::istream& ins, linevector_t <T> & v )
  {
  std::string s;
  std::getline( ins, s );
  std::istringstream iss( s );
  v.clear();
  copy( std::istream_iterator <T> ( iss ), std::istream_iterator <T> (), std::back_inserter( v ) );
  return ins;
  }

// Something useful for our example program
template <typename T>
std::ostream& operator << ( std::ostream& outs, const linevector_t <T> & v )
  {
  copy( v.begin(), v.end() - 1, std::ostream_iterator <T> ( outs, " " ) );
  return outs << v.back() << std::endl;
  }



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

void initialState(Graph * const graph)
{

	if(!graph->props.initialstate.compare("cold"))
	{
		graph->props.U.resize(graph->props.N);
		std::iota(graph->props.U.begin(), graph->props.U.end(), 0);

		graph->props.V.resize(graph->props.N);
		std::iota(graph->props.V.begin(), graph->props.V.end(), 0);

		graph->spins=std::vector<int>(graph->props.N,1);
    printf("starting in a cold state \n");
	}
	else if(!graph->props.initialstate.compare("random"))
	{
    graph->spins.resize(graph->props.N);
    graph->props.U.resize(graph->props.N);
    graph->props.V.resize(graph->props.N);

		randomSpinState(graph->spins, graph->props.mrng, graph->props.N);
		randomTotalOrder(graph->props.U, graph->props.N);
		randomTotalOrder(graph->props.V, graph->props.N);
    printf("starting in a random state \n");

	}
	else
	{
		std::ifstream f(graph->props.initialstate);
		f >> linevector( graph->spins );
  	f >> linevector( graph->props.U );
  	f >> linevector( graph->props.V );
  	f.close();
    printf("starting from %s", graph->props.initialstate.c_str());

		/*printvector(graph->spins,graph->props.N);
		uprintvector(graph->props.U,graph->props.N);
		uprintvector(graph->props.V,graph->props.N);*/
	}

}

void matprod(std::vector<std::vector<int>> &result, const Bitvector &imp, const int N)
{
  std::vector<std::vector<int>> prod;
  prod.resize(N);
  std::vector<int> zero(N,0);

    for(int j=0;j<N;j++)
      prod[j]=zero;

    for(int j=0;j<N;j++)
    {

      /// I need the product of the matrices. The labels actually don't change anything!
      for(int k=0;k<N;k++)
      {
        for(int i=k+1;i<N;i++)
        {
        if(imp[k].read(i))
          {
            prod[k][j]+=result[i][j];
          }
      }

    }
  }

  result=prod;
}

void measure_correlators(const std::vector<int> &spins,const Bitvector &link, const int N, double *corr, const std::vector<unsigned int> U)
{


  for(int i=0;i<N;i++) corr[i]=0;

//// here be an efficient way to get the product of all the link matrices. what I want to calculate is
//// <s^n> = sum_{i,j} s_i s_j L_{i j}^n  /// I guess this means I want a floating point array oh well
  std::vector<std::vector<int>> prod;
  std::vector<int> temp;
  std::vector<int> zero(N,0);

  prod.resize(N);
  for(int i=0;i<N;i++) /// initializing the product as a unit matrix
  {
    temp=zero;
    temp[i]=1;
    prod[i]=temp;
  }



  for(int i=0;i<N;i++)
  {
    matprod(prod,link,N);

    for(int j=0;j<N;j++)
    {
      for(int k=j;k<N;k++)
      {
      corr[i]+=prod[U[j]][U[k]]*spins[k]*spins[j];
      }
    }
  }

}


void IsingObservables(std::vector<int> &spins, Bitvector &adj, const int N, double & relcorr, double &magnetisation,std::vector<unsigned int> &U)
{
  magnetisation=0;
  relcorr=0;
  for(int i=0;i<N;i++)
  {
   magnetisation+=(double)spins[i];
   for(int j=0; j<N;j++)
    {

    if(adj[U[i]].read(U[j]))
    {
      relcorr+=spins[i]*spins[j];   /// what I need is spins^T.links.spins
    }

  }

  }
  relcorr/=N;
  magnetisation/=N;
//  std::cout<<magnetisation<<std::endl;
//  std::cout<<relcorr<<std::endl;
}

void randomSpinState(std::vector<int> &spin, MersenneRNG &mrng, const int N)
{
	#if UNIT_TEST
	assert (N > 0);
	#endif

	spin.resize(N);
	for (int i = 0; i < N; i++)
	{
		//spin[i]=rand()%2;
		//if(spin[i]==0) spin[i]=-1;
		spin[i] = (static_cast<int>(mrng.urng() * 2.0) << 1) - 1;
	}
}

/// this calculates the ising action
// algorithm is not optimized, but what gives
bool IsingAction(std::vector<int> &spins, double &Iaction, Bitvector &link, const int N, const double J,std::vector<unsigned int> &U)
{
	#if UNIT_TEST
	assert (spins.size() > 0);
	assert (link.size() > 0);
	assert (N > 0);
	#endif

	int IA = 0.0;
	for(int i = 0; i < N; i++) {
		int si = spins[i];
		for(int j = i + 1; j < N; j++)
			if(link[U[i]].read(U[j]))
				IA += si * spins[j];
	}
	Iaction = J * static_cast<double>(IA);

	return true;
}

//Measure Causal Set Action
//Algorithm has been optimized using minimal bitwise operations
//Requires the existence of the whole adjacency matrix
//This will calculate all cardinality intervals by construction
bool measureAction_v3(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon)
{
	#if UNIT_TEST
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
	#if UNIT_TEST
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
	#if UNIT_TEST
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
	for(int i = 0; i < N; i++)
		mat[i].printBitset();
}


void printvector(std::vector<int> &mat, const int N)
{
	for(int i = 0; i < N; i++)
		std::cout << mat[i] << " ";
	std::cout << std::endl;
}
void uprintvector(std::vector<unsigned int> &mat, const int N)
{
	for(int i = 0; i < N; i++)
		std::cout << mat[i] << " ";
	std::cout << std::endl;
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

void updateLinks(Bitvector &new_adj, Bitvector &new_link, const int N)
{
	#if UNIT_TEST
	assert (new_adj.size() > 0);
	assert (new_link.size() > 0);
	assert (N > 0);
	#endif

	for (int i = 0; i < N; i++)
		new_link[i].reset();

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			//Continue only if two elements are related
			if (!new_adj[i].read(j)) continue;

			//Check if their open Alexandroff set is size 0
			if (!new_adj[i].partial_vecprod(new_adj[j], i, j - i + 1)) {
				new_link[i].set(j);
				new_link[j].set(i);
			}
		}
	}
}
