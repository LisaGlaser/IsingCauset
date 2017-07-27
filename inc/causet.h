#ifndef CAUSET_H_
#define CAUSET_H_

//STL
#include <getopt.h>

//OpenMP
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

//Boost
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//Custom Files
#include <FastBitset.h>
#include <FastMath.h>
#include <stopwatch.h>
#include <printcolor.h>

//Local Files
#include "constants.h"

//Adjacency Matrix
typedef std::vector<FastBitset> Bitvector;

//Boost RNG
typedef boost::mt19937 Engine;
typedef boost::uniform_real<double> UDistribution;
typedef boost::variate_generator<Engine&, UDistribution> UGenerator;

struct MersenneRNG {
	MersenneRNG() : udist(0.0, 1.0), urng(eng, udist) {}

	Engine eng;
	UDistribution udist;
	UGenerator urng;
};

//Memory Resources
struct Memory {
	Memory() : used(0), max(0) {}

	size_t used;
	size_t max;
};

struct Flags {
	Flags() : print(false), print_edges(false), verbose(false), bench(false) {}

	bool print;
	bool print_edges;

	bool verbose;
	bool bench;
};


struct Properties {
	Properties() : flags(Flags()), N(0), sweeps(0), beta(0.0), epsilon(1.0), Jising(0.0), seed(12345L), graphID(0), mrng(MersenneRNG()), filename("") {}

	Flags flags;

	int N;
	int sweeps;

	std::vector<unsigned int> U;
	std::vector<unsigned int> V;

	double beta;
	double epsilon;
	double Jising=1.0;

	long seed;
	int graphID;

	MersenneRNG mrng;

	std::string filename;
};

struct Observables {
	Observables() : k(0.0), cardinalities(NULL), action(0.0), Iaction(0.0), action_data(NULL) {}

	float k;

	uint64_t *cardinalities;
	double action;
	double Iaction;

	float *action_data;
};

struct Graph {
	Graph() : props(Properties()), obs(Observables()) {}
	Graph(Properties _props) : props(_props), obs(Observables()) {}

	Properties props;
	Observables obs;
	Bitvector adj;				//Adjacency Matrix
	Bitvector new_adj;			//Updated Adjacency Matrix

	Bitvector link;				//Link Matrix
	Bitvector new_link;			//Updated Link Matrix

	std::vector<int> spins;  // Ising model spins
	std::vector<int> new_spins; //Updated spins

	std::vector<unsigned int> k_in;		//In-Degrees
	std::vector<unsigned int> k_out;	//Out-Degrees
};

//Algorithmic Performance
struct CausetPerformance {
	CausetPerformance() : sCauset(Stopwatch()), sMeasureAction(Stopwatch()) {}

	Stopwatch sCauset;
	Stopwatch sMeasureAction;
};

//Custom exception class used in this program
class CausetException : public std::exception
{
public:
	CausetException() : msg("Unknown Error!") {}
	explicit CausetException(char const * _msg) : msg(_msg) {}
	virtual ~CausetException() throw () {}
	virtual const char * what() const throw () { return msg; }

protected:
	char const * msg;
};

Properties parseArgs(int argc, char **argv);

bool init(Graph * const graph, Memory * const mem, CausetPerformance * const cp);

bool evolve(Graph * const graph, Memory * const mem, CausetPerformance * const cp);

bool printGraph(Graph * const graph);

void destroyGraph(Graph * const graph, size_t &used);

#endif
