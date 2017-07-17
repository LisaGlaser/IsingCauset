#include "subroutines.h"
#include "CuResources.h"

int main(int argc, char **argv)
{
	Memory mem = Memory();
	CausetPerformance cp = CausetPerformance();
	stopwatchStart(&cp.sCauset);

	Graph graph = Graph(parseArgs(argc, argv));

	int e_start = printStart((const char**)argv, 0);
	bool success = false;

	if (!init(&graph, &mem, &cp)) goto Exit;
	if (!evolve(&graph, &mem, &cp)) goto Free;
	if (graph.props.flags.print && !printGraph(&graph)) goto Free;

	Free:
	destroyGraph(&graph, mem.used);

	success = !mem.used;
	if (!success) fprintf(stderr, "WARNING: Memory leak detected!\n");

	Exit:
	stopwatchStop(&cp.sCauset);
	printFinish((const char**)argv, e_start, 0, success ? PASSED : FAILED);
	printf("Time: %5.6f sec\n", cp.sCauset.elapsedTime);
	printf("PROGRAM COMPLETED\n\n");
	fflush(stdout);

	return 0;
}

Properties parseArgs(int argc, char **argv)
{
	Properties props = Properties();

	int c, longIndex;
	//Single-character options
	static const char *optString = ":b:e:f:hn:pS:s:v";
	//Multi-character options
	static const struct option longOpts[] = {
		{ "beta",		required_argument,	NULL, 'b' },
		{ "epsilon",		required_argument,	NULL, 'e' },
		{ "filename",		required_argument,	NULL, 'f' },
		{ "help",		no_argument,		NULL, 'h' },
		{ "nodes",		required_argument,	NULL, 'n' },
		{ "print",		no_argument,		NULL, 'p' },
		{ "print-edges",	no_argument,		NULL,  0  },
		{ "seed",		required_argument,	NULL, 's' },
		{ "sweeps",		required_argument,	NULL, 'S' },
		{ "verbose",		no_argument,		NULL, 'v' },
		{ NULL,			0,			0,     0  }
	};

	try {
		while ((c = getopt_long(argc, argv, optString, longOpts, &longIndex)) != -1) {
			switch (c) {
			case 'b':	//Inverse temperature
				props.beta = atof(optarg);
				//if (props.beta <= 0.0 || props.beta >= 1.0)
				//	throw CausetException("Invalid argument for '--beta' parameter!");
				break;
			case 'e':	//Epsilon
				props.epsilon = atof(optarg);
				if (props.epsilon <= 0.0 || props.epsilon > 1.0)
					throw CausetException("Invalid argument for '--epsilon' parameter!");
				break;
			case 'f':	//Output filename
				props.filename = std::string(optarg);
				break;
			//case 'h' located at the end
			case 'n':	//Number of nodes
				props.N = atoi(optarg);
				if (props.N <= 0)
					throw CausetException("Invalid argument for '--nodes' parameter!");
				break;
			case 'p':	//Print simulation results to file
				props.flags.print = true;
				break;
			case 'S':	//Number of sweeps
				props.sweeps = atoi(optarg);
				if (props.sweeps <= 0)
					throw CausetException("Invalid argument for '--sweeps' parameter!");
				break;
			case 's':	//Random seed
				props.seed = atol(optarg);
				if (props.seed <= 0L)
					throw CausetException("Invalid argument for '--seed' parameter!");
				break;
			case 'v':	//Verbose output
				props.flags.verbose = true;
				break;
			case 0:
				if (!strcmp("print-edges", longOpts[longIndex].name))
					//Print edge list to file
					props.flags.print_edges = true;
				else {
					//Unrecognized options
					fprintf(stderr, "Option --%s is not recognized.\n", longOpts[longIndex].name);
					exit(4);
				}
				break;
			case 'h':
				//Print help menu
				printf("\nUsage  :  Evolution [options]\n\n");
				printf("Evolution Options...................\n");
				printf("====================================\n");
				printf("Flag:\t\t\tMeaning:\t\t\tSuggested Values:\n");
				printf("  -b, --beta\t\tInverse temperature\t\t0.3\n");
				printf("  -e, --epsilon\t\tSmearing Parameter\t\t(0, 1]\n");
				printf("  -f, --filename\tOutput filename\n");
				printf("  -h, --help\t\tDisplay this menu\n");
				printf("  -n, --nodes\t\tNumber of nodes\t\t\t100\n");
				printf("  -p, --print\t\tPrint results to file\n");
				printf("      --print-edges\tPrint edge list to file\n");
				printf("  -s, --seed\t\tRandom seed\t\t\t18100\n");
				printf("  -S, --sweeps\t\tNumber of sweeps\t\t1000\n");
				printf("  -v, --verbose\t\tVerbose output\n");
				printf("\n");

				printf("Report bugs to w.cunningham@northeastern.edu\n");
				printf("Bitbucket repository home page: <...>\n");
				exit(0);
			case ':':
				//Single-character flag needs an argument
				if (!!optopt)
					fprintf(stderr, "%s : option '-%c' requires an argument.\n", argv[0], optopt);
				else
					fprintf(stderr, "%s : option '%s' requires an argument.\n", argv[0], argv[optind-1]);
				exit(2);
			case '?':	//Unrecognized flag
			default:	//Default case
				if (!!optopt)
					fprintf(stderr, "%s : option '-%c' is not recognized.\n", argv[0], optopt);
				else
					fprintf(stderr, "%s : option '%s' is not recognized.\n", argv[0], argv[optind-1]);
				exit(3);
			}
		}
	} catch (CausetException c) {
		fprintf(stderr, "CausetException in %s on line %d: %s\n", __FILE__, __LINE__, c.what());
		exit(1);
	}

	//Initialize RNG
	if (props.seed == 12345L) {
		srand(time(NULL));
		props.seed = static_cast<long>(time(NULL));
	}

	props.mrng.urng.engine().seed(props.seed);
	props.mrng.urng.distribution().reset();

	for (int i = 0; i < RNG_THERMAL; i++)
		props.mrng.urng();

	return props;
}

bool init(Graph * const graph, Memory * const mem, CausetPerformance * const cp)
{	
	#if UNIT_TEST
	assert (graph != NULL);
	assert (mem != NULL);
	assert (cp != NULL);
	#endif

	#ifdef _OPENMP
	printf("\n\t[ *** OPENMP MODULE ACTIVE *** ]\n");
	#endif

	#ifdef AVX2_ENABLED
	printf("\n\t[ ***  AVX2 MODULE ACTIVE  *** ]\n");
	#endif

	printf("\nInitializing Network...\n");
	fflush(stdout);

	//randomTotalOrder(graph->props.U, graph->props.N);
	//randomTotalOrder(graph->props.V, graph->props.N);
	graph->props.U.resize(graph->props.N);
	std::iota(graph->props.U.begin(), graph->props.U.end(), 0);
	graph->props.V.resize(graph->props.N);
	std::iota(graph->props.V.begin(), graph->props.V.end(), 0);
	mem->used += sizeof(unsigned int) * graph->props.N * 2;
	/*printf("U: [ ");
	for (int i = 0; i < graph->props.N; i++)
		printf("%d ", graph->props.U[i]);
	printf("]\n");

	printf("V: [ ");
	for (int i = 0; i < graph->props.N; i++)
		printf("%d ", graph->props.V[i]);
	printf("]\n");*/

	graph->adj.reserve(graph->props.N);
	graph->new_adj.reserve(graph->props.N);
	for (int i = 0; i < graph->props.N; i++) {
		FastBitset fb(graph->props.N);
		graph->adj.push_back(fb);
		graph->new_adj.push_back(fb);
		mem->used += sizeof(BlockType) * fb.getNumBlocks() * 2;
	}

	updateRelations(graph->adj, graph->props.U, graph->props.V, graph->props.N);

	try {
		graph->obs.action_data = (float*)calloc(graph->props.sweeps, sizeof(float));
		if (graph->obs.action_data == NULL)
			throw std::bad_alloc();
		mem->used += sizeof(float) * graph->props.sweeps;
	} catch (std::bad_alloc) {
		fprintf(stderr, "Failed to allocate memory in %s at line %d.\n", __FILE__, __LINE__);
		return false;
	}

	printf("\tTask Completed.\n");
	fflush(stdout);

	return true;
}

bool evolve(Graph * const graph, Memory * const mem, CausetPerformance * const cp)
{
	printf("Evolving Network...\n");
	fflush(stdout);

	Bitvector workspace;
	uint64_t npairs = static_cast<uint64_t>(graph->props.N) * (graph->props.N - 1) / 2;
	int n = graph->props.N + graph->props.N % 2;
	uint64_t np = static_cast<uint64_t>(n) * (n - 1) / 2;
	double new_action = 0.0;
	uint64_t clone_length = graph->adj[0].getNumBlocks();
	unsigned int stdim = 2;
	double dS = 0.0;

	try {
		graph->obs.cardinalities = (uint64_t*)malloc(sizeof(uint64_t) * graph->props.N * omp_get_max_threads());
		if (graph->obs.cardinalities == NULL)
			throw std::bad_alloc();
		memset(graph->obs.cardinalities, 0, sizeof(uint64_t) * graph->props.N * omp_get_max_threads());
		mem->used += sizeof(uint64_t) * graph->props.N * omp_get_max_threads();

		//Allocate memory for workspace
		//This avoids many internal malloc/free inside parallel loop
		workspace.reserve(omp_get_max_threads());
		for (int i = 0; i < omp_get_max_threads(); i++) {
			FastBitset fb(static_cast<uint64_t>(graph->props.N));
			workspace.push_back(fb);
			mem->used += sizeof(BlockType) * fb.getNumBlocks();
		}
	} catch (std::bad_alloc) {
		fprintf(stderr, "Memory allocation failure in %s on line %d!\n", __FILE__, __LINE__);
		return false;
	}

	if (!measureAction_v3(graph->obs.cardinalities, graph->obs.action, graph->adj, workspace, stdim, graph->props.N, graph->props.epsilon))
		return false;

	std::ofstream data;
	std::stringstream sstm;
	sstm << "dat/int/";
	if (!graph->props.filename.compare(""))
		sstm << "intervals.evo.act.dat";
	else
		sstm << graph->props.filename;
	data.open(sstm.str().c_str());
	if (!data.is_open())
		return false;
	for (int s = 0; s < graph->props.sweeps; s++) {
		for (uint64_t k = 0; k < npairs; k++) {
			//Pick a random pair in U
			uint64_t v = static_cast<uint64_t>(graph->props.mrng.urng() * (np - 1)) + 1;

			int i = static_cast<int>(v / (n - 1));
			int j = static_cast<int>(v % (n - 1) + 1);
			int do_map = i >= j;
			i += do_map * ((((n >> 1) - i) << 1) - 1);
			j += do_map * (((n >> 1) - j) << 1);

			if (j == graph->props.N) continue;

			//Swap elements at indices i and j in U
			/*graph->props.U[i] ^= graph->props.U[j];
			graph->props.U[j] ^= graph->props.U[i];
			graph->props.U[i] ^= graph->props.U[j];*/

			std::vector<unsigned int> &W = graph->props.mrng.urng() < 0.5 ? graph->props.U : graph->props.V;
			W[i] ^= W[j];
			W[j] ^= W[i];
			W[i] ^= W[j];

			/*printf("\nU:");
			for (int i = 0; i < graph->props.N; i++)
				printf(" %d", graph->props.U[i]);
			printf("\n");
			printf("V:");
			for (int i = 0; i < graph->props.N; i++)
				printf(" %d", graph->props.V[i]);
			printf("\n");*/

			//Construct the new adjacency matrix
			//Optimize this later
			updateRelations(graph->new_adj, graph->props.U, graph->props.V, graph->props.N);

			if (!measureAction_v3(graph->obs.cardinalities, new_action, graph->new_adj, workspace, stdim, graph->props.N, graph->props.epsilon))
				return false;
			//if (!measureAction_v2(graph->obs.cardinalities, new_action, graph->new_adj, workspace, stdim, graph->props.N, graph->props.epsilon))
			//	return false;
			//if (!measureAction_v1(graph->obs.cardinalities, new_action, graph->props.U, graph->props.V, stdim, graph->props.N, graph->props.epsilon))
			//	return false;
			/*printf("\nC:");
			for (int i = 0; i < graph->props.N; i++)
				printf(" %d", (int)graph->obs.cardinalities[i]);
			printf("\n");
			if (!measureAction_v3(graph->obs.cardinalities, new_action, graph->new_adj, workspace, stdim, graph->props.N, graph->props.epsilon))
				return false;
			printf("D:");
			for (int i = 0; i < graph->props.N; i++)
				printf(" %d", (int)graph->obs.cardinalities[i]);
			printf("\n");*/

			dS = graph->props.beta * (new_action - graph->obs.action);
			if (dS < 0 || exp(-dS) > graph->props.mrng.urng()) {	//Accept change*/
				for (int m = 0; m < graph->props.N; m++)
					graph->new_adj[m].clone(graph->adj[m]);
				graph->obs.action = new_action;
			} else {	//Reject change
				W[i] ^= W[j];
				W[j] ^= W[i];
				W[i] ^= W[j];
			}

			//if (k == 5) return true;
		}
		graph->obs.action_data[s] = graph->obs.action;
		for (int i = 0; i < graph->props.N; i++)
			data << graph->obs.cardinalities[i] << " ";
		data << "\n";
		//printf("sweep [%d] action: %f\n", s, graph->obs.action);
	}
	data.flush();
	data.close();

	//Free Workspace
	mem->used -= sizeof(BlockType) * clone_length * omp_get_max_threads();
	workspace.clear();
	workspace.swap(workspace);

	printf("\tTask Completed.\n");
	fflush(stdout);

	return true;
}

bool printGraph(Graph * const graph)
{
	printf("Printing data to file.\n");

	std::ofstream os;
	std::stringstream sstm;
	sstm << "dat/act/";
	if (!graph->props.filename.compare(""))
		sstm << "action.evo.act.dat";
	else
		sstm << graph->props.filename;
	os.open(sstm.str().c_str());
	if (!os.is_open())
		return false;
	for (int i = 0; i < graph->props.sweeps; i++)
		os << graph->obs.action_data[i] << "\n";
	os.flush();
	os.close();

	/*os.open("dat/int/intervals.evo.act.dat");
	if (!os.is_open())
		return false;
	for (int i = 0; i < graph->props.N; i++)
		os << graph->obs.cardinalities[i] << "\n";
	os.flush();
	os.close();*/

	return true;
}

void destroyGraph(Graph * const graph, size_t &used)
{
	graph->props.U.clear();
	graph->props.U.swap(graph->props.U);

	graph->props.V.clear();
	graph->props.V.swap(graph->props.V);

	used -= sizeof(unsigned int) * graph->props.N * 2;

	used -= 2 * sizeof(BlockType) * graph->adj[0].getNumBlocks() * graph->props.N;
	graph->adj.clear();
	graph->adj.swap(graph->adj);

	graph->new_adj.clear();
	graph->new_adj.swap(graph->new_adj);

	free(graph->obs.cardinalities);
	graph->obs.cardinalities = NULL;
	used -= sizeof(uint64_t) * graph->props.N * omp_get_max_threads();

	free(graph->obs.action_data);
	graph->obs.action_data = NULL;
	used -= sizeof(float) * graph->props.sweeps;

	return;
}