#ifndef SUBROUTINES_H_
#define SUBROUTINES_H_

#include "causet.h"
#include "operations.h"
#include "CuResources.h"
#include <iterator>

void initialState(Graph * const graph);
void randomTotalOrder(std::vector<unsigned int> &U, const int N);

void randomSpinState(std::vector<int> &spin, MersenneRNG &mrng, const int N);
void measure_correlators(const std::vector<int> &spins, const Bitvector &link, const int N, double *corr, const std::vector<unsigned int> U);
bool IsingAction(std::vector<int> &spins, double &Iaction, Bitvector &link, const int N, const double J,std::vector<unsigned int> &U);
void IsingObservables(std::vector<int> &spins, Bitvector &link, const int N, double &relcorr, double &magnetisation,std::vector<unsigned int> &U);

bool measureAction_v3(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon);

bool measureAction_v2(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon);

bool measureAction_v1(uint64_t * const cardinalities, double &action, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const unsigned int stdim, const int N, const double epsilon);

void updateRelations(Bitvector &new_adj, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const int N);
void updateLinks(Bitvector &new_adj, Bitvector &new_link, const int N);

void matprod(std::vector<std::vector<int>> &result, Bitvector &input, int N);
void printmatrix(Bitvector &mat, const int N);
void printvector(std::vector<int> &mat, const int N);
void uprintvector(std::vector<unsigned int> &mat, const int N);
#endif
