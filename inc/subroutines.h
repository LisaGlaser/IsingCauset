#ifndef SUBROUTINES_H_
#define SUBROUTINES_H_

#include "causet.h"
#include "operations.h"
#include "CuResources.h"

void randomTotalOrder(std::vector<unsigned int> &U, const int N);

void randomSpinState(std::vector<int> &spin, const int N);

bool IsingAction(std::vector<int> &spins, double &Iaction, Bitvector &link, const int N, const double J);

bool measureAction_v3(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon);

bool measureAction_v2(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon);

bool measureAction_v1(uint64_t * const cardinalities, double &action, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const unsigned int stdim, const int N, const double epsilon);

void updateRelations(Bitvector &new_adj, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const int N);
void updateLinks(Bitvector &new_adj, Bitvector &new_link, const int N);

void printmatrix(Bitvector &mat, const int N);
void printvector(std::vector<int> &mat, const int N);
#endif
