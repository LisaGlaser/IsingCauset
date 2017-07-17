#ifndef SUBROUTINES_H_
#define SUBROUTINES_H_

#include "causet.h"
#include "operations.h"
#include "CuResources.h"

void randomTotalOrder(std::vector<unsigned int> &U, const int N);

bool measureAction_v3(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon);

bool measureAction_v2(uint64_t * const cardinalities, double &action, Bitvector &adj, Bitvector &workspace, const unsigned int stdim, const int N, const double epsilon);

bool measureAction_v1(uint64_t * const cardinalities, double &action, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const unsigned int stdim, const int N, const double epsilon);

void updateRelations(Bitvector &new_adj, const std::vector<unsigned int> U, const std::vector<unsigned int> V, const int N);

#endif
