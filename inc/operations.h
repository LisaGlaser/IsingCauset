#ifndef OPERATIONS_H_
#define OPERATIONS_H_

#include "causet.h"
#include "CuResources.h"

//Calculate the action from the abundancy intervals
inline double calcAction(const uint64_t * const cardinalities, const int stdim, const double epsilon)
{
	#if DEBUG
	assert (cardinalities != NULL);
	assert (stdim >= 2 && stdim <= 4);
	assert (epsilon > 0.0 && epsilon <= 1.0);
	#endif

	long double action = 0.0;

	if (epsilon < 1.0) {
		long double eps1 = static_cast<long double>(epsilon) / (1.0 - epsilon);
		long double eps2 = eps1 * eps1;
		long double eps3 = eps2 * eps1;
		long double epsi = 1.0;
		long double c3_1 = 27.0 / 8.0;
		long double c3_2 = 9.0 / 8.0;
		long double c4_1 = 4.0 / 3.0;
		long double ni;
		uint64_t i;

		for (i = 0; i < cardinalities[0] - 3; i++) {
			ni = static_cast<long double>(cardinalities[i+1]);
			if (stdim == 2)
				action += ni * epsi * (1.0 - 2.0 * eps1 * i + 0.5 * eps2 * i * (i - 1.0));
			else if (stdim == 3)
				action += ni * epsi * (1.0 - c3_1 * i * eps1 + c3_2 * i * (i - 1.0) * eps2);
			else if (stdim == 4)
				action += ni * epsi * (1.0 - 9.0 * eps1 * i + 8.0 * eps2 * i * (i - 1.0) - c4_1 * eps3 * i * (i - 1.0) * (i - 2.0));
			else
				action = NAN;
			epsi *= (1.0 - epsilon);
		}

		if (stdim == 2)
			action = 4.0 * epsilon * (cardinalities[0] - 2.0 * epsilon * action);
		else if (stdim == 3)
			action = pow(M_PI / (3.0 * sqrt(2.0)), 2.0 / 3.0) / GAMMA(5.0 / 3.0, STL) * (pow(epsilon, 2.0 / 3.0) * cardinalities[0] - pow(epsilon, 5.0 / 3.0) * action);
		else if (stdim == 4)
			action = (4.0 / sqrt(6.0)) * (sqrt(epsilon) * cardinalities[0] - POW(epsilon, 1.5, STL) * action);
		else
			action = NAN;
	} else {
		if (stdim == 2)
			action = 4.0 * (cardinalities[0] - 2.0 * (cardinalities[1] - 2.0 * cardinalities[2] + cardinalities[3]));
		else if (stdim == 3)
			action = pow(M_PI / (3.0 * sqrt(2.0)), 2.0 / 3.0) / GAMMA(5.0 / 3.0, STL) * (cardinalities[0] - cardinalities[1] + (27.0 / 8.0) * cardinalities[2] - (9.0 / 4.0) * cardinalities[3]);
		else if (stdim == 4)
			action = (4.0 / sqrt(6.0)) * (cardinalities[0] - cardinalities[1] + 9.0 * cardinalities[2] - 16.0 * cardinalities[3] + 8.0 * cardinalities[4]);
		else
			action = NAN;
	}

	return static_cast<double>(action);
}

#endif
