#include <armadillo>
#include <iostream>
#include <sstream>

#include "test_besselJ.h"

using namespace arma;

// Main function
int
main(int argc, char* argv[])
{
  // Arguments
  double zrMax = atof(argv[1]);
  double ziMax = atof(argv[2]);
  int    kMax  = atoi(argv[3]);
  int    N     = atoi(argv[4]);

  // Setup variables.
  colvec zr = linspace(0.0, zrMax, N);
  colvec zi = linspace(0.0, ziMax, N);

  // To store the tolerance achieved given the number of terms.
  mat resultsAdd(N * N, 4);

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      // Check the addition theorems
      resultsAdd(i * N + j, 0) = zr(i);
      resultsAdd(i * N + j, 1) = zi(j);
      resultsAdd(i * N + j, 2) = abs(besselJAddition1(std::complex<double>(zr(i), zi(j)), kMax) - 1.0);
      resultsAdd(i * N + j, 3) = abs(besselJAddition2(std::complex<double>(zr(i), zi(j)), N, kMax));
    }
  }

  std::stringstream ss;
  ss << "resultsAdd-kMax" << kMax << "-N" << N;
  resultsAdd.save(ss.str(), raw_ascii);

  return 0;
}
