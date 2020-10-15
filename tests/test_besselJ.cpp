#include "test_besselJ.h"

#include <sstream>

std::complex<double>
besselJAddition1(std::complex<double> z, int kMax)
{
  std::complex<double> sum(0.0, 0.0);

  sum += pow(besselJ(0, z), 2.0);

  for (int i = 1; i < kMax; i++) {
    sum += 2.0 * pow(besselJ(i, z), 2.0);
  }

  return sum;
}

std::complex<double>
besselJAddition2(std::complex<double> z, int N, int kMax)
{
  std::complex<double> sum(0.0, 0.0);

  for (int i = 0; i < 2 * N; i++) {
    sum += pow(-1.0, i) * besselJ(i, z) * besselJ(2 * N - i, z);
  }

  for (int i = 1; i < kMax; i++) {
    sum += 2.0 * besselJ(i, z) * besselJ(2 * N + i, z);
  }

  return sum;
}

std::complex<double>
besselJWronskian1(int order, std::complex<double> z)
{
  return besselJ(order + 1, z) * besselJ(-order, z) + besselJ(order, z) * besselJ(-(order + 1), z);
}
