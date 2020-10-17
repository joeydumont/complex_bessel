#include <armadillo>
#include <complex>
#include <complex_bessel.h>
#include <iostream>

using namespace sp_bessel;

int
main()
{
  int       maxValue = 500;
  arma::mat realAiry(2 * maxValue + 1, 2 * maxValue + 1);
  arma::mat imagAiry(2 * maxValue + 1, 2 * maxValue + 1);

  for (int i = 0; i < 2 * maxValue + 1; i++) {
    int m = i - maxValue;
    for (int j = 0; j < 2 * maxValue + 1; j++) {
      int                  n         = j - maxValue;
      std::complex<double> airyValue = airy(std::complex<double>(m, n));
      realAiry(i, j)                 = std::real(airyValue);
      imagAiry(i, j)                 = std::imag(airyValue);
    }
  }

  realAiry.save("realAiry.dat", arma::raw_ascii);
  imagAiry.save("imagAiry.dat", arma::raw_ascii);
  return 0;
}
