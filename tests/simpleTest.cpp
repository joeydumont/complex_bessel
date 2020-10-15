#include <complex_bessel.h>
//#include <complex_bessel_bits/sph_besselFunctions.h>
#include <iostream>

using namespace sp_bessel;

int
main(int argc, char* argv[])
{
  std::complex<double> z(1.0, 0.0);
  std::complex<double> nearOrigin(1.0e-10, 0);
  std::cout << besselJ(2, z) << std::endl;
  std::cout << besselJp(2, z) << std::endl;
  std::cout << besselJp(2, z, 2) << std::endl;
  std::cout << besselJp(2, z, 3) << std::endl;
  std::cout << besselY(2, z) << std::endl;
  std::cout << besselYp(2, z) << std::endl;
  std::cout << besselK(2, z) << std::endl;
  std::cout << besselKp(2, z, 2) << std::endl;
  std::cout << hankelH1(2, z) << std::endl;
  std::cout << hankelH1p(2, z) << std::endl;
  std::cout << sph_besselJ(-constants::pi, nearOrigin) << std::endl;
  std::cout << sph_besselY(-constants::pi, nearOrigin) << std::endl;
  std::cout << sph_hankelH1(-constants::pi, nearOrigin) << std::endl;
  std::cout << sph_hankelH2(-constants::pi, nearOrigin) << std::endl;
  return 0;
}
