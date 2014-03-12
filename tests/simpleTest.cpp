#include <complex_bessel.h>
#include <iostream>

using namespace sp_bessel;

int main(int argc, char* argv[])
{
	std::complex<double> z(1.0,0.0);
	std::cout << besselJ(2,z) << std::endl;
	std::cout << besselJp(2,z) << std::endl;
	std::cout << besselY(2,z) << std::endl;
	std::cout << besselYp(2,z) << std::endl;
	std::cout << besselK(2,z) << std::endl;
	std::cout << besselKp(2,z) << std::endl;
	std::cout << besselH1(2,z) << std::endl;
	std::cout << besselH1p(2,z) << std::endl;
	return 0;
}