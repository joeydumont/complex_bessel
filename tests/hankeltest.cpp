/* This example program computes the value
 * of the H_0(z) function of the first kind is computed in the
 * complex plane and the result is displayed and compared 
 * to a "known" result
 * 
 * Joey Dumont <joey.dumont@gmail.com>
 * Denis Gagnon <gagnon88@gmail.com>
 * 
 */


#include <complex>
#include <iostream>
#include <complex_bessel/sp_hankel.h>

using namespace std;

int main()
{

  sp_Hankel myHankel(1, 0.0);
  
  complex<double> Z(1.45,8.45);
  
  complex<double> Y = myHankel(Z);
  
  cout << "Computing the value of the Hankel function of "
  << "the first kind, order 0" << endl;
  cout << "at z = " << real(Z) << " + " << imag(Z) << "i" << endl;
  cout << endl;
  cout << "The result should be " << endl;
  cout << "5.7455e-05 - 2.1872e-06i" << endl;
 
  cout << endl;
  cout << "Result" << endl; 
  cout << "Real part: " << real(Y) << endl;
  cout << "Imaginary part: " << imag(Y) << endl;
  
  return 0;
  
}
 
