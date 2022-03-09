#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

const double M = 1; //mass for the particles (Kg)
const double K1 = 1; //elastic constant for the springs attached to the walls (Kg/s^2)

void everything(double K2, std::string datos);

int main (void)
{
  for (double K2 = 0.1; K2<=2.0; K2+=0.1) {
    if (K2==0.1) {everything(K2, "datosratio01.txt");}
    if (K2==0.5) {everything(K2, "datosratio05.txt");}
    if (K2==1.0) {everything(K2, "datosratio10.txt");}
    if (K2==1.5) {everything(K2, "datosratio15.txt");}
    if (K2==2.0) {everything(K2, "datosratio20.txt");}
  }
  return 0;
}

void everything(double K2, std::string datos)
{
  std::ofstream fout(datos);
  fout.precision(12); fout.setf(std::ios::scientific);
  Eigen::MatrixXd dif_eq;
  for(int exponent_of_two = 0; exponent_of_two<10;exponent_of_two++) {
    int N = std::pow(2,exponent_of_two);
    dif_eq = Eigen::MatrixXd::Constant(N, N, 0);
    for(int ii = 0; ii<N ; ii++) {
      dif_eq(ii,ii) = 2*K2;
      if (ii == 0 or ii == N-1) {
	dif_eq(ii,ii) = K1+K2;
      }
      if (ii<N-1) {
	dif_eq(ii,ii+1) = -K2;
	dif_eq(ii+1,ii) = -K2;
      }
    }
    dif_eq = dif_eq*(1/M);
    Eigen::VectorXcd eivals = dif_eq.eigenvalues();
    for(int value = 0; value<N; value ++) {
      fout<<N<<"\t";
      fout<<std::real(eivals(value,0))<<"\t";
      fout<<std::pow(std::real(eivals(value,0)), 0.5)<<std::endl;
    }
  }
  fout.close();
}
