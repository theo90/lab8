// Relativistic oscillator
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
const int dim = 2;
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,
            double* const k1, double* const k2, double* const k3, double* const k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
	double dx = 0.1,x;
	const double L = 15;
  double y0[dim];
  double yn[dim];
  double k1[dim], k2[dim], k3[dim], k4[dim];

  for(double x0 = 0.1; x0< 5; x0+=0.1){
     y0[0] = x0;
     y0[1] = 0;
     x = 0;


  	while(x<L)
  	{
  		x += dx;
  		RKstep(yn, y0, x, dx, k1, k2, k3, k4);

      // Check whether dp/dt changed sign in last step
      if ( y0[1]>0 && yn[1]<0 ) break;

      // otherwise continue
      for(int i=0; i<dim; i++) y0[i] = yn[i];
  	}

    // find exactly the point at which dp/dt = 0
    // remember, up to here y0[1] > 0
    double thetaL = 0, thetaR = 1, theta = 0.5;
    double h = 1;
    double b1, b2, b4;
    while( abs(h) > 1e-8){
      b1 = theta - 3.*theta*theta/2 + 2./3*pow(theta,3);
      b2 = pow(theta,2) - 2./3 * pow(theta,3);
      b4 = - pow(theta,2)/2 + 2./3 * pow(theta,3);

      // calculate interpolated dp/dt at point theta
      h = y0[1] + dx * ( b1 * k1[1] + b2*(k2[1] +k3[1]) + b4*k4[1] );

      if (h>0)
         thetaL = theta;
      else
         thetaR = theta;

      theta = (thetaR + thetaL) * 0.5;
    }

    out << x0 << "\t" << x + theta * dx - dx << endl;
  }
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx,
            double* const k1, double* const k2, double* const k3, double* const k4)
{


  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
  double h = y0[0];
  y0[0] = y0[1];
  y0[1] = - h/sqrt(1+h*h);
}
