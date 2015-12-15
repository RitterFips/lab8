// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx, double* k1, double* k2, double* k3, double* k4);
void inter(double* k1, double* k2,double* k3, double* k4, double p0, double& pn, double theta, double x);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
	ofstream bisec("bisection");
  const int dim = 2;
	double dx = 0.1,x=0;
	const double L = 20;
	
	double y0[dim] = {1.0, 0.0};
	double yn[dim];
	int n = 0;
	double k1[dim], k2[dim], k3[dim], k4[dim];
	double theta = 0.5;
	double p0 = y0[1];
	double pn = 0;
	double thetar = 0;
	double thetal = 1;
	double a = 4;
	double tol = 1e-10;
	
  //double temp[dim] = {0.0, 0.0};
  out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << endl;
	for(int i = 0; i < 50; i++ ){
	y0[0] = 0.0 + i*0.1;
	y0[1] = 0.0;
	n = 0;
	double bla = y0[0];
	x = 0;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx, k1, k2, k3, k4);
		

	if(yn[1]*y0[1] < 0){
	  n++;
	  if(n == 2)break;
	}
	        for(int i=0; i<dim; i++) y0[i] = yn[i];
		out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	
	}
	
	p0 = y0[1];
	pn = 0;
	thetar = 0;
	thetal = 1;
	a = 4;
	tol = 1e-10;
	cout << y0[1] << endl;
	while(a > tol){
	  
	  theta = (thetar+thetal)/2.0;
	  inter(k1,k2,k3,k4,p0,pn,theta,dx);
	  //cout << x << "\t" << theta << "\t" << pn << endl;
	  if(pn < 0)
	    thetal = theta;
	  else
	    thetar = theta;
	  a = abs(pn);
	}
	out.close();
	cout << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	double T = x+theta*dx; 
	bisec << bla << "\t" << T << endl;
	}
	bisec.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;
	

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
	double y[2] = { y0[0], y0[1]};

	y0[0] = y[1];
	y0[1] = -y[0]/sqrt(1+y[0]*y[0]);
}	

void inter(double* k1, double* k2,double* k3, double* k4, double p0, double& pn, double theta, double dx){
  double b1, b2, b3, b4;
  b1 = theta - 3*pow(theta,2)/2 + 2*pow(theta,3)/3;
  b2 = pow(theta,2)- 2*pow(theta,3)/3;
  b3 = b2;
  b4 = - pow(theta,2)/2 + 2*pow(theta,3)/3;
  pn = p0 + dx * (b1*k1[1] + b2*k2[1] + b3*k3[1] + b4*k4[1]);
}