#include <functional>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
//performs Romberg integration over f from a to b
//computes until the error is less than tolerance or until R_n,n is found, whichever comes first
double romberg(std::function<double(double)> f, double a, double b, int n, double tolerance)
{
  double h = b-a;
  //we only keep two rows of the extrapolation table in memory at a time
  std::vector<double> currRow, lastRow;
  lastRow.push_back(0.5*h*(f(a)+f(b))); //R_1,1
  for(int i = 2; i <= n; i++)
  {
    currRow.clear();
    double sum = 0;
    for(int k = 1; k <= pow(2,i-2); k++)
      sum += f(a+(k-0.5)*h);
    currRow.push_back(0.5*(lastRow[0] + h*sum));
    for(int j = 1; j < i; j++)
      currRow.push_back(currRow[j-1] + (currRow[j-1]-lastRow[j-1])/(pow(4,j)-1));
    h *= 0.5;
    if(fabs(currRow.back() - lastRow.back()) < tolerance)
    {
      std::cout << "finished with n = " << i << std::endl;
      return currRow.back();
    }
    lastRow = currRow;
  }
  return currRow.back();
}

int main()
{
  std::cout.precision(10);
  double tolerance = 1e-6;

  //part c
  std::function<double(double)> fc = [](double x){
    return 2/(x*x-4);
  };
  std::cout << "5c:\n" << romberg(fc, 0, 0.35, 10, tolerance) << std::endl;

  //part e
  std::function<double(double)> fe = [](double x){
    return exp(3*x)*sin(2*x);
  };
  std::cout << "5e:\n" << romberg(fe, 0, M_PI/4, 10, tolerance) << std::endl;
  return 0;
}
