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
  int maxN = 10000;

  //part a
  std::function<double(double)> integrand1 = [](double x){
    return 1/(1+x*x*x*x);
  };
  std::function<double(double)> integrand2 = [](double t){
    return (t*t)/(1+t*t*t*t);
  };

  double integral1 = romberg(integrand1, 0, 1, maxN, tolerance);
  double integral2 = romberg(integrand2, 0, 1, maxN, tolerance);

  std::cout << "First integral = " << integral1 << std::endl;
  std::cout << "Second integral = " << integral2 << std::endl;
  std::cout << "Sum = " << integral1 + integral2 << std::endl;


  return 0;
}
