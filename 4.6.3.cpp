#include <functional>
#include <iostream>
#include <cmath>
#include <unordered_map>
#include <string>
#include <sstream>

double adaptiveQuad(std::function<double(double)> & f, double a, double b, double tolerance, std::unordered_map<std::string, double> & lookup)
{
  //simpson's method with memoization for repeated calls
  std::function<double(double, double)> simpson = [&](double x0, double x1){
    std::stringstream ss;
    ss << std::hexfloat << x0 << '|' << x1;
    std::string key = ss.str();
    if(lookup.count(key))
      return lookup[key];
    double h = (x1-x0)/2;
    double result = (h/3)*(f(x0) + 4*f(x0+h) +f(x1));
    lookup[key] = result;
    return result;
  };

  double mid = (a+b)/2;
  double whole = simpson(a, b);
  double left = simpson(a, mid);
  double right = simpson(mid, b);
  double diff = left+right-whole;
  if(fabs(diff) < 10*tolerance) //10*tolerance is used because that's what is used in the Burden text
    return left + right;

  return adaptiveQuad(f, a, mid, tolerance/2, lookup) + adaptiveQuad(f, mid, b, tolerance/2, lookup);
}

//integrates f from a to b using adaptive quadrature until tolerance is reached
double adaptiveQuad(std::function<double(double)> & f, double a, double b, double tolerance)
{
  //create a memoization table
  std::unordered_map<std::string, double> lookup;
  return adaptiveQuad(f, a, b, tolerance, lookup);
}


int main()
{
  //part a
  double a = 1;
  double b = 1.5;
  double tolerance = 1e-3;
  std::function<double(double)> f = [](double x){
    return x*x*log(x);
  };
  double integral = adaptiveQuad(f, a, b, tolerance);
  std::cout.precision(10);
  std::cout << "integral from " << a << " to " << b << " of (x^2)lnx = " << integral << std::endl;
  return 0;
}
