#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>

//integrates from a to b until error is less than tolerance
//performs adaptive quadrature with simpson's rule
//this version is slower, but simpler
//because it recalculates values that are used multiple times
double adaptiveQuadSlow(std::function<double(double)> & f, double a, double b, double tolerance)
{
  //integrates f from x0 to x1 using simpson's rule
  std::function<double(double, double)> simpson = [&](double x0, double x1){
    double h = (x1-x0)/2;
    return (h/3)*(f(x0) + 4*f(x0+h) +f(x1));
  };

  double mid = (a+b)/2;
  double whole = simpson(a, b);
  double left = simpson(a, mid);
  double right = simpson(mid, b);
  double diff = left+right-whole;
  if(fabs(diff) < 10*tolerance) //10*tolerance is used because that's what is used in the Burden text
    return left + right;

  return adaptiveQuadSlow(f, a, mid, tolerance/2) + adaptiveQuadSlow(f, mid, b, tolerance/2);
}

//integrates from a to b until error is less than tolerance
//performs adaptive quadrature with simpson's rule
double adaptiveQuad(std::function<double(double)> & f, double a, double b, double tolerance)
{
  //integrates f from x0 to x1 using simpson's rule
  //f_x0 is f evaluated at x0
  //f_x1 is f  evaluated at x1
  auto simpson = [&](double x0, double x1, double f_x0, double f_x1){
    double mid = (x0 + x1)/2;
    double f_mid = f(mid);
    double result =  (fabs(x1-x0)/6)*(f_x0 + 4*f_mid + f_x1);
    return std::make_tuple(mid, f_mid, result);
  };

  //performs adaptive quadrature recursively
  //function evaluations and intermediate integrals are passed through arguments to avoid recalcuating values
  std::function<double (double, double, double, double, double, double, double, double)> helper =
  [&](double x0, double x1, double f_x0, double f_x1, double tol, double whole, double mid, double f_mid){
    double left_mid, f_left_mid, left, right_mid, f_right_mid, right;
    std::tie(left_mid, f_left_mid, left) = simpson(x0, mid, f_x0, f_mid);
    std::tie(right_mid, f_right_mid, right) = simpson(mid, x1, f_mid, f_x1);
    double diff = left + right - whole;
    if (fabs(diff) < 10*tol)
      return left + right;
    return helper(x0, mid, f_x0, f_mid, tol/2, left, left_mid, f_left_mid)
         + helper(mid, x1, f_mid, f_x1, tol/2, right, right_mid, f_right_mid);
  };

  double f_a = f(a);
  double f_b = f(b);
  double mid, f_mid, whole;
  std::tie(mid, f_mid, whole) = simpson(a, b, f_a, f_b);
  return helper(a, b, f_a, f_b, tolerance, whole, mid, f_mid);
}

int main()
{
  double tolerance = 1e-4;
  tolerance = 1e-5;
  std::function<double(double)> integrandC = [](double w){
    return cos(M_PI/2*w*w);
  };
  std::function<double(double)> integrandS = [](double w){
    return sin(M_PI/2*w*w);
  };

  std::cout << "t\tc(t)\t\ts(t)\n";
  for(double t = 0.1; t <= 1.0; t+=0.1)
  {
    std::cout << t << '\t'
              << std::setw(8) << adaptiveQuad(integrandC, 0, t, tolerance) << '\t'
              << std::setw(8) << adaptiveQuad(integrandS, 0, t, tolerance) << std::endl;
  }
  return 0;
}
