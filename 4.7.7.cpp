#include <functional>
#include <iostream>
#include <cmath>

//integrates f from a to b using Gaussian Quadrature with n=5
double gaussQuad(std::function<double(double)> f, double a, double b)
{
    const double roots[] = {
      0.9061798459,
      0.5384693101,
      0.0000000000,
      -0.5384693101,
      -0.9061798459
    };
    const double coeff[] = {
      0.2369268850,
      0.4786286705,
      0.5688888889,
      0.4786286705,
      0.2369268850
    };

    //change of variable from x to t so we can integrate from -1 to 1
    std::function<double(double)> g = [&](double t){
      return f(((b-a)*t+b+a)/2)*(b-a)/2;
    };

    double sum = 0;
    for(int i = 0; i < 5; i++)
      sum += coeff[i]*g(roots[i]);
    return sum;
}

int main()
{
  //part a
  double a = 1;
  double b = 1.5;
  std::function<double(double)> f = [](double x){
    return x*x*log(x);
  };
  double integral = gaussQuad(f, a, b);
  std::cout.precision(10);
  std::cout << "integral from " << a << " to " << b << " of (x^2)lnx = " << integral << std::endl;

  //part d
  a = 0;
  b = M_PI/4;
  f = [](double x){
    return x*x*sin(x);
  };
  integral = gaussQuad(f, a, b);
  std::cout << "integral from " << a << " to " << b << " of (x^2)sinx = " << integral << std::endl;
  return 0;
}
