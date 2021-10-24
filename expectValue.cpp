#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>

//integrates f from 0 to infinity
//uses Gauss-Laguerre Quadrature with numPoints
double gaussLaguerre(std::function<double(double)> f, int numPoints)
{
  std::vector<double> x, w;
  switch (numPoints)
  {
    case 3:
      x = {
        4.15774556783479E-01,
        2.29428036027904E+00,
        6.28994508293748E+00
      };
      w = {
        7.11093009929173E-01,
        2.78517733569241E-01,
        1.03892565015861E-02
      };
      break;
    case 7:
      x = {
        1.93043676560362E-01,
        1.02666489533919E+00,
        2.56787674495075E+00,
        4.90035308452648E+00,
        8.18215344456286E+00,
        1.27341802917978E+01,
        1.93957278622625E+01
      };
      w = {
        4.09318951701273E-01,
        4.21831277861720E-01,
        1.47126348657505E-01,
        2.06335144687169E-02,
        1.07401014328075E-03,
        1.58654643485641E-05,
        3.17031547899558E-08
      };
      break;
    case 15:
      x = {
        9.33078120172818E-02,
        4.92691740301884E-01,
        1.21559541207095E+00,
        2.26994952620374E+00,
        3.66762272175144E+00,
        5.42533662741355E+00,
        7.56591622661307E+00,
        1.01202285680191E+01,
        1.31302824821757E+01,
        1.66544077083300E+01,
        2.07764788994488E+01,
        2.56238942267288E+01,
        3.14075191697539E+01,
        3.85306833064860E+01,
        4.80260855726858E+01
      };
      w = {
        2.18234885940085E-01,
        3.42210177922883E-01,
        2.63027577941712E-01,
        1.26425818105934E-01,
        4.02068649210009E-02,
        8.56387780361183E-03,
        1.21243614721425E-03,
        1.11674392344252E-04,
        6.45992676202291E-06,
        2.22631690709627E-07,
        4.22743038497937E-09,
        3.92189726704109E-11,
        1.45651526407312E-13,
        1.48302705111330E-16,
        1.60059490621113E-20
      };

      break;
  }

  double sum = 0.0;
  for(int i = 0; i < numPoints; i++)
      sum +=  w[i] * f(x[i]) * exp(x[i]);
  return sum;
}

int main()
{
  std::function<double(double)> R = [](double x){
    double Z = 1.0;
    double a0 = 0.529e-10;
    double r = x*a0;
    double a = 1.0 / (4.0 * sqrt(2.0*M_PI))*pow(Z/a0, 1.5);
    double b = 2.0 - Z/a0*r;
    double c = exp(-(Z/a0*r)/2.0);
    double psi = a*b*c;
    return psi*psi*r*r*r;
  };

  std::cout.precision(15);
  double integral = gaussLaguerre(R, 3);
  std::cout << "3-Point Gauss-Laguerre" << std::endl;
  std::cout << "\tValue of integral is " << integral << std::endl;
  std::cout << "\tValue of expectation value is " << (integral * 4 * M_PI) << std::endl;
  integral = gaussLaguerre(R, 7);
  std::cout << "7-Point Gauss-Laguerre" << std::endl;
  std::cout << "\tValue of integral is " << integral << std::endl;
  std::cout << "\tValue of expectation value is " << (integral * 4 * M_PI) << std::endl;
  integral = gaussLaguerre(R, 15);
  std::cout << "15-Point Gauss-Laguerre" << std::endl;
  std::cout << "\tValue of integral is " << integral << std::endl;
  std::cout << "\tValue of expectation value is " << (integral * 4 * M_PI) << std::endl;
}
