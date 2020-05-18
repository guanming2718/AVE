#ifndef TOOLS_HPP
#define TOOLS_HPP
#include "header.hpp"
using complex_array = std::vector<complex<double> >;
class Tools{
  public:
    static complex_array solve_quartic(double a0,double a1,double a2,double a3,double a4);
    static complex_array solve_cubic(double a0,double a1,double a2,double a3);
    static complex_array solve_quadratic(double a0,double a1,double a2);
    static std::vector<vec2> union_eqs(coeff &p1,coeff &p2);  
    static std::vector<double> roots_filter(complex_array& roots);
    
    static double eval_coeff(coeff& p,double x,double y);
    
    static bool is_real(complex<double> num){
      if (std::abs(num) < 1e-6){
        // number = 0
        return true;
      }
      else if(num.imag() == 0.0){
        return true;
      }
      else if(num.real() == 0.0){
        return false;
      }
      else if (num.imag()/num.real() < 1e-3){
        return true;
      } 
    return false;
   }
};

#endif