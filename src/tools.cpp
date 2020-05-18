#include "tools.hpp"

// do not use key word static when implementing a static function outside the header

double Tools::eval_coeff(coeff& p,double x,double y){
    return p[0]*x*x + p[1]*y*y + p[2]*x*y + p[3]*x + p[4]*y +p[5];
}
std::vector<double> Tools::roots_filter(complex_array& roots){
    std::vector<double> real_roots;
    // filter out complex roots
    for (auto const& root:roots){
        if (Tools::is_real(root)){
            real_roots.push_back(root.real());
        }
    }
    return real_roots;

}
std::vector<vec2> Tools::union_eqs(coeff &p0,coeff &p1){
    //  union of 4th order polynomial(Bezout's Theorem) 
    //  A_1*x^2 + B_1*y^2 + C_1*xy + D_1*x + E_1*y + F =0 and 
    //  A_0*x^2 + B_0*y^2 + C_0*xy + D_0*x + E_0*y + F =0
    //  with p0 = (A_0,B_0,C_0,D_0,E_0)
    //  check https://www.math.unm.edu/~vageli/papers/FLEX/IntersectionOfEllipses.pdf for detail
    //  convert to matrix form [x,y][a00,a01; a01,a11][x,y]^T +[b0,b1][x,y]^T + c = 0
    double a00_0,a01_0,a11_0,b0_0,b1_0,c_0; //matrix elements for ellipse number 0 
    double a00_1,a01_1,a11_1,b0_1,b1_1,c_1; //matrix elements for ellipse number 1
    a00_0 = p0[0]; //a00 = A
    a00_1 = p1[0];
    a11_0 = p0[1]; // a11 = B
    a11_1 = p1[1]; 
    a01_0 = 0.5*p0[2]; // a01 = 0.5*C
    a01_1 = 0.5*p1[2]; 
    b0_0 = p0[3]; // b0 = D
    b0_1 = p1[3];
    b1_0 = p0[4]; // b1 = E
    b1_1 = p1[4]; 
    c_0 = p0[5]; // c=F
    c_1 = p1[5]; // c=F
    //aux variavble
    double v[11];
    v[0] = 2.0*(a00_0*a01_1 - a00_1*a01_0);
    v[1] = a00_0*a11_1 - a00_1*a11_0;
    v[2] = a00_0*b0_1 - a00_1*b0_0;
    v[3] = a00_0*b1_1 - a00_1*b1_0;
    v[4] = a00_0*c_1 - a00_1*c_0;
    v[5] = 2.0*(a01_0*a11_1 - a01_1*a11_0);
    v[6] = 2.0*(a01_0*b1_1 - a01_1*b1_0);
    v[7] = 2.0*(a01_0*c_1 - a01_1*c_0);
    v[8] = a11_0*b0_1 - a11_1*b0_0;
    v[9] = b0_0*b1_1 - b0_1*b1_0;
    v[10]= b0_0*c_1 - b0_1*c_0;
    // polynomial R(y) =u0 + u1*y + u2*y^2+ u3*y^3 + u4*y^4
    double u[5];
    u[0] = v[2]*v[10] - v[4]*v[4];
    u[1] = v[0]*v[10] + v[2]*(v[7] + v[9]) - 2.0*v[3]*v[4];
    u[2] = v[0]*(v[7] + v[9]) + v[2]*(v[6]- v[8]) - v[3]*v[3] - 2.0*v[1]*v[4];
    u[3] = v[0]*(v[6] - v[8]) + v[2]*v[5] - 2.0*v[1]*v[3];
    u[4] = v[0]*v[5] - v[1]*v[1];

    std::vector<vec2> points;
    
    //get complex roots
    complex_array y_complex= Tools::solve_quartic(u[4],u[3],u[2],u[1],u[0]);
    std::vector<double> y = Tools::roots_filter(y_complex);
    
    //  solve x from A_0*x^2 + (C_0*y + D0)*x + (B_0*y^2+ E_0*y + F) =0 
    for (auto const& y0:y){
        double a = p0[0],b = p0[2]*y0 + p0[3],c = p0[1]*y0*y0 + p0[4]*y0 + p0[5];
        complex_array x_complex = Tools::solve_quadratic(a,b,c);
        std::vector<double> x_test = Tools::roots_filter(x_complex);
        for (auto const& x0:x_test){
            //  test A_0*x^2 + B_0*y^2 + C_0*xy + D_0*x + E_0*y + F =0
            if (abs(Tools::eval_coeff(p1,x0,y0)) < 1e-3){
                points.push_back({x0,y0});  
            }
        }  
    }
    //return intersection points 
    return points;
   
}

complex_array Tools::solve_quartic(double a,double b,double c,double d,double e){
    //solve equation ax^4 + bx^3 + cx^2 +dx +e = 0
    // https://en.wikipedia.org/wiki/Quartic_function  section Solution methods
    complex_array roots;
    if (a==0.0){
        return Tools::solve_cubic(b,c,d,e);
    }
    if ((b==0.0) && (c==0.0) &&(d==0.0) &&(e==0.0)){
        roots.push_back(0.0);
        return roots;
    }
    // convert to standard form x^4 + alpha3*x^3 + alpha2*x^2 +alpha1*x + alpha0 = 0
    double alpha[4] = {e/a,d/a,c/a,b/a};
    // convert to depressed equation 
    // y^4 + py^2 +qy +r = 0 with x = y-0.25*a3
    // the value of p,q,r are from https://mathworld.wolfram.com/QuarticEquation.html
    double p = alpha[2] - 3.0/8.0*alpha[3]*alpha[3];
    double q = alpha[1] - 0.5*alpha[2]*alpha[3] + 0.125*pow(alpha[3],3.0);
    double r = alpha[0] - 0.25*alpha[1]*alpha[3] + 1.0/16.0*alpha[2]*alpha[3]*alpha[3] - 3.0/256.0*pow(alpha[3],4.0);
    // solve the cubic equation 8m^3 + 8p*m^2 + (2p^2-8r)m - q^2 = 0
    complex_array cubic_roots = Tools::solve_cubic(8.0,8.0*p,2.0*p*p-8.0*r,-q*q);
    // chose a nonzero root
    complex<double> m,sqt_m;
    for(auto const& root: cubic_roots) {
        if (abs(root) > 1e-8){
            m = root;
            break;
        }
    }
    sqt_m = sqrt(m);
    roots.push_back(-0.25*alpha[3] + 0.5*( sqrt(2.0)*sqt_m + sqrt(-(2.0*p+2.0*m + sqrt(2.0)*q/sqt_m))));
    roots.push_back(-0.25*alpha[3] + 0.5*( sqrt(2.0)*sqt_m - sqrt(-(2.0*p+2.0*m + sqrt(2.0)*q/sqt_m))));
    roots.push_back(-0.25*alpha[3] - 0.5*( sqrt(2.0)*sqt_m + sqrt(-(2.0*p+2.0*m - sqrt(2.0)*q/sqt_m))));
    roots.push_back(-0.25*alpha[3] - 0.5*( sqrt(2.0)*sqt_m - sqrt(-(2.0*p+2.0*m - sqrt(2.0)*q/sqt_m))));
    return roots;
}

complex_array Tools::solve_cubic(double a,double b,double c,double d){
    // solve ax^3 + bx^2 + cx + d = 0
    // general cubic formula https://en.wikipedia.org/wiki/Cubic_equation 
    complex_array roots;
    if (a==0.0){
        return Tools::solve_quadratic(b,c,d);
    }
    if ((b==0.0) && (c==0.0) &&(d==0.0)){
        roots.push_back(0.0);
        return roots;
    }
    complex<double> delta0 = b*b-3.0*a*c;
    complex<double> delta1 = 2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d;
    complex<double> sqt = sqrt(pow(delta1,2.0) - 4.0*pow(delta0,3.0));
    complex<double> C; 
    if (abs(delta1 + sqt) > 1e-8){
        C = pow(0.5*(delta1 + sqt),1.0/3.0);
    } 
    else{
        // change sign if delta2 if C==0
        C = pow(0.5*(delta1 - sqt),1.0/3.0);
    }
    if (abs(C) < 1e-8){
        roots.push_back(-1.0/(3.0*a)*b);
        return roots;
    }
    complex<double> zeta(-0.5,0.5*std::sqrt(3.0)),zeta2 = zeta*zeta;
    roots.push_back(-1.0/(3.0*a)*(b + C + delta0/C) );
    roots.push_back(-1.0/(3.0*a)*(b + zeta*C + delta0/(zeta*C)) );
    roots.push_back(-1.0/(3.0*a)*(b + zeta2*C + delta0/(zeta2*C)) );

    return roots;
}

complex_array Tools::solve_quadratic(double a,double b,double c){
    // solve quadratice equation ax^2 + b^x + c = 0
    std::vector<complex<double> > roots;
    if (a==0.0){
        roots.push_back(-c/b);
        return roots;
    }
    // convert coefficients to complex number
    complex<double> delta = b*b - 4.0*a*c;
    roots.push_back(0.5/a*(-b + sqrt(delta)));
    roots.push_back(0.5/a*(-b - sqrt(delta)));
    return roots;
}

