#include "geometry.hpp"
#include "tools.hpp"

coeff Ellipse::eqn_coeff(){
        // return the array for the coefficients Ax^2 + By^2 + Cxy + Dx + Ey + F = 0
        // transform in to global xy coordinates
        double s = sin(theta), c = cos(theta);
        double T1 = a*a*s*s + b*b*c*c;
        double T2 = 2.0*(b*b*s*c-a*a*s*c);
        double T3 = b*b*s*s + a*a*c*c;
        double A = T1;
        double B = T3;
        double C = T2;
        double D = -2.0*T1*center[0] - T2*center[1];
        double E = -2.0*T3*center[1] - T2*center[0];
        double F = T1*center[0]*center[0] + T2*center[0]*center[1] + T3*center[1]*center[1] - a*a*b*b; 
        coeff poly = {A,B,C,D,E,F};
        return poly;
}

std::vector<vec2 > Ellipse::intersection(Geometry& geo){
    coeff c1=this->eqn_coeff();
    coeff c2= geo.eqn_coeff();
    std::vector<vec2> points = Tools::union_eqs(c1,c2);
    // get the box packing this ellipse
    double x_max = center[0] + a*cos(theta),x_min = center[1] - a*cos(theta);
    double y_max = center[0] + a*sin(theta),y_min = center[1] - a*sin(theta);
    if (x_min>x_max) std::swap(x_min,x_max);
    if (y_min>y_max) std::swap(y_min,y_max);
    std::vector<bool> added ={false}; 
    std::vector<vec2> checked_points;
    // eliminates points that are too close
    for (unsigned i=0;i<points.size();i++){
        if (added[i]) continue;
        checked_points.push_back(points[i]);
        added[i] = true;
        for (unsigned j=i+1;j<points.size();j++){
            if (added[j]) continue;
            // distance between two points
            double d = (points[i] - points[j]).norm(); 
            if (d<0.1*b) added[j] = true;
        }
    }   
    return checked_points;
}
 