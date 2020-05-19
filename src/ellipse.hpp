#ifndef ELLIPSE_HPP
#define ELLIPSE_HPP
#include"header.hpp"
#include"geometry.hpp"
class Ellipse:Geometry{
public:
    double theta; // orientation
    double a,b; // a --semi-major axis,b--semi-minor axis
    vec<double,2> center = {0.0,0.0};// center of this ellipse
    double Ma, Mp; //Ma-active torque(chirality) Mp-passive torque exerted on this object
    vec2 Fa,Fp;//active force(self-propulsion)  Fp--passive force exerted on this object 
    Ellipse() = default;
    Ellipse(double a,double b,double theta,vec<double,2>& center){
        this->a = a;
        this->b = b;
        this->center = center;
        this->theta = theta;
        this->tag = "ellipse";
    }
    void set_center(double x,double y){
        center[0] = x;
        center[1] = y;
    }
    
    void set_axes(double a,double b){
        this->a = a;
        this->b = b;
    }
    
    void set_angle(double t){theta = t;}

    coeff eqn_coeff();
    
    friend std::ostream& operator<<(std::ostream&,Ellipse&);
};

std::ostream& operator<<(std::ostream &os,Ellipse &e){
    std::vector<double> poly = e.polynomial();
    os<<"Ellipse "<< e.index<<std::endl;
    os<<"Orientational angle : "<<e.theta<<std::endl;
    os<<"Center : ("<<e.center[0]<<","<<e.center[1]<<")"<<std::endl;
    os<<"Semi-major axies : "<<e.a<<"  Semi-minor axies : "<<e.b<<std::endl;
    os<<"Equation : ["<<poly[0];
    for (size_t i = 1;i<6;i++){
        os<<","<<poly[i];
    }
    os<<"] * [x^2,y^2,xy,x,y,1]^T=0"<<std::endl;
    
    return os;
}

#endif