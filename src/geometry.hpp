#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include"header.hpp"
class Geometry{ //an abstract class
public:
    unsigned index;
    std::string tag;
//member functions
    void set_index(unsigned i) {index = i;}
    unsigned get_index(){return index;}                           
    coeff eqn_coeff();    
    // pure virtual functions
    virtual std::vector<vec2 > intersection(Geometry&)=0; // intersection points of two geometry object
    virtual vec2 tangent_vector(vec2&)=0;
    friend std::ostream& operator<<(std::ostream&,Geometry&);
};

class Ellipse:Geometry{
public:
    //geometry parameter
    double theta; // orientation
    double a,b; // a --semi-major axis,b--semi-minor axis
    double c; // distance from the center to the focal points
    vec2 center = {0.0,0.0};// center of this ellipse
    vec2 f1,f2;//two focal points    

    
    Ellipse() = default;
    Ellipse(double a,double b,double theta,vec2& center){
        this->a = a;
        this->b = b;
        this->center = center;
        this->theta = theta;
        this->tag = "ellipse";
        c = sqrt(a*a-b*b);
        vec2 t = {cos(theta),sin(theta)};
        f1[0] = center[0] - c*t[0];
        f1[1] = center[0] - c*t[0];
        f2[0] = center[0] + c*t[0];
        f2[1] = center[0] + c*t[0];
    }
    void set_center(double x,double y){
        center[0] = x;
        center[1] = y;
    }
    
    void set_axes(double a,double b){
        this->a = a;
        this->b = b;
    }
    std::vector<vec2 > intersection(Geometry&);

    vec2 tangent_vector(vec2& p){    //this point p should be on this object
       vec<double,6> c = this->eqn_coeff();
       vec2 t;
       // tx = 2By+Cx+E,ty=-(2Ax+Cy+D)
       t[0] = 2.0*c[1]*p[1] + c[2]*p[0] + c[4];
       t[1] = -(2.0*c[0]+c[2]*p[1]+c[3]);
       return t; 
    }

    void set_angle(double t){theta = t;}
    coeff eqn_coeff();
    
    bool inside_ellipse(vec2 p){
        double d1 = (p - f1).norm();
        double d2 = (p - f2).norm();
        if (d1+d2>2.0*c) return false;
        return true;
    }
    
    friend std::ostream& operator<<(std::ostream&,Ellipse&);
};

std::ostream& operator<<(std::ostream &os,Ellipse &e){
    coeff poly = e.eqn_coeff();
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
