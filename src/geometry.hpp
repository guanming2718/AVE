#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include"header.hpp"
class Geometry{ //a template class
public:
    unsigned index;
    std::string tag;
//member functions
    void set_index(unsigned i) {index = i;}
    unsigned get_index(){return index;}                           
    std::vector<double> eqn_coeff();    
    std::vector<vec2 > intersection(Geometry&)=0; // intersection points of two geometry object
    vec2 tangent_vector(vec2&)=0;
    friend std::ostream& operator<<(std::ostream&,Geometry&);
};
#endif
