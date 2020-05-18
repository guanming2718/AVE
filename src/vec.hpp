#ifndef VEC_HPP
#define VEC_HPP
#include<vector>
#include<stdexcept>
#include<array>
#include<iostream>
using std::vector;
template<class T,size_t D> // T for type, L is the length of the vector
class vec{
    public:
        std::array<T,D> arr;    
        vec() = default;
        vec& operator+=(const vec& v){ 
            //vec1 += vec2 equals to vec1.operator+=(vec2)
            for (size_t i=0;i<D;i++){
                arr[i] += v.arr[i];
            }
            return *this;
        }

        vec& operator-=(const vec& v){ 
            //vec1 += vec2 equals to vec1.operator+=(vec2)
            for (size_t i=0;i<D;i++){
                arr[i] -= v.arr[i];
            }
            return *this;
        }
        
        vec operator+(const vec& v) const{
           vec ret;
           for(size_t i = 0; i<D; ++i)
             ret[i] = arr[i] + v.arr[i];
           return ret;
        }

        vec operator-(const vec& v) const{
           vec ret;
           for(size_t i = 0; i<D; ++i)
             ret[i] = arr[i] - v.arr[i];
           return ret;
        }
        vec& operator=(vec& v){
            for (size_t i=0;i<D;i++){
                arr[i] = v.arr[i];
            }
            return *this;
        }
        
        size_t size(){
            return D;
        }
        
        T& operator[](size_t i){
            if (i < 0 || i>=D){
                throw std::out_of_range("vec index out of range");
            }else{
                return arr[i];
            }
        }

        T poly_value(T x){
          double pow_x = 1.0;
          T sum = 0.0;
          for (size_t i=0;i<D;i++){
              // pow_x = x^i
              sum += arr[i]*pow_x;
              pow_x *= x;
          }
          return sum;
        }

        T norm(){
           T sum_sqr = 0.0;
           for (size_t i=0;i<D;i++){
            sum_sqr += arr[i]*arr[i];
           }
           return std::sqrt(sum_sqr);
        }
        
        template<class T0,size_t D0> 
        friend std::ostream& operator<<(std::ostream&,const vec<T0,D0>& );
        
};
//end of the class 

// external functions
template<class T,size_t D> 
std::ostream& operator<<(std::ostream& os,const vec<T,D> &v){ 
    //cout<< vec equals to  operator<<(cout,vec)
    // print [arr[0],arr[1],arr[2] ...]
    os << '[';
    if (D > 0){
        os << v.arr[0];
    }
    for (size_t i =1;i<D;i++){
        os << ',' << v.arr[i];
    }
    os << ']';
    return os;
}
#endif
