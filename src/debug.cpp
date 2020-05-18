#include"debug.hpp"
void test_quadratic(){
    // two complex roots
    complex_array roots = Tools::solve_quadratic(1.0,2.0,3.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;
    // two real roots
    roots = Tools::solve_quadratic(2.0,10.0,1.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;
    // two coinciding roots
    roots = Tools::solve_quadratic(1.0,2.0,1.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;

}

void test_cubic(){
    // three real roots
    complex_array roots;
    roots = Tools::solve_cubic(2.0,-30.0,136.0,-168.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;
    
    // three coinciding roots
    roots = Tools::solve_cubic(1.0,-3.0,3.0,-1.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;
    // two complex roots + one real root
    roots = Tools::solve_cubic(10.0,34.0,47.0,90.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;

    // two coinciding roots + one real root
    roots = Tools::solve_cubic(1.0,-9.0,24.0,-20.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;
}

void test_quatic(){
    // four different roots
    complex_array roots;
    roots = Tools::solve_quartic(1.0,-10.0,35.0,-50.0,24.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;

    // four roots
    roots = Tools::solve_quartic(10.0,9.0,8.0,7.0,6.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;

    // two double roots
    roots = Tools::solve_quartic(1.0,-12.0,52.0,-96.0,64.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;

    // third order roots + 1 root
    roots = Tools::solve_quartic(1.0,-15.0,81.0,-189.0,162.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;

    // third order roots + 1 root
    roots = Tools::solve_quartic(-36,0.0,-37.6,0,-9.0);
    for (auto root:roots){
        std::cout<<root<<' ';
    }
    std::cout<<std::endl;

    
}

void test_union(){
    coeff p1 = {1.0,1.0,0.0,0.0,0.0,-1.0};// x^2 + y^2 = 1 
    coeff p2 = {0.,0.,0.,1.0,1.0,0.0};// x+y = 0
    vector<vec2> pts = Tools::union_eqs(p1,p2);
    for (auto pt :pts){
        std::cout<<'('<<pt<<") ";
    }
    std::cout<<std::endl<<std::endl;
    

    coeff p3 = {1.0,1.0,0.0,0.0,0.0,-1.0};// x^2 + y^2 = 5
    coeff p4 = {1.0,1.0,0.0,-2.0,0.0,0.0};// (x-1)^2 + y^2 = 0
    pts = Tools::union_eqs(p3,p4);
    for (auto pt :pts){
        std::cout<<'('<<pt<<") ";
    }
    std::cout<<std::endl<<std::endl;
   
    coeff p5 = {1.0,1.0,0.0,0.0,0.0,-5.0};// x^2 + y^2 = 5
    coeff p6 = {0.2,1.0,0.0,3.0,0.0,0.0};// 0.2x^2 + y^2 +3x = 1
    pts = Tools::union_eqs(p5,p6);
    for (auto pt :pts){
        std::cout<<'('<<pt<<") ";
    }
    std::cout<<std::endl<<std::endl;

    coeff p7 = {2.0,1.0,0.0,4.0,0.0,-5.0};//  2x^2 +y ^2 +4x= 5 
    coeff p8 = {0.2,1.0,0.0,3.0,0.0,0.0};// 0.2x^2 + y^2 +3x = 1
    pts = Tools::union_eqs(p7,p8);
    for (auto pt :pts){
        std::cout<<'('<<pt<<") ";
    }
    std::cout<<std::endl<<std::endl;

}