#include "catch.hpp"
#include "proj2.h"


TEST_CASE("Testing max a(i,j)"){
    int n=3;
    double pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    mat a = zeros<mat>(n,n);
    vec r(n);
    r(0)=h;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j){
                a(i,j)=2/(h*h)+r(i)*r(i);
            }
            else if (i==j+1 or i==j-1){
                a(i,j)=-1/(h*h);
            } 
        }
    }
    int p=0;
    int q=0;
    double apq=0;
    find_max(a,p,q,apq,n);
    REQUIRE(p==2);
    REQUIRE(q==1);
    REQUIRE(apq==Approx(-0.09));
    
}

TEST_CASE("Testing eigenvalues"){
    
}
TEST_CASE("Testing eigenvector orthogonality"){
    
}
