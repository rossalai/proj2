/* Project: 905 Project 2
 * File:   main.cpp
 * Author: alaina ross
 *
 * Created on February 9, 2017, 4:52 PM
 */

#include <iostream>
#include<fstream>
//#include "../armadillo"
#include </opt/local/include/armadillo.h>
#include<time.h>
using namespace std;
using namespace arma;

//template functions
void print_vals(mat , vec , vec ,int);

// performs gauss elimination
// for general matrix
int main(int argc, char** argv) {
    int n;
    cout<<"n: ";
    cin >> n;
    double pmin=0;
    double pmax=10;      

    double h = (pmax-pmin)/(double(n));
    clock_t begin,finish;
    
    mat a = zeros<mat>(n,n);
    vec b = zeros<vec>(n);
    vec v = zeros<vec>(n);
    vec r(n);
    
    //initialize x values
    r(0)=0;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }

    
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j)
                a(i,j)=(2/(h*h))+r(i)*r(i);
            else if (i==j+1 or i==j-1){
                a(i,j)=-1/(h*h);
                a(i,j)=1+j-i*i;
            }
            
        }
        //b(i)=h*h*100*exp(-10*(x(i)));
    }
    
    if(n<=10){
        //cout<<"Before LU decomp"<<endl;
        print_vals(a,b,v,n);
    }
    
    double max=0;
    int p=0;
    int q=0;
    begin=clock();
    
    //find maximum non-diag matrix elements
    for (int i=0;i<n;i++){
         for (int j=0;j<n;j++){
            if(i!=j && abs(a(i,j))>max){
                max=a(i,j);
                p=i;
                q=j;
            }
         }
    }
  
//    Unit test for max(a(i,j))    
    cout<<endl<<"Testing max non-diagonal matrix element"<<endl;
    cout<<"Results should be: max = -6 p,q = 3,2"<<endl;
    cout<<"max: "<<max<<endl;
    cout<<"p,q: "<<p<<","<<q<<endl<<endl;
    
    
    finish=clock();
    
    cout<<"Total CPU time (sec) : "<<((double)finish-(double)begin)/
            CLOCKS_PER_SEC<<endl;
    
    
//    vec u(n);
//    ofstream outfile("output_general.txt");
//    ofstream error("error_general.txt");
//    ofstream actual("output_analytic.txt");
//    outfile<<"0.0 0.0"<<endl;
//    actual<<"0.0 0.0"<<endl;
//    for (int i=0;i<n;i++){
//        outfile<<x(i)<<" "<<v(i)<<endl;
//        u(i)=1-(1-exp(-10))*x(i)-exp(-10*x(i));
//        actual<<x(i)<<" "<<u(i)<<endl;
//        error<<x(i)<<" "<<log10(abs((v(i)-u(i))/u(i)))<<endl;
//    }
//    outfile<<"1.0 0.0"<<endl;
//    actual<<"1.0 0.0"<<endl;
//    outfile.close();
//    actual.close();
//    error.close();
    
    
    return 0;
}

void print_vals(mat A, vec b, vec v,int n){
    cout<<"A: ";
    for (int i=0;i<n;i++){
        if(i>0){
            cout<<"   ";
         }
        for (int j=0;j<n;j++){
            cout<<A(i,j)<<" ";
        }
        cout<<endl;
    }
    cout<<"b: ";
    for (int i=0;i<n;i++){
        cout<<b(i)<<" ";
    }
    cout<<endl;
    cout<<"v: ";
    for (int i=0;i<n;i++){
        cout<<v(i)<<" ";
    }
    cout<<endl;
    
}