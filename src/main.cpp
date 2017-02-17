/* Project: 905 Project 2
 * File:   main.cpp
 * Author: alaina ross
 *
 * Created on February 9, 2017, 4:52 PM
 */

#include <iostream>
#include<fstream>
#include<math.h>
//#include "../armadillo"
#include </opt/local/include/armadillo.h>
#include<time.h>
using namespace std;
using namespace arma;

//template functions
void print_vals(mat , vec , vec ,int);

// performs jacobi algorithm
// to find eigenvalues/vectors
int main(int argc, char** argv) {
    //FOR UNIT TESTING USE n = 3 AND pmax = 10
    int n;
    cout<<"n: ";
    cin >> n;
    
    double convergence=0.01;
    double pmin=0;
    double pmax=10;      
    double h = (pmax-pmin)/(double(n));
    clock_t start,end;
    
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
                a(i,j)=2/(h*h)+r(i);
            else if (i==j+1 or i==j-1){
                a(i,j)=-1/(h*h);
            }
            
        }
        //b(i)=h*h*100*exp(-10*(x(i)));
    }
    
    if(n<=10){
        cout<<"Before diagonalization"<<endl;
        print_vals(a,b,v,n);
        cout<<endl;
    }
    
    //off diag all same value to start
    //pick last as first maximum
    int p=n-1;
    int q=n-2;
    double app=a(p,p);
    double aqq=a(q,q);
    double apq=a(p,q);
    double aip=0;
    double aiq=0;
    int count=1;        //count of iterations
    double tau=0;
    double t=0;         //tan(theta)
    double s=0;         //sin(theta)
    double c=0;         //cos(theta)

    
    start=clock();
    
    
    while(abs(apq)>convergence){
    //while(count<=5){
        //find maximum non-diag matrix elements
        if(count>1){
            apq=0;
            for (int i=0;i<n;i++){
                 for (int j=0;j<n;j++){
                    if(i!=j && abs(a(i,j))>=abs(apq)){
                        apq=a(i,j);
                        //cout<<apq<<endl;
                        p=i;
                        q=j;
                    }
                 }
            }
        }
        
//        cout<<"p: "<<p<<" q: "<<q<<endl;

        //unit test for max(a(i,j))    
        if(n==3 && count==1 && pmax==10){  
            cout<<endl<<"Testing max non-diagonal matrix element"<<endl;
            cout<<"Results should be: max = -0.09 p,q = 2,1"<<endl;
            cout<<"Results are: max = "<<apq<<" p,q = "<<p<<","<<q<<endl;
            if(abs(-0.09-apq)<0.001 && p==2 && q==1)
                cout<<"Test passed"<<endl<<endl;
            else
                cout<<"Test failed"<<endl<<endl;
        }

        //calculate sin(theta) and cos(theta)
        aqq=a(q,q);
        app=a(p,p);
        tau=(aqq-app)/(2*apq);
        if(tau>0)
            t=1/(tau+sqrt(1+tau*tau));
        else
            t=1/(-tau+sqrt(1+tau*tau));   
        c=1/sqrt(1+t*t);
        s=c*t;

        //calculate new matrix elements
        for(int i=0;i<n;i++){
            if(i!=p && i!=q){
                aip=a(i,p);
                aiq=a(i,q);
                a(i,p)=aip*c-aiq*s;
                a(i,q)=aiq*c+aip*s;
            }
        }
        a(p,p)=app*c*c-2*apq*c*s-aqq*s*s;
        a(q,q)=app*s*s+2*apq*c*s+aqq*c*c;
        a(p,q)=0;
        
//        cout<<"tau = "<<tau<<" t = "<<t<<" s = "<<s<<" c = "<<c<<endl;
//        print_vals(a,b,v,n);
//        cout<<endl;
//        cout<<"p: "<<p<<" q: "<<q<<endl<<endl;
//        cout<<abs(apq)<<endl<<endl;
        count++;
    }
    
    end=clock();
    
    if(n<=10){
        cout<<"After diagonalization"<<endl;
        print_vals(a,b,v,n);
        cout<<endl;
    }

    //unit test for eigenvalues)    
    if(n==3 && pmax==10){
        cout<<endl<<"Testing eigenvalues"<<endl;
        cout<<"Results should be: 0.1776, 3.513, 6.848"<<endl;
        cout<<"Results are: "<<a(0,0)<<", "<<a(1,1)<<", "<<a(2,2)<<endl;
        if(abs(0.1776-a(0,0))<0.01 && abs(3.513-a(1,1))<0.01 
                && abs(6.848-a(2,2))<0.015){
            cout<<"Test passed"<<endl<<endl;
        }
        else
            cout<<"Test failed"<<endl<<endl;
    }
    
    cout<<"Diagonalization took "<<count<<" iterations"<<endl;
    cout<<"CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    
    
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
//    cout<<"b: ";
//    for (int i=0;i<n;i++){
//        cout<<b(i)<<" ";
//    }
//    cout<<endl;
//    cout<<"v: ";
//    for (int i=0;i<n;i++){
//        cout<<v(i)<<" ";
//    }
//    cout<<endl;
    
}