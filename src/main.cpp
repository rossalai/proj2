/* Project: 905 Project 2
 * File:   main.cpp
 * Author: alaina ross
 *
 * Created on February 9, 2017, 4:52 PM
 */

#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
//#include "../armadillo"
//#include </opt/local/include/armadillo.h>
#include<time.h>

#include "proj2.h"

using namespace std;
using namespace arma;

int main(){
    int interact=0;
//    cout<<"interacting? (0=no 1=yes): ";
//    cin>>interact;
    int n;
    cout<<"n: ";
    cin >> n;
    cout<<endl;
    jacobi(n,interact);
    return 0;
}

//template functions
//void print_vals(mat , mat ,int ,double );

// performs jacobi algorithm
// to find eigenvalues/vectors
// FOR UNIT TESTING USE n = 3 AND pmax = 10
int jacobi(int n, int interact) {
    cout.precision(4);
    double conv=0.001,wr=0.01, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    double aip=0, aiq=0, vpi=0, vqi=0;
    double tau=0, t=0, s=0, c=0;//tan(theta), sin(theta), cos(theta)    
    int count=1;                //count of iterations
    int count_old=count-10;     //keep track of every 10th iteration    
    int p=n-1, q=n-2;           //off diag all same value to start
                                //pick last as first maximum
    clock_t start, end;
    
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    
    //initialize x values
    r(0)=h;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }
    
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j && interact==0){
                a(i,j)=2/(h*h)+r(i)*r(i);
                v(i,j)=1;
            }
            else if (i==j && interact==1){
                a(i,j)=2/(h*h)+wr*wr+1/(r(i)*r(i));
                v(i,j)=1;
            }
            else if (i==j+1 or i==j-1){
                a(i,j)=-1/(h*h);
            } 
        }
    }
    
    if(n<=10){
        cout<<"Before diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    double app=a(p,p);
    double aqq=a(q,q);
    double apq=a(p,q);
    
    start=clock();
    
    while(abs(apq)>conv){
        if(count>1){
            apq=0;
            find_max(a,p,q,apq,n);
        }
//    //while(count<=5){
//        //find maximum non-diag matrix elements
//        if(count>1){
//            apq=0;
//            for (int i=0;i<n;i++){
//                 for (int j=0;j<n;j++){
//                    if(i!=j && abs(a(i,j))>=abs(apq)){
//                        apq=a(i,j);
//                        p=i;
//                        q=j;
//                    }
//                 }
//            }
//        }
//        cout<<"p: "<<p<<" q: "<<q<<endl;

        //unit test for max(a(i,j))    
        if(n==3 && count==1 && pmax==10){  
            cout<<"Testing max non-diagonal matrix element"<<endl;
            cout<<"Results should be: max = -0.09 p,q = 2,1"<<endl;
            cout<<setprecision(2)<<"Results are: max = "<<apq<<" p,q = "
                    <<p<<","<<q<<endl;
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

        //calculate new matrix elements and vectors
        for(int i=0;i<n;i++){
            if(i!=p && i!=q){
                aip=a(i,p);
                aiq=a(i,q);
                a(i,p)=aip*c-aiq*s;
                a(p,i)=aip*c-aiq*s;
                a(i,q)=aiq*c+aip*s;
                a(q,i)=aiq*c+aip*s;
            }
            vpi=v(p,i);
            vqi=v(q,i);
            v(p,i)=c*vpi-s*vqi;
            v(q,i)=c*vqi+s*vpi;
        }
        a(p,p)=app*c*c-2*apq*c*s-aqq*s*s;
        a(q,q)=app*s*s+2*apq*c*s+aqq*c*c;
        a(p,q)=0;
        a(q,p)=0;
                
        
//        cout<<"tau = "<<tau<<" t = "<<t<<" s = "<<s<<" c = "<<c<<endl;
//        print_vals(a,v,n,convergence);
//        cout<<endl;
//        cout<<"p: "<<p<<" q: "<<q<<endl<<endl;
//        cout<<abs(apq)<<endl<<endl;
        
        count++;
    }
    
    end=clock();
    
    //if(n<=10 && n!=3){
    if(n<=10){
        cout<<"After diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    //unit test for eigenvalues  
    if(n==3 && pmax==10 && interact==0){
        cout<<"Testing eigenvalues"<<endl;
        cout<<"Results should be: 11.2909, 44.6241, 100.18"<<endl;
        cout<<fixed<<"Results are: "<<a(0,0)<<", "<<a(1,1)<<", "<<a(2,2)<<endl;
        if(abs(11.2909-a(0,0))<0.01 && abs(44.6241-a(1,1))<0.01 
                && abs(100.18-a(2,2))<0.01){
            cout<<"Test passed"<<endl<<endl;
        }
        else
            cout<<"Test failed"<<endl<<endl;
    }
    
    //unit test for eigenvectors    
    if(n==3 && pmax==10 && interact==0){
        cout<<"Testing eigenvectors"<<endl;
        cout<<"Results should be:"<<endl;
        cout<<"v0:  0.9999  0.0081 0.0000"<<endl;
        cout<<"v1: -0.0081  0.9999 0.0027"<<endl; 
        cout<<"v2:  0.0000 -0.0027 0.9999"<<endl;
        cout<<"Results are: "<<endl;
        for (int i=0;i<n;i++){
            cout<<"v"<<i<<": ";
            for (int j=0;j<n;j++){
                if(abs(v(i,j))>conv)
                    cout<<fixed<<v(i,j)<<" ";
                else cout<<"0.0000 ";
            }
            cout<<endl;
        }
        if(abs(0.999-v(0,0))<0.01 && abs(0.008-v(0,1))<0.01 && abs(0-v(0,02))<0.01
                && abs(-0.0081-v(1,0))<0.01 && abs(0.999-v(1,1))<0.01
                && abs(0.0027-v(1,2))<0.01 && abs(0-v(2,0))<0.01
                && abs(-0.0027-v(2,1))<0.01 && abs(0.999-v(2,2))<0.01){
            cout<<"Test passed"<<endl<<endl;
        }
        else
            cout<<"Test failed"<<endl<<endl;
    }
    
    //unit test for eigenvectors    
    if(n==3 && pmax==10 && interact==0){
        cout<<"Testing orthogonality"<<endl;
        cout<<"Results should be: v0*v1 = 0, v0*v0 = 1"<<endl;
        double dot1=v(0,0)*v(1,0)+v(0,1)*v(1,1)+v(0,2)*v(1,2);
        double dot2=v(0,0)*v(0,0)+v(0,1)*v(0,1)+v(0,2)*v(0,2);
        cout<<"Results are: v0*v1 = "<<dot1<<", v0*v0 = "<<dot2<<endl;
        if(abs(dot1)<0.01 && abs(dot2-1)<0.01)
            cout<<"Test passed"<<endl<<endl;
        else
            cout<<"Test failed"<<endl<<endl;
    }
    
    cout<<"Diagonalization took "<<count<<" iterations"<<endl;
    cout<<scientific<<"CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    
    //print eigenvectors to files
    for(int i=0;i<n;i++){
        string filename="eigen_";
        filename+=to_string(i);
        filename+=".txt";
        ofstream outfile(filename);
        if(abs(a(i,i))<15&& a(i,i)>0){
            outfile<<"# "<<a(i,i)<<endl<<endl;
            for(int j=0;j<n;j++){
                outfile<<r(j)<<"   "<<r(j)*v(i,j)<<endl;
            }
            outfile.close();
        }
    }

    for(int i=0;i<n;i++){
        if(abs(a(i,i))<15 && a(i,i)>0){
             cout<<i<<" "<<fixed<<a(i,i)<<endl;
        }
    }
    
    
    return 0;
}

//find maximum non-diag matrix elements
void find_max(mat a,int& p,int& q,double& apq,int n){
    for (int i=0;i<n;i++){
         for (int j=0;j<n;j++){
            if(i!=j && abs(a(i,j))>=abs(apq)){
                apq=a(i,j);
                p=i;
                q=j;
            }
         }
    }
}

void print_vals(mat A, mat v,int n,double conv){
    cout<<"A: ";
    for (int i=0;i<n;i++){
        if(i>0){
            cout<<"   ";
         }
        for (int j=0;j<n;j++){
            if(abs(A(i,j))>conv)
                cout<<fixed<<A(i,j)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }
    for (int i=0;i<n;i++){
        cout<<"v"<<i<<": ";
        for (int j=0;j<n;j++){
            if(abs(v(i,j))>conv)
                cout<<fixed<<v(i,j)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }  
}