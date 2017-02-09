/* Project: 905 Project 1
 * File:   main_general.cpp
 * Author: alaina ross
 *
 * Created on January 12, 2017, 4:16 PM
 */

#include <iostream>
#include<fstream>
#include <armadillo>
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

    double h = 1/(double(n)+1);
    clock_t begin,finish;
    
    mat a = zeros<mat>(n,n);
    vec b(n);
    vec v = zeros<vec>(n);
    vec x(n);
    
    //initialize x values
    x(0)=h;
    for (int i=1; i<n ;i++)
        x(i)=x(i-1)+h;

    
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j)
                a(i,j)=2;
            else if (i==j+1 or i==j-1)
                a(i,j)=-1;
        }
        b(i)=h*h*100*exp(-10*(x(i)));
    }
    
    if(n<=10){
        cout<<"Before Gauss elimination"<<endl;
        print_vals(a,b,v,n);
    }
    
    begin=clock();
    
      //LU decomposition    
//    mat L=trimatl(a);
//    mat U=trimatu(a);
//    vec y=solve(trimatl(L),b);
//    v=solve(trimatu(U),y);
    
    
    //do forward elimination
    double temp=0;
    for (int k=0;k<=n-2;k++){
        for (int i=k+1;i<=n-1;i++){
            temp=a(i,k)/a(k,k);
            for ( int j=k;j<=n-1;j++){
                a(i,j)-=a(k,j)*temp;
            }
            b(i)-=temp*b(k);
            if(abs(b(i))<=pow(10,-14))
                b(i)=0;
        }
    }

    //do back substitution
    if(abs(a(n-1,n-1))<=pow(10,-14)){
        v(n-1)=1;
        a(n-1,n-1)=0;
    }
    else{
        v(n-1)=b(n-1)/a(n-1,n-1);
    }
    double sum=0;
    for (int i=n-2;i>=0;i--){
        for (int j=n-1;j>i;j--){
            sum+=a(i,j)*v(j);
        }
        v(i)=(b(i)-sum)/a(i,i);
        if(abs(v(i))<=pow(10,-14))
            v(i)=0;
        sum=0;
    }
    
    finish=clock();
    
    if(n<=10){
        cout<<"After Gauss elimination"<<endl;
        print_vals(a,b,v,n); 
    }
    
    cout<<"Total CPU time (sec) : "<<((double)finish-(double)begin)/
            CLOCKS_PER_SEC<<endl;
    
    
    vec u(n);
    ofstream outfile("output_general.txt");
    ofstream error("error_general.txt");
    ofstream actual("output_analytic.txt");
    outfile<<"0.0 0.0"<<endl;
    actual<<"0.0 0.0"<<endl;
    for (int i=0;i<n;i++){
        outfile<<x(i)<<" "<<v(i)<<endl;
        u(i)=1-(1-exp(-10))*x(i)-exp(-10*x(i));
        actual<<x(i)<<" "<<u(i)<<endl;
        error<<x(i)<<" "<<log10(abs((v(i)-u(i))/u(i)))<<endl;
    }
    outfile<<"1.0 0.0"<<endl;
    actual<<"1.0 0.0"<<endl;
    outfile.close();
    actual.close();
    error.close();
    
    
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