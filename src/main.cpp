/* Project: 905 Project 2
 * File:   main.cpp
 * Author: alaina ross
 *
 * Created on February 9, 2017, 4:52 PM
 */

#include <cstdlib>
#include </opt/local/include/armadillo.h>

#include "jacobi.h"

using namespace std;
using namespace arma;

/*
 * 
 */
int main(int argc, char** argv) {
    int n;
    int interact=0;
    cout<<"n: ";
    cin >> n;
    cout<<"interacting? (0=no 1=yes): ";
    cin>>interact;
    cout<<endl;
    vec r(n);
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    double conv=0.0001,wr=5.0, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    n=n-1;
    initialize(n,h,a,r,v,interact,wr);
    jacobi(n,interact,conv,wr,a,r,v);
    vector<double>eigenval=get_eigenvals(a,n);
    mat eigenvec(3,n);
    eigenvec=get_eigenvecs(a,v,n);
    //print eigenvectors to files
    for(int i=0;i<3;i++){
        string filename="eigen_";
        filename+=to_string(i);
        filename+=".txt";
        ofstream outfile(filename);
        outfile<<"# "<<eigenval[i]<<endl<<endl;
        outfile<<"0   0"<<endl;
        for(int j=0;j<n;j++){
            outfile<<r(j)<<"   "<<eigenvec(i,j)<<endl;
        }
        outfile<<(r(n-1)+h)<<"   "<<0<<endl;
        outfile.close();
    }
    
    cout<<endl;
    for(int i=0;i<3;i++){
        cout<<fixed<<i<<" "<<eigenval[i]<<endl;
    }
    
/*    //uncomment to use armadillo
    clock_t start, end;
    initialize(n,h,a,r,v,interact,wr);
    vec eigval;
    mat eigvec;
    
    start=clock();
    eig_sym(eigval,eigvec,a);
    end=clock();
    
    cout<<scientific<<"CPU TIME (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    */

    return 0;
}

