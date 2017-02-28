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
    cout<<endl;
//    cout<<"interacting? (0=no 1=yes): ";
//    cin>>interact;
    vec r(n);
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    double conv=0.001,wr=0.01, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    
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
        for(int j=0;j<n;j++){
            outfile<<r(j)<<"   "<<r(j)*eigenvec(i,j)<<endl;
        }
        outfile.close();
    }
    
   /*for(int i=0;i<n;i++){
       // cout<<i<<" "<<fixed<<a(i,i)<<endl;
        cout<<i<<" "<<eigen[i]<<endl;
    }*/

    return 0;
}

