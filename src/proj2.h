/* 
 * File:   proj2.h
 * Author: alaina
 *
 * Created on February 21, 2017, 2:04 PM
 */

#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include </opt/local/include/armadillo.h>
#include<time.h>

#ifndef PROJ2_H
#define	PROJ2_H
using namespace std;
using namespace arma;

void print_vals(mat,mat,int,double);
int jacobi(int, int);
void find_max(mat,int&,int&,double&,int);

#endif /* PROJ2_H */

