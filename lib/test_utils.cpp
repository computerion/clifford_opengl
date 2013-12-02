#include <ctime> /* used to seed random generator */
#include "e3ga.h"
#include "test_utils.h"
#include <iostream>
#include <cmath>
#include <string>
#include <stdlib.h>

using namespace std;
using namespace e3ga;

mv **e, **be, **w, **sw; //The e's, bar e's, w's, and w stars
double **coords; //Used for construction of the above arrays and for deallocation

extern int ORDER;

void init_test(){
	srand(time(NULL));
	e = new mv*[4];
	be = new mv*[4];
	w = new mv*[4];
	sw = new mv*[4];
	coords = new double*[16];

	for(int i=0;i<4;i++){
		double *coord1 = new double[8];
		double *coord2 = new double[8];
		double *coord3 = new double[8];
		double *coord4 = new double[8];
		for(int j=0;j<8;j++){
			coord1[j] = 0;
			coord2[j] = 0;
			coord3[j] = 0;
			coord4[j] = 0;
		}
		coord1[i] = 1;
		coord1[4+i] = 1;
		coord2[i] = 1;
		coord2[4+i] = -1;
		coord3[i] = 1;
		coord4[4+i] = 1;

		coords[i] = coord1;
		coords[4+i] = coord2;
		coords[8+i] = coord3;
		coords[12+i] = coord4;

		e[i] = new mv(GRADE_1, coord1);
		be[i] = new mv(GRADE_1, coord2);
		w[i] = new mv(GRADE_1, coord3);
		sw[i] = new mv(GRADE_1, coord4);
	}

	cout<<"Beginning Test\n"<<endl;
}

void cleanup(){
	for(int i=0;i<4;i++){
		delete e[i];
		delete be[i];
		delete w[i];
		delete sw[i];
	}

	for(int i=0;i<16;i++){
		delete [] coords[i];
	}

	delete [] e;
	delete [] be;
	delete [] coords;
	delete [] w;
	delete [] sw;
	cout<<"Done"<<endl;
}

mv neg(mv vector){
	mv zero = mv(0);
	return subtract(zero, vector);
}

int randIndex() {
	return rand()%4;
}

double randDouble(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void compare(mv expected, mv result){
	mv diff = subtract(expected, result);
	double largeCoord = diff.largestCoordinate();
	if(largeCoord >= 5){
		cout<<"-- Calculation-Off\nExpected: "<<expected.toString()<<"\nGot: "<<result.toString()<<"\nDifference: "<<diff.toString()<<endl<<endl;
	}
	else{
		cout<<"--- Success"<<endl;
	}
}

mv scale_Left(int i, double theta){
	mv K = gp(*e[i], *be[i]);
	return cosh(mv(theta),ORDER) - gp(sinh(mv(theta),ORDER), K);
}

mv scale_right(int i, double theta){
	mv K = gp(*e[i], *be[i]);
	return cosh(mv(theta),ORDER) + gp(sinh(mv(theta),ORDER), K);
}

mv rotate_left(double theta, int i, int j){
	mv eiej = gp(*e[i], *e[j]);
	mv beibej = gp(*be[i], *be[j]);
	mv B = subtract(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), eiej));
	mv B_inv = add(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), beibej));
	return gp(B, B_inv);
}

mv rotate_right(double theta, int i, int j){
	mv eiej = gp(*e[i], *e[j]);
	mv beibej = gp(*be[i], *be[j]);
	mv A = add(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), eiej));
	mv A_inv = subtract(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), beibej));
	return gp(A, A_inv);
}

mv scshear_left(double theta, int i, int j){
	mv eibej = gp(*e[i], *be[j]);
	mv beiej = gp(*be[i], *e[j]);
	mv B = subtract(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), eibej));
	mv B_inv = add(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), beiej));
	return gp(B,B_inv);
}

mv scshear_right(double theta, int i, int j){
	mv eibej = gp(*e[i], *be[j]);
	mv beiej = gp(*be[i], *e[j]);
	mv A = add(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), eibej));
	mv A_inv = subtract(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), beiej));
	return gp(A,A_inv);
}

mv reflect_left(int i){
	return gp(neg(*be[i]), *e[i]);
}

mv reflect_right(int i){
	return gp(*e[i], *be[i]);
}

mv shear_left(int i, int j, double t){
	return 1 -  t * gp(*sw[i], *w[j]);
}

mv shear_right(int i, int j, double t){
	return 1 + t * gp(*sw[i], *w[j]);
}

void test_3_1_Scaling(int index, double theta, mv p, mv expected){
	mv A_inv = scale_Left(index, theta);
	mv A = scale_right(index, theta);
	mv result = gp(A_inv, gp(p, A));
	compare(expected, result);
}

void test_3_2_Rotation(mv vector, double theta, int i, int j, mv expected){
	mv right = rotate_right(theta, i, j);
	mv left = rotate_left(theta, i ,j);

	mv result = gp(left, gp(vector, right));
	compare(expected, result);
}

void test_3_3_Shear(mv vector, double theta, int i, int j, mv expected){
	mv right = scshear_right(theta, i, j);
	mv left = scshear_left(theta, i, j);

	mv result = gp(left, gp(vector, right));
	compare(expected, result);
}

void test_4_1_Reflection_Lemma_1(int i, int j){
	mv left = reflect_left(i);
	mv right = reflect_right(j);
	mv result = gp(left, gp(*w[j], right));
	compare(*w[j], result);

	left = reflect_left(i);
	right = reflect_right(i);
	result = gp(left, gp(*sw[j], right));
	compare(*sw[j], result);
}

void test_4_1_Reflection_Lemma_2(int i){
	mv left = reflect_left(i);
	mv right = reflect_right(i);
	mv result = gp(left, gp(*w[i], right));
	compare(neg(*sw[i]), result);

	left = reflect_left(i);
	right = reflect_right(i);
	result = gp(left, gp(*sw[i], right));
	compare(neg(*sw[i]), result);
}

void test_4_2_Shear(int i, int j, double t, mv vector, mv expected){
	mv left = shear_left(i,j,t);
	mv right = shear_right(i,j,t);
	mv result = gp(left, gp(vector, right));
	compare(expected, result);
}

void test_5_1_Classical_Shear(int i, int j, double t, mv p,  mv expected){
	mv left = shear_left(i,j,t);
	mv right = shear_right(i,j,t);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_2_Uniform_Scaling(int i, double t, mv p,  mv expected){
	mv right = scale_right(i,t);
	mv left = scale_Left(i,t);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_3_Translation(int i, double t, mv p,  mv expected){
	mv left = shear_left(0,i,t);
	mv right = shear_right(0,i,t);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_4_Rotation(int i, int j, double t, mv p,  mv expected){
	mv right = rotate_right(t, i, j);
	mv left = rotate_left(t, i, j);;
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_5_Reflection(int i, mv p,  mv expected){
	mv right = reflect_left(i);
	mv left = reflect_right(i);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}