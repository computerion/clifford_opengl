#include <ctime> /* used to seed random generator */
#include "e3ga.h"
#include <iostream>
#include <cmath>
#include <string>
#include <stdlib.h>

using namespace std;
using namespace e3ga;

void init_test();
void cleanup();
mv neg(mv vector);
int randIndex();
double randDouble(double fMin, double fMax);
void compare(mv expected, mv result);

void test_3_1(int iterations);
void test_3_2(int iterations);
void test_3_3(int iterations);
void test_K_scaling(int index, double theta, mv p, mv expected);
mv K_left_rotor(int index, double theta);
mv K_right_rotor(int index, double theta);
void test_E_rotation(mv vector, double theta, int i, int j, mv expected);
mv E_right_rotor(double theta, int i, int j);
mv E_left_rotor(double theta, int i, int j);
void test_F_shear(mv vector, double theta, int i, int j, mv expected);
mv F_right_shear(double theta, int i, int j);
mv F_left_shear(double theta, int i, int j);
void test_4_1(int iterations);
void test_4_1_Lemma_1(int i, int j);
void test_4_1_Lemma_2(int i);
mv test_4_1_left_rotor(int i);
mv test_4_1_right_rotor(int i);
void test_4_2(int iterations);
void test_4_2_Shear(int i, int j, double t, mv vector, mv expected);
mv test_4_2_left_rotor(int i, int j, double t);
mv test_4_2_right_rotor(int i, int j, double t);

void test_5_1(int iterations);
void test_5_1_Classical_Shear(int i, int j, double t, mv p,  mv expected);
void test_5_2(int iterations);
void test_5_2_Uniform_Scaling(int i, double t, mv p,  mv expected);
void test_5_3(int iterations);
void test_5_3_Translation(int i, double t, mv p,  mv expected);
void test_5_4(int iterations);
void test_5_4_Rotation(int i, int j, double t, mv p,  mv expected);
void test_5_5(int iterations);
void test_5_5_Reflection(int i, mv p,  mv expected);
void test_composite(void);

mv **e, **be, **w, **sw; //The e's, bar e's, w's, and w stars
double **coords; //Used for construction of the above arrays and for deallocation

int ORDER = 12;

int NUM_3_1 = 5;
int NUM_3_1_RAND_MIN = 3;
int NUM_3_1_RAND_MAX = 10;
int NUM_3_1_RAND_DOUBLE_MIN = 0;
int NUM_3_1_RAND_DOUBLE_MAX = 2;

int NUM_3_2 = 5;
int NUM_3_2_RAND_MIN = 0;
int NUM_3_2_RAND_MAX = 10;

int NUM_3_3 = 5;
int NUM_3_3_RAND_MIN = 2;
int NUM_3_3_RAND_MAX = 6;

int NUM_4_1 = 5;

int NUM_4_2 = 5;
int NUM_4_2_RAND_MIN = 4;
int NUM_4_2_RAND_MAX = 10;

int NUM_5_1 = 5;
int NUM_5_1_RAND_MIN = 4;
int NUM_5_1_RAND_MAX = 10;

int NUM_5_2 = 5;
int NUM_5_2_RAND_MIN = 3;
int NUM_5_2_RAND_MAX = 10;
int NUM_5_2_RAND_DOUBLE_MIN = 0;
int NUM_5_2_RAND_DOUBLE_MAX = 2;

int NUM_5_3 = 5;
int NUM_5_3_RAND_MIN = 3;
int NUM_5_3_RAND_MAX = 10;
int NUM_5_3_RAND_DOUBLE_MIN = 0;
int NUM_5_3_RAND_DOUBLE_MAX = 15;

int NUM_5_4 = 5;
int NUM_5_4_RAND_MIN = 3;
int NUM_5_4_RAND_MAX = 10;
int NUM_5_4_RAND_DOUBLE_MIN = 0;
int NUM_5_4_RAND_DOUBLE_MAX = 15;

int NUM_5_5 = 5;
int NUM_5_5_RAND_MIN = 3;
int NUM_5_5_RAND_MAX = 10;

int main(int argc, char *argv[]) {

	init_test();
	
	test_3_1(NUM_3_1);
	test_3_2(NUM_3_2);
	test_3_3(NUM_3_3);

	test_4_1(NUM_4_1);
	test_4_2(NUM_4_2);

	test_5_1(NUM_5_1);
	test_5_2(NUM_5_2);
	test_5_3(NUM_5_3);
	test_5_4(NUM_5_4);
	test_5_5(NUM_5_5);

	test_composite();
	
	cleanup();

	return 0;
}

void test_3_1(int iterations){
	if (iterations == 0) return;
	cout<<"3.1 Ki = ei.bei: non-uniform scaling in the wi direction:"<<endl;

	for(int i=0;i<iterations;i++){
		double theta = randDouble(NUM_3_1_RAND_MIN, NUM_3_1_RAND_MAX);
		double theta2 = randDouble(NUM_3_1_RAND_DOUBLE_MIN, NUM_3_1_RAND_DOUBLE_MAX);
		int i1 = randIndex();
		int i2 = randIndex();
		while(i2 == i1){
			i2 = randIndex();
		}
		mv expected = gp(mv(exp(2*theta)), *w[i1]);
		cout<<"Testing non-uniform scaling in wi direction with rotor exp(theta ei bei)";
		test_K_scaling(i1, theta, *w[i1], expected);

		expected = gp(mv(exp(-2*theta2)), *sw[i1]);
		cout<<"Testing non-uniform scaling in swi direction with rotor exp(theta ei bei)";
		test_K_scaling(i1, theta2, *sw[i1], expected);

		cout<<"Testing non-uniform scaling in wj, i!=j, direction with rotor exp(theta ei bei)";
		test_K_scaling(i1, theta, *w[i2], *w[i1]);
		cout<<"Testing non-uniform scaling in swj, i!=j, direction with rotor exp(theta ei bei)";
		test_K_scaling(i1, theta, *sw[i2], *sw[i1]);
	}	

	cout<<"3.1 Test Done\n"<<endl;
}

void test_K_scaling(int index, double theta, mv p, mv expected){
	mv A_inv = K_left_rotor(index, theta);
	mv A = K_right_rotor(index, theta);
	mv result = gp(A_inv, gp(p, A));
	compare(expected, result);
}

mv K_left_rotor(int i, double theta){
	mv K = gp(*e[i], *be[i]);
	return cosh(mv(theta),ORDER) - gp(sinh(mv(theta),ORDER), K);
}

mv K_right_rotor(int i, double theta){
	mv K = gp(*e[i], *be[i]);
	return cosh(mv(theta),ORDER) + gp(sinh(mv(theta),ORDER), K);
}

void test_3_2(int iterations){
	if (iterations == 0) return;
	cout<<"3.2 Ei = ei.ej - bei.bej: simple rotation in the wi.wj-plane:"<<endl;

	for(int i=0;i<iterations;i++){
		//Testing e_i, e_j, w_i, w_j
		int i1 = randIndex();
		int i2 = randIndex();
		while(i2 == i1){
			i2 = randIndex();
		}
		int k = randIndex();
		while(k==i1 || k==i2){
			k = randIndex();
		}
		mv expected;
		double ran = randDouble(NUM_3_2_RAND_MIN, NUM_3_2_RAND_MAX);
		mv psin(sin(2*ran));
		mv nsin(-sin(2*ran));
		mv pcos(cos(2*ran));

		expected = add(gp(pcos, *w[i1]), gp(psin, *w[i2]));
		cout<<"Testing wi with rotor exp(theta(eiej - beibej))";
		test_E_rotation(*w[i1], ran, i1, i2, expected);

		expected = add(gp(nsin, *w[i1]), gp(pcos, *w[i2]));
		cout<<"Testing wj with rotor exp(theta(eiej - beibej))";
		test_E_rotation(*w[i2], ran, i1, i2, expected);

		cout<<"Testing swi with rotor exp(theta(eiej - beibej))";
		expected = add(gp(pcos, *w[i1]), gp(psin, *w[i2]));
		test_E_rotation(*sw[i1], ran, i1, i2, expected);

		cout<<"Testing swj with rotor exp(theta(eiej - beibej))";
		expected = add(gp(nsin, *sw[i1]), gp(pcos, *sw[i2]));
		test_E_rotation(*sw[i2], ran, i1, i2, expected);

		cout<<"Testing wk. k != i,j with rotor exp(theta(eiej - beibej))";
		test_E_rotation(*w[k], ran, i1, i2, *w[k]);

		cout<<"Testing swk. k != i,j with rotor exp(theta(eiej - beibej))";
		test_E_rotation(*sw[k], ran, i1, i2, *sw[k]);
	}	

	cout<<"3.2 Test Done\n"<<endl;
}

void test_E_rotation(mv vector, double theta, int i, int j, mv expected){
	mv right = E_right_rotor(theta, i, j);
	mv left = E_left_rotor(theta, i ,j);

	mv result = gp(left, gp(vector, right));
	compare(expected, result);
}

mv E_right_rotor(double theta, int i, int j){
	mv eiej = gp(*e[i], *e[j]);
	mv beibej = gp(*be[i], *be[j]);
	mv A = add(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), eiej));
	mv A_inv = subtract(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), beibej));
	return gp(A, A_inv);
}

mv E_left_rotor(double theta, int i, int j){
	mv eiej = gp(*e[i], *e[j]);
	mv beibej = gp(*be[i], *be[j]);
	mv B = subtract(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), eiej));
	mv B_inv = add(cos(mv(theta),ORDER), gp(sin(mv(theta),ORDER), beibej));
	return gp(B, B_inv);
}

void test_3_3(int iterations){
	if (iterations == 0) return;
	cout<<"3.3 Fi = ei.bej - bei.ej: scissors shear in the wi.wj-plane:"<<endl;

	for(int i=0;i<iterations;i++){
		//Testing e_i, e_j, w_i, w_j
		int i1 = randIndex();
		int i2 = randIndex();
		while(i2 == i1){
			i2 = randIndex();
		}
		int k = randIndex();
		while(k==i1 || k==i2){
			k = randIndex();
		}
		mv expected;
		double ran = randDouble(NUM_3_3_RAND_MIN, NUM_3_3_RAND_MAX);
		mv psin(sinh(2*ran));
		mv pcos(cosh(2*ran));

		expected = add(gp(pcos, *w[i1]), gp(psin, *w[i2]));
		cout<<"Testing wi with rotor exp(theta(eibej - beiej))";
		test_F_shear(*w[i1], ran, i1, i2, expected);

		expected = add(gp(psin, *w[i1]), gp(pcos, *w[i2]));
		cout<<"Testing wj with rotor exp(theta(eibej - beiej))";
		test_F_shear(*w[i2], ran, i1, i2, expected);

		expected = subtract(gp(pcos, *sw[i1]), gp(psin, *sw[i2]));
		cout<<"Testing swi with rotor exp(theta(eibej - beiej))";
		test_F_shear(*sw[i1], ran, i1, i2, expected);

		expected = subtract(gp(pcos, *sw[i2]), gp(psin, *sw[i1]));
		cout<<"Testing swj with rotor exp(theta(eibej - beiej))";
		test_F_shear(*sw[i2], ran, i1, i2, expected);

		cout<<"Testing wk k!=i,j with rotor exp(theta(eibej - beiej))";
		test_F_shear(*w[k], ran, i1, i2, *w[k]);

		cout<<"Testing swk k!=i,j with rotor exp(theta(eibej - beiej))";
		test_F_shear(*sw[k], ran, i1, i2, *sw[k]);
	}	

	cout<<"3.3 Test Done\n"<<endl;
}

void test_F_shear(mv vector, double theta, int i, int j, mv expected){
	mv right = F_right_shear(theta, i, j);
	mv left = F_left_shear(theta, i, j);

	mv result = gp(left, gp(vector, right));
	compare(expected, result);
}

mv F_right_shear(double theta, int i, int j){
	mv eibej = gp(*e[i], *be[j]);
	mv beiej = gp(*be[i], *e[j]);
	mv A = add(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), eibej));
	mv A_inv = subtract(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), beiej));
	return gp(A,A_inv);
}

mv F_left_shear(double theta, int i, int j){
	mv eibej = gp(*e[i], *be[j]);
	mv beiej = gp(*be[i], *e[j]);
	mv B = subtract(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), eibej));
	mv B_inv = add(cosh(mv(theta),ORDER), gp(sinh(mv(theta),ORDER), beiej));
	return gp(B,B_inv);
}

void test_4_1(int iterations){
	if (iterations == 0) return;
	cout<<"4.1 Reflection with Ri = ei.bei"<<endl;

	for(int i=0;i<iterations;i++){
		//Testing e_i, e_j, w_i, w_j
		int i1 = randIndex();
		int i2 = randIndex();
		while(i2 == i1){
			i2 = randIndex();
		}
		cout<<"Testing wj and swj with Rotor eibei";
		test_4_1_Lemma_1(i1, i2);
		cout<<"Testing wi and swi with Rotor eibei";
		test_4_1_Lemma_2(i1);
	}	

	cout<<"4.1 Test Done\n"<<endl;
}

void test_4_1_Lemma_1(int i, int j){
	mv left = test_4_1_left_rotor(i);
	mv right = test_4_1_right_rotor(j);
	mv result = gp(left, gp(*w[j], right));
	compare(*w[j], result);

	left = test_4_1_left_rotor(i);
	right = test_4_1_right_rotor(i);
	result = gp(left, gp(*sw[j], right));
	compare(*sw[j], result);
}

void test_4_1_Lemma_2(int i){
	mv left = test_4_1_left_rotor(i);
	mv right = test_4_1_right_rotor(i);
	mv result = gp(left, gp(*w[i], right));
	compare(neg(*sw[i]), result);

	left = test_4_1_left_rotor(i);
	right = test_4_1_right_rotor(i);
	result = gp(left, gp(*sw[i], right));
	compare(neg(*sw[i]), result);
}

mv test_4_1_left_rotor(int i){
	return gp(neg(*be[i]), *e[i]);
}

mv test_4_1_right_rotor(int i){
	return gp(*e[i], *be[i]);
}

void test_4_2(int iterations){
	if (iterations == 0) return;
	cout<<"4.2 Standard Form of Shear swi.wj"<<endl;

	for(int q=0;q<iterations;q++){
		int i = randIndex();
		int j = randIndex();
		int k = randIndex();
		double t = randDouble(NUM_4_2_RAND_MIN, NUM_4_2_RAND_MAX);
		while(i == j){
			j = randIndex();
		}
		while(k == i || k == j){
			k = randIndex();
		}
		mv tvec = mv(t);
		mv expected = add(*w[i], gp(tvec, *w[j]));
		cout<<"Testing wi with rotor (1+t(swi wj))";
		test_4_2_Shear(i,j,t,*w[i], expected);
		cout<<"Testing wj with rotor (1+t(swi wj))";
		test_4_2_Shear(i,j,t,*w[j], *w[j]);
		cout<<"Testing wk with rotor (1+t(swi wj))";
		test_4_2_Shear(i,j,t,*w[k], *w[k]);
		cout<<"Testing swi with rotor (1+t(swi wj))";
		test_4_2_Shear(i,j,t,*sw[i], *sw[i]);
		expected = subtract(*sw[j], gp(tvec,*sw[i]));
		cout<<"Testing wj with rotor (1+t(swi wj))";
		test_4_2_Shear(i,j,t,*sw[j], expected);
		cout<<"Testing swk with rotor (1+t(swi wj))";
		test_4_2_Shear(i,j,t,*sw[k], *sw[k]);		
	}	

	cout<<"4.2 Test Done\n"<<endl;
}

void test_4_2_Shear(int i, int j, double t, mv vector, mv expected){
	mv left = test_4_2_left_rotor(i,j,t);
	mv right = test_4_2_right_rotor(i,j,t);
	mv result = gp(left, gp(vector, right));
	compare(expected, result);
}

mv test_4_2_left_rotor(int i, int j, double t){
	return 1 -  t * gp(*sw[i], *w[j]);
}

mv test_4_2_right_rotor(int i, int j, double t){
	return 1 + t * gp(*sw[i], *w[j]);
}

void test_5_1(int iterations){
	if (iterations == 0) return;
	cout<<"5.1 Classical Shear"<<endl;
	double c[4];

	for(int q=0;q<iterations;q++){
		int i = randIndex();
		int j = randIndex();
		c[0] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		c[1] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		c[2] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		c[3] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		double t = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		while(i == j){
			j = randIndex();
		}

		mv p = c[0] * *w[0] + c[1] * *w[1] + c[2] * *w[2] + c[3] * *w[3];
		cout<<"i="<<i<<" j="<<j<<" t="<<t<<" p="<<c[0]<<"w0 + "<<c[1]<<"w1 + "<<c[2]<<"w2 + "<<c[3]<<"w3"<<endl;
		cout<<"Testing shear with rotor (1+t(swi wj))";
		test_5_1_Classical_Shear(i, j, t, p,  p + c[i] * (t* *w[j]) );
	}	

	cout<<"5.1 Test Done\n"<<endl;
}

void test_5_1_Classical_Shear(int i, int j, double t, mv p,  mv expected){
	mv right = test_4_2_right_rotor(i,j,t);
	mv left = test_4_2_left_rotor(i,j,t);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_2(int iterations){
	if (iterations == 0) return;
	cout<<"5.2 Uniform and non-uniform scaling"<<endl;
	double c[4];

	for(int q=0;q<iterations;q++){
		int i = randIndex();
		c[0] = randDouble(NUM_5_2_RAND_MIN, NUM_5_2_RAND_MAX);
		c[1] = randDouble(NUM_5_2_RAND_MIN, NUM_5_2_RAND_MAX);
		c[2] = randDouble(NUM_5_2_RAND_MIN, NUM_5_2_RAND_MAX);
		c[3] = randDouble(NUM_5_2_RAND_MIN, NUM_5_2_RAND_MAX);
		double t = randDouble(NUM_5_2_RAND_DOUBLE_MIN, NUM_5_2_RAND_DOUBLE_MAX);

		mv p = c[0] * *w[0] + c[1] * *w[1] + c[2] * *w[2] + c[3] * *w[3];
		cout<<"i="<<i<<" t="<<t<<" p="<<c[0]<<"w0 + "<<c[1]<<"w1 + "<<c[2]<<"w2 + "<<c[3]<<"w3"<<endl;
		cout<<"Testing scaling with rotor exp(theta ei bei)";
		mv expected = p - c[i] * *w[i];
		expected = expected + c[i] * exp(2*t) * *w[i];
		test_5_2_Uniform_Scaling(i, t, p, expected);
	}	

	cout<<"5.2 Test Done\n"<<endl;
}

void test_5_2_Uniform_Scaling(int i, double t, mv p,  mv expected){
	mv right = K_right_rotor(i,t);
	mv left = K_left_rotor(i,t);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_3(int iterations){
	if (iterations == 0) return;
	cout<<"5.3 Translation"<<endl;
	double c[4];

	for(int q=0;q<iterations;q++){
		int i = randIndex();
		while(i==0){
			i = randIndex();
		}
		c[0] = randDouble(NUM_5_3_RAND_MIN, NUM_5_3_RAND_MAX);
		c[1] = randDouble(NUM_5_3_RAND_MIN, NUM_5_3_RAND_MAX);
		c[2] = randDouble(NUM_5_3_RAND_MIN, NUM_5_3_RAND_MAX);
		c[3] = randDouble(NUM_5_3_RAND_MIN, NUM_5_3_RAND_MAX);
		double t = randDouble(NUM_5_3_RAND_DOUBLE_MIN, NUM_5_3_RAND_DOUBLE_MAX);

		mv p = c[0] * *w[0] + c[1] * *w[1] + c[2] * *w[2] + c[3] * *w[3];
		cout<<"i="<<i<<" t="<<t<<" p="<<c[0]<<"w0 + "<<c[1]<<"w1 + "<<c[2]<<"w2 + "<<c[3]<<"w3"<<endl;
		cout<<"Testing translation with rotor 1+ (t sw0 wi)";
		mv expected = p + c[0] * t * *w[i];
		test_5_3_Translation(i, t, p, expected);
	}	

	cout<<"5.3 Test Done\n"<<endl;
}

void test_5_3_Translation(int i, double t, mv p,  mv expected){
	mv right = test_4_2_right_rotor(0,i,t);
	mv left = test_4_2_left_rotor(0,i,t);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_4(int iterations){
	if (iterations == 0) return;
	cout<<"5.4 Rotation"<<endl;
	double c[4];

	for(int q=0;q<iterations;q++){
		int i = randIndex();
		int j = randIndex();
		c[0] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		c[1] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		c[2] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		c[3] = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		double t = randDouble(NUM_5_1_RAND_MIN, NUM_5_1_RAND_MAX);
		while(i == j){
			j = randIndex();
		}

		mv p = c[0] * *w[0] + c[1] * *w[1] + c[2] * *w[2] + c[3] * *w[3];
		cout<<"i="<<i<<" j="<<j<<" t="<<t<<" p="<<c[0]<<"w0 + "<<c[1]<<"w1 + "<<c[2]<<"w2 + "<<c[3]<<"w3"<<endl;
		cout<<"Testing rotation with rotor exp( theta ( eiej - beibej))";
		mv expected = p - c[i] * *w[i] - c[j] * *w[j];
		expected += c[i] * (cos(2*t) * *w[i] + sin(2*t) * *w[j]);
		expected += c[j] * (-sin(2*t) * *w[i] + cos(2*t) * *w[j]);
		test_5_4_Rotation(i, j, t, p,  expected);
	}	

	cout<<"5.4 Test Done\n"<<endl;
}

void test_5_4_Rotation(int i, int j, double t, mv p,  mv expected){
	mv right = E_right_rotor(t, i, j);
	mv left = E_left_rotor(t, i, j);;
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_5_5(int iterations){
	if (iterations == 0) return;
	cout<<"5.5 Rotation"<<endl;
	double c[4];

	for(int q=0;q<iterations;q++){
		int i = randIndex();
		c[0] = randDouble(NUM_5_5_RAND_MIN, NUM_5_5_RAND_MAX);
		c[1] = randDouble(NUM_5_5_RAND_MIN, NUM_5_5_RAND_MAX);
		c[2] = randDouble(NUM_5_5_RAND_MIN, NUM_5_5_RAND_MAX);
		c[3] = randDouble(NUM_5_5_RAND_MIN, NUM_5_5_RAND_MAX);

		mv p = c[0] * *w[0] + c[1] * *w[1] + c[2] * *w[2] + c[3] * *w[3];
		cout<<"i="<<i<<" p="<<c[0]<<"w0 + "<<c[1]<<"w1 + "<<c[2]<<"w2 + "<<c[3]<<"w3"<<endl;
		cout<<"Testing reflection with rotor ei bei";
		mv expected = p - 2 * c[i] * *w[i];
		test_5_5_Reflection(i, p,  expected);
	}	

	cout<<"5.5 Test Done\n"<<endl;
}

void test_5_5_Reflection(int i, mv p,  mv expected){
	mv right = test_4_1_right_rotor(i);
	mv left = test_4_1_left_rotor(i);
	mv result = gp(left, gp(p, right));
	compare(expected, result);
}

void test_composite(void) {
	cout<<"Composite Test"<<endl;
    mv shearleft;
    mv shearright;
    mv rotateleft;
    mv rotateright;
    double theta=M_PI/3.;
    mv compositeleft;
    mv compositeright;

    shearleft = test_4_2_left_rotor(1,2,2.);
    shearright = test_4_2_right_rotor(1,2,2.);
    rotateright = E_right_rotor(theta, 1, 3);
    rotateleft = E_left_rotor(theta, 1,3);

    compositeleft = gp(shearleft,rotateleft);
    compositeright = gp(rotateright,shearright);

    cout << "shearleft = " << shearleft.toString() << endl;
    cout << "shearright = " << shearright.toString() << endl;

    cout << "compositeleft = " << compositeleft.toString() << endl;
    cout << "compositeright = " << compositeright.toString() << endl;

    mv r = gp(compositeleft,gp(*w[1],compositeright));
    cout << "composite on w1 = " << r.toString() << endl;

    r = gp(compositeleft,gp(*w[2],compositeright));
    cout << "composite on w2 = " << r.toString() << endl;

    r = gp(compositeleft,gp(*w[3],compositeright));
    cout << "composite on w3 = " << r.toString() << endl;
    cout<<"Test Done\n"<<endl;
}

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
