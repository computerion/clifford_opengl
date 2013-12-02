#ifndef TEST_UTILS
#define TEST_UTILS

#include "e3ga.h"
#include <stdlib.h>

using namespace e3ga;

extern mv **e, **be, **w, **sw; //The e's, bar e's, w's, and w stars
extern double **coords; //Used for construction of the above arrays and for deallocation

void init_test(void);
void cleanup(void);
mv neg(mv vector);
int randIndex();
double randDouble(double fMin, double fMax);

mv scale_Left(int index, double theta);
mv scale_right(int index, double theta);

mv rotate_left(double theta, int i, int j);
mv rotate_right(double theta, int i, int j);

mv scshear_left(double theta, int i, int j);
mv scshear_right(double theta, int i, int j);

mv reflect_left(int i);
mv reflect_right(int i);

mv shear_left(int i, int j, double t);
mv shear_right(int i, int j, double t);

void test_3_1_Scaling(int index, double theta, mv p, mv expected);
void test_3_2_Rotation(mv vector, double theta, int i, int j, mv expected);
void test_3_3_Shear(mv vector, double theta, int i, int j, mv expected);
void test_4_1_Reflection_Lemma_1(int i, int j);
void test_4_1_Reflection_Lemma_2(int i);
void test_4_2_Shear(int i, int j, double t, mv vector, mv expected);
void test_5_1_Classical_Shear(int i, int j, double t, mv p,  mv expected);
void test_5_2_Uniform_Scaling(int i, double t, mv p,  mv expected);
void test_5_3_Translation(int i, double t, mv p,  mv expected);
void test_5_4_Rotation(int i, int j, double t, mv p,  mv expected);
void test_5_5_Reflection(int i, mv p,  mv expected);

#endif