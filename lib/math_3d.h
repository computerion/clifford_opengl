#ifndef MATH_3D_H
#define	MATH_3D_H

#include <stdio.h>
#ifdef WIN32
#define _USE_MATH_DEFINES 
#include <cmath>
#else
#include <math.h>
#endif
#include "e3ga.h"
#include <iostream>
using namespace std;

using namespace e3ga;
extern mv **e, **be, **w, **sw; //The e's, bar e's, w's, and w stars
extern double **coords; //Used for construction of the above arrays and for deallocation

#define ToRadian(x) (float)(((x) * M_PI / 180.0f))
#define ToDegree(x) (float)(((x) * 180.0f / M_PI))

struct Vector3f
{
    float x;
    float y;
    float z;

    Vector3f()
    {
    }

    Vector3f(float _x, float _y, float _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    Vector3f& operator+=(const Vector3f& r)
    {
        x += r.x;
        y += r.y;
        z += r.z;

        return *this;
    }

    mv to_mv(){
        return x * *w[1] + y * *w[2] + z * *w[3];
    }

    void applyRotors(mv left, mv right){
        mv vectorForm = to_mv();
        mv product = gp(left, gp(vectorForm, right));
        cout<<x<<" "<<y<<" "<<z<<endl;
        const double *coord = product.getC();
        x = (float) coord[1];
        y = (float) coord[2];
        z = (float) coord[3];
        cout<<product.toString()<<endl<<vectorForm.toString()<<endl;
    }

    Vector3f& operator-=(const Vector3f& r)
    {
        x -= r.x;
        y -= r.y;
        z -= r.z;

        return *this;
    }

    Vector3f& operator*=(float f)
    {
        x *= f;
        y *= f;
        z *= f;

        return *this;
    }

    Vector3f Cross(const Vector3f& v) const;

    Vector3f& Normalize();

    void Print() const
    {
        printf("(%.02f, %.02f, %.02f", x, y, z);
    }
};

inline Vector3f operator-(const Vector3f& l, const Vector3f& r)
{
    Vector3f Ret(l.x - r.x,
                 l.y - r.y,
                 l.z - r.z);

    return Ret;
}

inline Vector3f operator*(const Vector3f& l, float f)
{
    Vector3f Ret(l.x * f,
                 l.y * f,
                 l.z * f);

    return Ret;
}


class Matrix4f
{
public:
    float m[4][4];

    Matrix4f()
    {        
    }


    inline void InitIdentity()
    {
        m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = 0.0f;
        m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = 0.0f;
        m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = 0.0f;
        m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
    }

    inline Matrix4f operator*(const Matrix4f& Right) const
    {
        Matrix4f Ret;

        for (unsigned int i = 0 ; i < 4 ; i++) {
            for (unsigned int j = 0 ; j < 4 ; j++) {
                Ret.m[i][j] = m[i][0] * Right.m[0][j] +
                              m[i][1] * Right.m[1][j] +
                              m[i][2] * Right.m[2][j] +
                              m[i][3] * Right.m[3][j];
            }
        }

        return Ret;
    }

    void updateVector(Vector3f& Left, Vector3f& Right)
    {
        Left.x = m[0][0] * Right.x +
                m[0][1] * Right.y +
                m[0][2] * Right.z +
                m[0][3] * 1 ;
        Left.y = m[1][0] * Right.x +
                m[1][1] * Right.y +
                m[1][2] * Right.z +
                m[1][3] * 1;
        Left.z = m[2][0] * Right.x +
                m[2][1] * Right.y +
                m[2][2] * Right.z +
                m[2][3] * 1;
    }

    void InitScaleTransform(float ScaleX, float ScaleY, float ScaleZ);
    void InitRotateTransform(float RotateX, float RotateY, float RotateZ);
    void InitTranslationTransform(float x, float y, float z);
    void InitCameraTransform(const Vector3f& Target, const Vector3f& Up);
    void InitPersProjTransform(float FOV, float Width, float Height, float zNear, float zFar);
};

void setPseudoPerspective(Vector3f &vertex);


#endif	/* MATH_3D_H */

