#ifndef VMATH_H
#define VMATH_H

#include <cmath>

const double PI = 3.14159265358979323846;
const double SEC24H = 24*60*60.0;
const double DEG2RAD = PI / 180.0;
const double SEC2RAD = 2*PI / SEC24H;
const double RAD2DEG = 180.0 / PI;

template<typename T>
T clamp(T x, T xmin, T xmax)
{
    return x < xmin ? xmin : x > xmax ? xmax : x;
}
template<typename T>
T lint(T x, T x0, T x1, T y0, T y1)
{
    if(x0 > x1){
        std::swap(x0, x1);
    }
    return
        (x < x0) ? y0 :
        (x > x1) ? y1 :
        y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}
template<typename T>
T lerp(T x, T y, T a)
{
    return (1 - a) * x + a * y;
}

template<typename T>
struct Vec4
{
    using Element = T;
    Element x, y, z, w;
    Vec4(Element x_, Element y_, Element z_, Element w_)
        : x(x_), y(y_), z(z_), w(w_)
    {}
};
using Vec4F = Vec4<float>;
using Vec4D = Vec4<double>;

template<typename T>
struct Mat4
{
    using Element = T;
    // elements(row col)
    Element e00,e10,e20,e30, e01,e11,e21,e31, e02,e12,e22,e32, e03,e13,e23,e33;
    Mat4(Element v00,Element v10,Element v20,Element v30, Element v01,Element v11,Element v21,Element v31, Element v02,Element v12,Element v22,Element v32, Element v03,Element v13,Element v23,Element v33)
        : e00(v00),e10(v10),e20(v20),e30(v30), e01(v01),e11(v11),e21(v21),e31(v31), e02(v02),e12(v12),e22(v22),e32(v32), e03(v03),e13(v13),e23(v23),e33(v33)
    {}
    static Mat4 identity()
    {
        return Mat4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
    }
    static Mat4 rotZ(Element rad)
    {
        return Mat4(std::cos(rad), std::sin(rad), 0, 0,
                    -std::sin(rad), std::cos(rad), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);
    }
    static Mat4 rotY(Element rad)
    {
        return Mat4(std::cos(rad), 0, -std::sin(rad), 0,
                    0, 1, 0, 0,
                    std::sin(rad), 0, std::cos(rad), 0,
                    0, 0, 0, 1);
    }
    static Mat4 rotX(Element rad)
    {
        return Mat4(1, 0, 0, 0,
                    0, std::cos(rad), std::sin(rad), 0,
                    0, -std::sin(rad),std::cos(rad), 0,
                    0, 0, 0, 1);
    }
    static Mat4 perspective(Element fovYDeg, Element screenW, Element screenH, Element nearZ, Element farZ)
    {
        const Element fovYRad = fovYDeg * DEG2RAD;
        const Element aspectRatio = screenH / screenW;

        const Element h = 1.0 / std::tan(fovYRad / 2.0);
        const Element w = h * aspectRatio;
        const Element zNearFar = nearZ - farZ;
        return Mat4(
            w, 0, 0, 0,
            0, h, 0, 0,
            0, 0, (farZ+nearZ)/zNearFar, -1,
            0, 0, 2*nearZ*farZ/zNearFar, 0);
    }
    static Mat4 translate(Element dx, Element dy, Element dz)
    {
        return Mat4(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            dx, dy, dz, 1);
    }
    static Mat4 scale(Element sx, Element sy, Element sz)
    {
        return Mat4(
            sx, 0, 0, 0,
            0, sy, 0, 0,
            0, 0, sz, 0,
            0, 0, 0, 1);
    }
};
using Mat4F = Mat4<float>;
using Mat4D = Mat4<double>;
template<typename T>
Mat4<T> operator*(const Mat4<T> &lhs, const Mat4<T> &rhs)
{
    return Mat4<T>(
        lhs.e00*rhs.e00 + lhs.e01*rhs.e10 + lhs.e02*rhs.e20 + lhs.e03*rhs.e30,
        lhs.e10*rhs.e00 + lhs.e11*rhs.e10 + lhs.e12*rhs.e20 + lhs.e13*rhs.e30,
        lhs.e20*rhs.e00 + lhs.e21*rhs.e10 + lhs.e22*rhs.e20 + lhs.e23*rhs.e30,
        lhs.e30*rhs.e00 + lhs.e31*rhs.e10 + lhs.e32*rhs.e20 + lhs.e33*rhs.e30,

        lhs.e00*rhs.e01 + lhs.e01*rhs.e11 + lhs.e02*rhs.e21 + lhs.e03*rhs.e31,
        lhs.e10*rhs.e01 + lhs.e11*rhs.e11 + lhs.e12*rhs.e21 + lhs.e13*rhs.e31,
        lhs.e20*rhs.e01 + lhs.e21*rhs.e11 + lhs.e22*rhs.e21 + lhs.e23*rhs.e31,
        lhs.e30*rhs.e01 + lhs.e31*rhs.e11 + lhs.e32*rhs.e21 + lhs.e33*rhs.e31,

        lhs.e00*rhs.e02 + lhs.e01*rhs.e12 + lhs.e02*rhs.e22 + lhs.e03*rhs.e32,
        lhs.e10*rhs.e02 + lhs.e11*rhs.e12 + lhs.e12*rhs.e22 + lhs.e13*rhs.e32,
        lhs.e20*rhs.e02 + lhs.e21*rhs.e12 + lhs.e22*rhs.e22 + lhs.e23*rhs.e32,
        lhs.e30*rhs.e02 + lhs.e31*rhs.e12 + lhs.e32*rhs.e22 + lhs.e33*rhs.e32,

        lhs.e00*rhs.e03 + lhs.e01*rhs.e13 + lhs.e02*rhs.e23 + lhs.e03*rhs.e33,
        lhs.e10*rhs.e03 + lhs.e11*rhs.e13 + lhs.e12*rhs.e23 + lhs.e13*rhs.e33,
        lhs.e20*rhs.e03 + lhs.e21*rhs.e13 + lhs.e22*rhs.e23 + lhs.e23*rhs.e33,
        lhs.e30*rhs.e03 + lhs.e31*rhs.e13 + lhs.e32*rhs.e23 + lhs.e33*rhs.e33);
}
template<typename T>
Vec4<T> operator*(const Mat4<T> &lhs, const Vec4<T> &rhs)
{
    return Vec4<T>(
        lhs.e00*rhs.x + lhs.e01*rhs.y + lhs.e02*rhs.z + lhs.e03*rhs.w,
        lhs.e10*rhs.x + lhs.e11*rhs.y + lhs.e12*rhs.z + lhs.e13*rhs.w,
        lhs.e20*rhs.x + lhs.e21*rhs.y + lhs.e22*rhs.z + lhs.e23*rhs.w,
        lhs.e30*rhs.x + lhs.e31*rhs.y + lhs.e32*rhs.z + lhs.e33*rhs.w);
}

#endif//VMATH_H
