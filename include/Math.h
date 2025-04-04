#pragma once

#include <cmath>
#include <iostream>


/**
 * CREDIT: Dmitry V. Sokolov
 * REPO: https://github.com/ssloy/tinyrenderer/tree/a175be75a8a9a773bdfae7543a372e3bc859e02f
 */


template <class t> struct Vec2
{
	union 
	{
		struct { t u, v; };
		struct { t x, y; };
		t raw[2];
	};

	// Constructors

	Vec2() : u(0), v(0) {}
	Vec2(t _u, t _v) : u(_u), v(_v) {}

	// Operators

	inline Vec2<t> operator +(const Vec2<t> &V) const { return Vec2<t>(u+V.u, v+V.v); }
	inline Vec2<t> operator -(const Vec2<t> &V) const { return Vec2<t>(u-V.u, v-V.v); }
	inline Vec2<t> operator *(float f)          const { return Vec2<t>(u*f, v*f); }

	// For debuging

	template <class > friend std::ostream& operator<<(std::ostream& s, Vec2<t>& v);
};


template <class t> struct Vec3 
{
	union 
	{
		struct { t x, y, z; };
		struct { t ivert, iuv, inorm; };
		t raw[3];
	};

	// Constructors

	Vec3() : x(0), y(0), z(0) {}
	Vec3(t _s) : x(_s), y(_s), z(_s) {}
	Vec3(t _x, t _y, t _z) : x(_x), y(_y), z(_z) {}

	// Operators

	inline Vec3<t> operator ^(const Vec3<t> &v) const { return Vec3<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
	inline Vec3<t> operator +(const Vec3<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z); }
	inline Vec3<t> operator -(const Vec3<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z); }
	inline Vec3<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f); }
	inline t       operator *(const Vec3<t> &v) const { return x*v.x + y*v.y + z*v.z; }

	// Methods

	float     length () const { return std::sqrt(x*x+y*y+z*z); }
	Vec3<t> & normalize(t l=1) { *this = (*this)*(l/length()); return *this; }

	// For debugging

	template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
};


typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec3<int>   Vec3i;


template <class t> std::ostream& operator<<(std::ostream& s, Vec2<t>& v) {
	s << "(" << v.x << ", " << v.y << ")\n";
	return s;
}


template <class t> std::ostream& operator<<(std::ostream& s, Vec3<t>& v) {
	s << "(" << v.x << ", " << v.y << ", " << v.z << ")\n";
	return s;
}


struct Mat3x3f
{
    // row-major
    float mat[3][3];

	// Constructors

    Mat3x3f() {}
    Mat3x3f(Vec3f x, Vec3f y, Vec3f z) 
    {
        mat[0][0] = x.x;
        mat[1][0] = x.y;
        mat[2][0] = x.z;

        mat[0][1] = y.x;
        mat[1][1] = y.y;
        mat[2][1] = y.z;

        mat[0][2] = z.x;
        mat[1][2] = z.y;
        mat[2][2] = z.z;
    }

	// Operators

	Vec3f operator*(Vec3f v)
	{
		return Vec3f(
			v.x * mat[0][0] + v.y * mat[0][1] + v.z * mat[0][2],
			v.x * mat[1][0] + v.y * mat[1][1] + v.z * mat[1][2],
			v.x * mat[2][0] + v.y * mat[2][1] + v.z * mat[2][2]
		);
	}

	// Methods

	float determinant()
	{
		return (
			mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
			mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) + 
			mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0])
		);
	}

	Vec3f col(int i)
	{
		return Vec3f(mat[0][i], mat[1][i], mat[2][i]);
	}
};


Vec3f get_triangle_normal(Vec3f a, Vec3f b, Vec3f c);
float clampedf(float a, float min, float max);
Vec3f clampedVec3f(Vec3f v, float min, float max);