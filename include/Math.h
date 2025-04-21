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
	inline float   operator *(const Vec2<t> &V) const { return u*V.u + v*V.v; }
	inline bool	   operator ==(const Vec2<t> &V) const { return x == V.x && y == V.y; }

	// Methods

	float     length () const { return std::sqrt(x*x+y*y); }
	Vec2<t> & normalize(t l=1) { *this = (*this)*(l/length()); return *this; }


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
	Vec3(t _s) : x(_s), y(_s), z(_s) {} // WARNING: will implicitly cast t to Vec3f
	Vec3(t _x, t _y, t _z) : x(_x), y(_y), z(_z) {}

	// Operators

	inline Vec3<t> operator ^(const Vec3<t> &v) const { return Vec3<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
	inline Vec3<t> operator +(const Vec3<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z); }
	inline Vec3<t> operator -(const Vec3<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z); }
	inline Vec3<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f); } // should not have this operator causes bugs, do something about this!!!
	inline t       operator *(const Vec3<t> &v) const { return x*v.x + y*v.y + z*v.z; }
	inline bool    operator ==(const Vec3<t> &v) const { return x == v.x && y == v.y && z == v.z; }

	// Methods

	float     length () const { return std::sqrt(x*x+y*y+z*z); }
	Vec3<t> & normalize(t l=1) { *this = (*this)*(l/length()); return *this; }
	Vec2<t> xy() const { return Vec2<t>(x, y); }
	Vec3<t> static hadamard_product(const Vec3<t>& a, const Vec3<t>& b) { return Vec3<t>(a.x*b.x, a.y*b.y, a.z*b.z); }

	// For debugging

	template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
};


template <class t> struct Vec4 
{
	union 
	{
		struct { t x, y, z, w; };
		t raw[4];
	};

	// Constructors

	Vec4() : x(0), y(0), z(0), w(0) {}
	Vec4(Vec3<t> v, t w) : x(v.x), y(v.y), z(v.z), w(w) {}
	Vec4(t _x, t _y, t _z, t _w) : x(_x), y(_y), z(_z), w(_w) {}

	// Operators

	inline Vec4<t> operator +(const Vec4<t> &v) const { return Vec4<t>(x+v.x, y+v.y, z+v.z, w+v.w); }
	inline Vec4<t> operator -(const Vec4<t> &v) const { return Vec4<t>(x-v.x, y-v.y, z-v.z, w-v.w); }
	inline Vec4<t> operator *(float f)          const { return Vec4<t>(x*f, y*f, z*f, w*f); } // should not have this operator 

	// Methods

	Vec3<t> xyz() const { return Vec3<t>(x, y, z); }
	// float     length () const { return std::sqrt(x*x+y*y+z*z); }
	// Vec3<t> & normalize(t l=1) { *this = (*this)*(l/length()); return *this; }

	// For debugging

	// template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
};


typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec3<int>   Vec3i;
typedef Vec4<float> Vec4f;


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

	Mat3x3f operator*(float s)
	{
		return Mat3x3f(col(0) * s, col(1) * s, col(2) * s);
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

	Mat3x3f inv()
	{
		return adj() * (1.0f / determinant());
	}

	Mat3x3f adj()
	{
		return cofactor().transposed();
	}

	Mat3x3f cofactor()
	{
		float a = mat[0][0];
		float b = mat[0][1];
		float c = mat[0][2];

		float d = mat[1][0];
		float e = mat[1][1];
		float f = mat[1][2];

		float g = mat[2][0];
		float h = mat[2][1];
		float i = mat[2][2];

		return Mat3x3f(
			Vec3f(e*i - f*h, c*h - b*i, b*h - c*e),
			Vec3f(f*g - d*i, a*i - c*g, b*g - a*h),
			Vec3f(d*f - e*g, c*d - a*f, a*e - b*d)
		);
	}

	Mat3x3f transposed()
	{
		Mat3x3f trans;

		trans.mat[0][0] = mat[0][0];
        trans.mat[0][1] = mat[1][0];
        trans.mat[0][2] = mat[2][0];

        trans.mat[1][0] = mat[0][1];
        trans.mat[1][1] = mat[1][1];
        trans.mat[1][2] = mat[2][1];

        trans.mat[2][0] = mat[0][2];
        trans.mat[2][1] = mat[1][2];
        trans.mat[2][2] = mat[2][2];

		return trans;
	}
};


struct Mat4x4f
{
	Vec4f cols[4]; // this might have been bad idea

	// Constructors

    Mat4x4f() 
	{
		cols[0] = Vec4f(1.0f, 0.0f, 0.0f, 0.0f);
		cols[1] = Vec4f(0.0f, 1.0f, 0.0f, 0.0f);
		cols[2] = Vec4f(0.0f, 0.0f, 1.0f, 0.0f);
		cols[3] = Vec4f(0.0f, 0.0f, 0.0f, 1.0f);
	}
    Mat4x4f(Vec4f x, Vec4f y, Vec4f z, Vec4f w, bool transpose = false)
	{
		if (!transpose)
		{
			cols[0] = x; cols[1] = y; cols[2] = z; cols[3] = w;
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				cols[i].x = x.raw[i];
				cols[i].y = y.raw[i];
				cols[i].z = z.raw[i];
				cols[i].w = w.raw[i];
			}
		}
	}

	// Operators

	Vec4f operator*(Vec4f v)
	{
		return (cols[0] * v.x) + (cols[1] * v.y) + (cols[2] * v.z) + (cols[3] * v.w); 
	}


	Mat4x4f operator*(Mat4x4f m)
	{
		Mat4x4f prod;

		for (int i = 0; i < 4; i++)
		{
			prod.cols[i] = (this->cols[0] * m.cols[i].x) + (this->cols[1] * m.cols[i].y) + (this->cols[2] * m.cols[i].z) + (this->cols[3] * m.cols[i].w);   
		}

		return prod;
	}

	// Methods

	Mat3x3f truncated()
	{
		return Mat3x3f(cols[0].xyz(), cols[1].xyz(), cols[2].xyz());
	}
};

float cross(Vec2f a, Vec2f b);

Vec3f reflect(Vec3f n, Vec3f l);
Vec3f get_triangle_normal(Vec3f a, Vec3f b, Vec3f c);

float max(float a, float b);
float min(float a, float b);
float clampedf(float a, float min, float max);
Vec3f clampedVec3f(Vec3f v, float min, float max);
Vec2f clampedVec2f(Vec2f v, float min, float max);

Mat4x4f get_transformation(Vec3f pos, float scale);
Mat4x4f look_at(Vec3f pos, Vec3f at, Vec3f up);

float power(float a, int b);
Vec3f component_wise_product_vec3f(Vec3f a, Vec3f b);
Vec2f component_wise_product_vec2f(Vec2f a, Vec2f b);

Vec3f interpolate_barycentric_vec3f(Vec3f a, Vec3f b, Vec3f c, Vec3f weights);
Vec2f interpolate_barycentric_vec2f(Vec2f a, Vec2f b, Vec2f c, Vec3f weights);
float interpolate_barycentric_f    (float a, float b, float c, Vec3f weights);