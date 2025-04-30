#include <cassert>
#include "Math.h"



float cross(Vec2f a, Vec2f b)
{
    return a.x * b.y - a.y * b.x;
}


float min(float a, float b)
{
    return (a < b) ? a : b;
}


float max(float a, float b)
{
    return (a > b) ? a : b;
}


// (a)^b
float power(float a, int b)
{
    float result = 1.0f;
    int exp = (b >= 0) ? b : -b;

    while (exp > 0)
    {
        result *= a;
        exp--;
    }
    
    return (b >= 0) ? result : (1.0f / result);
}


// l is incoming vector
// n is surface normal
Vec3f reflect(Vec3f n, Vec3f l)
{
    // Credit: https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
    return (l - (n * 2 * (l * n))).normalize(); // Optimize: might not need to normalize if n and l are unit
}


Vec3f triangle_normal(Vec3f a, Vec3f b, Vec3f c)
{
    Vec3f v1 = b - a;
    Vec3f v2 = c - a;

    return (v1 ^ v2).normalize();
}


Vec3f clampedVec3f(Vec3f v, float min, float max)
{
    return Vec3f(clampedf(v.x, min, max), clampedf(v.y, min, max), clampedf(v.z, min, max));
}


Vec2f clampedVec2f(Vec2f v, float min, float max)
{
    return Vec2f(clampedf(v.x, min, max), clampedf(v.y, min, max));
}


float clampedf(float a, float min, float max)
{
    if (a < min)
    {
        return min;
    }
    else if (a > max)
    {
        return max;
    }
    else
    {
        return a;
    }
}


Vec3f barycentric_vec3f(Vec3f a, Vec3f b, Vec3f c, Vec3f weights)
{
    return (a * weights.x) + (b * weights.y) + (c * weights.z);
}


Vec2f barycentric_vec2f(Vec2f a, Vec2f b, Vec2f c, Vec3f weights)
{
    return (a * weights.x) + (b * weights.y) + (c * weights.z);
}


float barycentric_f(float a, float b, float c, Vec3f weights)
{
    return (a * weights.x) + (b * weights.y) + (c * weights.z);
}


float radians(float degree)
{
    return degree * (PI / 180.0f);
}


Mat4x4f look_at(Vec3f pos, Vec3f at, Vec3f up)
{
    Vec3f z = (pos - at).normalize(); // 'backwards' of camera
    Vec3f x = (up ^ z).normalize();   // right of camera
    Vec3f y = (z ^ x).normalize();    // up of camera

    Mat4x4f rotation_transposed = Mat4x4f(
        Vec4f(x, 0.0f),
        Vec4f(y, 0.0f),
        Vec4f(z, 0.0f),
        Vec4f(0.0f, 0.0f, 0.0f, 1.0f), true
    );

    Mat4x4f translation_inv = Mat4x4f();
    translation_inv.cols[3] = Vec4f(pos * -1.0f, 1.0f);

    return rotation_transposed * translation_inv;
}


Mat4x4f perspective(float t, float r, float b, float l, float n, float f)
{
    /**
     * CREDIT: Song Ho Ahn
     * SOURCE: https://www.songho.ca/opengl/gl_projectionmatrix.html
     */

    return Mat4x4f(
        (2.0f*n)/(r-l),           0.0f,  (r+l)/(r-l),              0.0f,
                  0.0f, (2.0f*n)/(t-b),  (t+b)/(t-b),              0.0f,
                  0.0f,           0.0f, -(f+n)/(f-n), (-2.0f*f*n)/(f-n),
                  0.0f,           0.0f,        -1.0f,              0.0f
    );
}

Mat4x4f orthographic(float t, float r, float b, float l, float n, float f)
{
    /**
     * CREDIT: Song Ho Ahn
     * SOURCE: https://www.songho.ca/opengl/gl_projectionmatrix.html
     */

    return Mat4x4f(
        2.0f/(r-l),       0.0f,        0.0f, -(r+l)/(r-l),
              0.0f, 2.0f/(t-b),        0.0f, -(t+b)/(t-b),
              0.0f,       0.0f, -2.0f/(f-n), -(f+n)/(f-n),
              0.0f,       0.0f,        0.0f,         1.0f
    );
}


Mat4x4f rotation_x(float theta)
{
    return Mat4x4f(
        1,          0,           0, 0,
        0, cos(theta), -sin(theta), 0,
        0, sin(theta),  cos(theta), 0,
        0,          0,           0, 1
    );
}


Mat4x4f rotation_y(float theta)
{
    return Mat4x4f(
         cos(theta), 0, sin(theta), 0,
                  0, 1,          0, 0,
        -sin(theta), 0, cos(theta), 0,
                  0, 0,          0, 1
    );
}


Mat4x4f rotation_z(float theta)
{
    return Mat4x4f(
        cos(theta), -sin(theta), 0, 0,
        sin(theta),  cos(theta), 0, 0,
                 0,           0, 1, 0,
                 0,           0, 0, 1
    );
}


Mat4x4f translation(Vec3f v)
{
    return Mat4x4f(
        1, 0, 0, v.x,
        0, 1, 0, v.y,
        0, 0, 1, v.z,
        0, 0, 0, 1
    );
}


Mat4x4f scale(Vec3f scale)
{
    return Mat4x4f(
        scale.x,       0,       0, 0,
              0, scale.y,       0, 0,
              0,       0, scale.z, 0,
              0,       0,       0, 1
    );
}