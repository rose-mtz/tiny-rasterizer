#include "Math.h"


Vec3f get_triangle_normal(Vec3f a, Vec3f b, Vec3f c)
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


Mat4x4f get_transformation(Vec3f pos, float scale)
{
    Mat4x4f mat;
    mat.cols[0] = mat.cols[0] * scale;
    mat.cols[1] = mat.cols[1] * scale;
    mat.cols[2] = mat.cols[2] * scale;
    mat.cols[3] = Vec4f(pos, 1.0f);

    return mat;
}


Mat4x4f look_at(Vec3f pos, Vec3f at, Vec3f up)
{
    Vec3f z = (pos - at).normalize(); // 'backwards' of camera
    Vec3f x = (up ^ z).normalize(); // right of camera
    Vec3f y = (z ^ x).normalize(); // up of camera

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