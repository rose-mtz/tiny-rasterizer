#include <vector>
#include "Math.h"

struct TrianglePixel
{
    Vec2i pixel;
    Vec3f barycentric;
};

struct LinePixel
{
    Vec2i pixel;
    float t;
};

struct QuadPixel
{
    Vec2i pixel;
    Vec4f barycentric;
    Vec2f uv;
};

std::vector<LinePixel>     rasterize_line    (Vec2f p0, Vec2f p1, float thickness, Vec2i device);
std::vector<QuadPixel>     rasterize_quad    (Vec2f A, Vec2f B, Vec2f C, Vec2f D,  Vec2i device);
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c,           Vec2i device);