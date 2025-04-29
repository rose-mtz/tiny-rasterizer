#include "Rasterize.h"


std::vector<QuadPixel> rasterize_quad(Vec2f A, Vec2f B, Vec2f C, Vec2f D, Vec2i device)
{
    std::vector<QuadPixel> pixels;

    // Bounding box (pixel grid coordinate system)
    Vec2i min;
    Vec2i max;
    min.x = std::min(A.x, std::min(B.x, std::min(C.x, D.x)));
    min.y = std::min(A.y, std::min(B.y, std::min(C.y, D.y)));
    max.x = std::max(A.x, std::max(B.x, std::max(C.x, D.x)));
    max.y = std::max(A.y, std::max(B.y, std::max(C.y, D.y)));
    // Clip box to be within pixel grid
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, device.x - 1);
    max.y = std::min(max.y, device.y - 1);

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f P (col + 0.5f, row + 0.5f);
            float u = -1.0f;
            float v = -1.0f;

            float a = cross(A - B, A - B - D + C);
            float b = cross(P - A, A - B - D + C) + cross(A - B, D - A);
            float c = cross(P - A, D - A);

            if (a == 0) 
            {
                // Solve for u (degenerates to linear equation)
                u = -c / b;
            }
            else
            {
                // Solve for u (quadratic equation)
                float discriminant = (b * b) - (4 * a * c);
                if (discriminant < 0) continue;

                float u_plus = (-b + std::sqrt(discriminant)) / (2 * a);
                float u_minus = (-b - std::sqrt(discriminant)) / (2 * a);

                if (u_plus  >= 0.0f && u_plus  <= 1.0f) u = u_plus;
                if (u_minus >= 0.0f && u_minus <= 1.0f) u = u_minus;
            }

            // Solve for v
            float numerator_x = (P.x - A.x + (u * A.x) - (u * B.x));
            float numerator_y = (P.y - A.y + (u * A.y) - (u * B.y));
            float denominator_x = (-A.x + (u * A.x) - (u * B.x) + D.x - (u * D.x) + (u * C.x));
            float denominator_y = (-A.y + (u * A.y) - (u * B.y) + D.y - (u * D.y) + (u * C.y));
            v = (denominator_x != 0.0f) ? numerator_x / denominator_x : numerator_y / denominator_y;

            // Outside quad check
            if (u < 0.0f || u > 1.0f || v < 0.0f || v > 1.0f) continue;

            float alpha  = (1.0f - u) * (1.0f - v);
            float beta   = u * (1.0f - v);
            float gamma  = u * v;
            float lambda = v * (1.0f - u);

            QuadPixel pix;
            pix.pixel  = Vec2i(col, row);
            pix.barycentric = Vec4f(alpha, beta, gamma, lambda);
            pix.uv = Vec2f(u, v);

            pixels.push_back(pix);
        }
    }

    return pixels;
}


// Rasterize line
// Points must be in device coordinates
// Thickness is in pixels
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness, Vec2i device)
{   
    std::vector<LinePixel> pixels;

    // Bounding box (includes thickness)
    Vec2i min, max;
    min.x = std::min(p0.x, p1.x) - thickness;
    min.y = std::min(p0.y, p1.y) - thickness;
    max.x = std::max(p0.x, p1.x) + thickness;
    max.y = std::max(p0.y, p1.y) + thickness;

    // Clip bounding box to device screen
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, device.x - 1);
    max.y = std::min(max.y, device.y - 1);

    Vec2f dir = Vec2f(p1.x - p0.x, p1.y - p0.y).normalize();
    float line_length = (p1 - p0).length();
    float half_thickness = thickness / 2.0f;

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f p (col + 0.5f, row + 0.5f);
            float projected_distance = (p - p0) * dir;

            if (projected_distance < 0.0f || projected_distance > line_length)
            {
                continue;
            }

            float perp_distance = ((p - p0) - (dir * projected_distance)).length();
            float t = projected_distance / line_length;

            if (perp_distance <= half_thickness)
            {
                pixels.push_back(LinePixel {Vec2i(col, row), t});
            }
        }
    }

    return pixels;
}


// Rasterize a triangle
// Points must be in device coordinates
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c, Vec2i device)
{
    std::vector<TrianglePixel> pixels;

    // Bounding box for triangle (pixel grid coordinate system)
    Vec2i min;
    Vec2i max;
    min.x = std::min(a.x, std::min(b.x, c.x));
    min.y = std::min(a.y, std::min(b.y, c.y));
    max.x = std::max(a.x, std::max(b.x, c.x));
    max.y = std::max(a.y, std::max(b.y, c.y));
    // Clip box to be within pixel grid
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, device.x - 1);
    max.y = std::min(max.y, device.y - 1);

    Mat3x3f D (Vec3f(a.x, a.y, 1.0f), Vec3f(b.x, b.y, 1.0f), Vec3f(c.x, c.y, 1.0f));
    float determinant = D.determinant();

    // Degenerate triangle check
    if (std::abs(determinant) < 1e-5f) return pixels;

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f pixel_center (col + 0.5f, row + 0.5f);

            // Cramer's rule
            Vec3f B (pixel_center.x, pixel_center.y, 1.0f);
            float alpha = Mat3x3f(B,        D.col(1), D.col(2)).determinant() / determinant;
            float beta  = Mat3x3f(D.col(0), B,        D.col(2)).determinant() / determinant;
            float gamma = Mat3x3f(D.col(0), D.col(1), B       ).determinant() / determinant;

            // Inside triangle check
            if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f)
            {
                pixels.push_back(TrianglePixel {Vec2i(col, row), Vec3f(alpha, beta, gamma)});
            }
        }
    }

    return pixels;
}