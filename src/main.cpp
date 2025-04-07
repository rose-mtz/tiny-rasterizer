#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <limits>
#include "tgaimage.h"
#include "model.h"

const Vec3f WHITE (1.0f, 1.0f, 1.0f);
const Vec3f RED   (1.0f, 0.0f, 0.0f);
const Vec3f GREEN (0.0f, 1.0f, 0.0f);
const Vec3f BLUE  (0.0f, 0.0f, 1.0f);
const Vec3f BLACK (0.0f, 0.0f, 0.0f);

const TGAColor TGA_WHITE = TGAColor(255, 255, 255, 255);
const TGAColor TGA_RED   = TGAColor(255, 0,   0,   255);
const TGAColor TGA_GREEN = TGAColor(0,   255, 0,   255);

const int   device_width  = 1000;
const int   device_height = 1000;
const float virtual_screen_width = 0.65;
const float virtual_screen_height = 0.65;

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


TGAColor to_tgacolor(Vec3f color);
Vec3f to_vec3fcolor(TGAColor color);
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness = 1.0f);
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c);

int main()
{   
    TGAImage image(device_width, device_height, TGAImage::RGB);
    TGAImage texture;
    texture.read_tga_file("./obj/african_head_diffuse.tga");
    texture.flip_vertically();
    Model model ("./obj/african_head.obj");

    float** z_buffer = new float*[device_height];
    for (int i = 0; i < device_height; i++)
    {
        z_buffer[i] = new float[device_width];
        for (int j = 0; j < device_width; j++)
        {
            z_buffer[i][j] = -std::numeric_limits<float>::max();
        }
    }

    for (int i = 0; i < model.nfaces(); i++)
    {
        std::vector<int> face = model.face(i);
        Vec3f v0 = model.vert(face[0]);
        Vec3f v1 = model.vert(face[3]);
        Vec3f v2 = model.vert(face[6]);
        Vec3f n0 = model.norm(face[2]);
        Vec3f n1 = model.norm(face[5]);
        Vec3f n2 = model.norm(face[8]);
        Vec2f uv0 = model.uv(face[1]);
        Vec2f uv1 = model.uv(face[4]);
        Vec2f uv2 = model.uv(face[7]);

        // Back facing triangle check
        if (get_triangle_normal(v0, v1, v2) * Vec3f(0.0f, 0.0f, 1.0f) <= 0.0f) // this is fine for now because the normals are the same in local, global, and camera
        {
            continue;
        }

        // Transformations: local --> global --> camera --> virtual screen --> device screen

        Vec3f v0_global = v0;
        Vec3f v1_global = v1;
        Vec3f v2_global = v2;
        Vec3f v0_camera = v0_global + Vec3f(0.0f, 0.0f, -2.0f);
        Vec3f v1_camera = v1_global + Vec3f(0.0f, 0.0f, -2.0f);
        Vec3f v2_camera = v2_global + Vec3f(0.0f, 0.0f, -2.0f);
        Vec3f v0_virtual = v0_camera;
        Vec3f v1_virtual = v1_camera;
        Vec3f v2_virtual = v2_camera;
        Vec3f v0_device = Vec3f((v0_virtual.x * ((device_width / virtual_screen_width) / 2.0f)) + (device_width / 2.0f), (v0_virtual.y * ((device_height / virtual_screen_height) / 2.0f)) + (device_height / 2.0f), v0_virtual.z);
        Vec3f v1_device = Vec3f((v1_virtual.x * ((device_width / virtual_screen_width) / 2.0f)) + (device_width / 2.0f), (v1_virtual.y * ((device_height / virtual_screen_height) / 2.0f)) + (device_height / 2.0f), v1_virtual.z);
        Vec3f v2_device = Vec3f((v2_virtual.x * ((device_width / virtual_screen_width) / 2.0f)) + (device_width / 2.0f), (v2_virtual.y * ((device_height / virtual_screen_height) / 2.0f)) + (device_height / 2.0f), v2_virtual.z);
        
        // Shading / colors

        Vec3f rnd_color_1;
        rnd_color_1.x = (std::rand() % 101) / 100.0f;
        rnd_color_1.y = (std::rand() % 101) / 100.0f;
        rnd_color_1.z = (std::rand() % 101) / 100.0f;
        Vec3f rnd_color_2;
        rnd_color_2.x = (std::rand() % 101) / 100.0f;
        rnd_color_2.y = (std::rand() % 101) / 100.0f;
        rnd_color_2.z = (std::rand() % 101) / 100.0f;
        Vec3f rnd_color_3;
        rnd_color_3.x = (std::rand() % 101) / 100.0f;
        rnd_color_3.y = (std::rand() % 101) / 100.0f;
        rnd_color_3.z = (std::rand() % 101) / 100.0f;

        Vec3f normal = get_triangle_normal(v0, v1, v2); // this is fine for now because the normals are the same in local, global, and camera
        Vec3f light_dir = Vec3f(0.0f, 0.0f, 1.0f);
        float ambient = 0.01f;
        float diffuse_gouraud = normal * light_dir;

        Vec3f flat_gray_shaded_gouraud = clampedVec3f(Vec3f(diffuse_gouraud + ambient), 0.0f, 1.0f);
        Vec3f rnd_color_1_shaded_gouraud = clampedVec3f((rnd_color_1 * diffuse_gouraud) + Vec3f(ambient), 0.0f, 1.0f); 
        Vec3f rnd_color_2_shaded_gouraud = clampedVec3f((rnd_color_2 * diffuse_gouraud) + Vec3f(ambient), 0.0f, 1.0f); 
        Vec3f rnd_color_3_shaded_gouraud = clampedVec3f((rnd_color_3 * diffuse_gouraud) + Vec3f(ambient), 0.0f, 1.0f);

        std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0_device.x, v0_device.y), Vec2f(v1_device.x, v1_device.y), Vec2f(v2_device.x, v2_device.y));
        for (int i = 0; i < raster_triangle.size(); i++)
        {
            TrianglePixel pixel = raster_triangle[i];
            float interpolated_depth = (v0_device.z * pixel.barycentric.x) + (v1_device.z * pixel.barycentric.y) + (v2_device.z * pixel.barycentric.z);

            if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
            {
                Vec3f interpolated_normal = ((n0 * pixel.barycentric.x) + (n1 * pixel.barycentric.y) + (n2 * pixel.barycentric.z)).normalize();
                Vec2f interpolated_uv = clampedVec2f((uv0 * pixel.barycentric.x) + (uv1 * pixel.barycentric.y) + (uv2 * pixel.barycentric.z), 0.0f, 1.0f);
                
                float diffuse = interpolated_normal * light_dir;
                Vec3f gray_shaded = clampedVec3f(Vec3f(diffuse + ambient), 0.0f, 1.0f);
                Vec3f interpolated_color = clampedVec3f((rnd_color_1 * pixel.barycentric.x) + (rnd_color_2 * pixel.barycentric.y) + (rnd_color_3 * pixel.barycentric.z), 0.0f, 1.0f);
                Vec3f interpolated_color_shaded = clampedVec3f((interpolated_color * diffuse) + ambient, 0.0f, 1.0f);
                Vec3f single_rnd_color_shaded = clampedVec3f((rnd_color_1 * diffuse) + ambient, 0.0f, 1.0f);
                Vec3f single_color_shaded = clampedVec3f((RED * diffuse) + ambient, 0.0f, 1.0f);

                Vec2i texture_pixel = Vec2i((texture.get_width() - 1) * interpolated_uv.x, (texture.get_height() - 1) * interpolated_uv.y);
                Vec3f texture_color = to_vec3fcolor(texture.get(texture_pixel.x, texture_pixel.y));

                Vec3f texture_color_shaded = clampedVec3f((texture_color * diffuse) + Vec3f(ambient), 0.0f, 1.0f);
                Vec3f texture_color_shaded_gouraud = clampedVec3f((texture_color * diffuse_gouraud) + Vec3f(ambient), 0.0f, 1.0f);

                // image.set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(texture_color_shaded_gouraud));
                // image.set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(texture_color_shaded));
                // image.set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(interpolated_color_shaded));
                // image.set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(single_rnd_color_shaded));
                image.set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(single_color_shaded));
                z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
            }
        }

        // draw_line(v0_device, v1_device, image, RED);
        // draw_line(v1_device, v2_device, image, RED);
        // draw_line(v2_device, v0_device, image, RED);
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    return 0;
}


Vec3f to_vec3fcolor(TGAColor color)
{
    return Vec3f(color.r, color.g, color.b) * (1.0f / 255.0f);
}


TGAColor to_tgacolor(Vec3f color)
{
    return TGAColor(int(color.x * 255.99), int(color.y * 255.99), int(color.z * 255.99), 255);
}


// Rastersizes line
// Points must be in device coordinates
// Thickness is in pixels
// KNOWN BUG: bounding box clips thickness.
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness)
{
    std::vector<LinePixel> pixels;

    // Bounding box
    Vec2i min, max;
    min.x = std::min(p0.x, p1.x);
    min.y = std::min(p0.y, p1.y);
    max.x = std::max(p0.x, p1.x);
    max.y = std::max(p0.y, p1.y);
    // Clip bounding box
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, device_width - 1);
    max.y = std::min(max.y, device_height - 1);

    Vec2f dir = Vec2f(p1.x - p0.x, p1.y - p0.y).normalize();
    float line_length = (p1 - p0).length();
    float half_thickness = thickness / 2.0f;

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f p (col + 0.5f, row + 0.5f);
            float projected_distance = (p - p0) * dir;
            float perp_distance = ((p - p0) - (dir * projected_distance)).length();

            if (perp_distance <= half_thickness)
            {
                pixels.push_back(LinePixel {Vec2i(col, row), (1.0f - (projected_distance / line_length))});
            }
        }
    }

    return pixels;
}


// Rastersizes a triangle
// Points must be in device coordinates
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c)
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
    max.x = std::min(max.x, device_width - 1);
    max.y = std::min(max.y, device_height - 1);

    Mat3x3f D (Vec3f(a.x, a.y, 1.0f), Vec3f(b.x, b.y, 1.0f), Vec3f(c.x, c.y, 1.0f));
    float determinant = D.determinant();

    if (std::abs(determinant) < 1e-5f)
    {
        std::cout << "degenerate triangle, you little bastard!\n"; // for testing, haven't seen it log once yet
        return pixels;
    };

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