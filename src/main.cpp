#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <limits>
#include "tgaimage.h"
#include "Mesh.h"

// const Vec3f WHITE (1.0f, 1.0f, 1.0f);
// const Vec3f RED   (1.0f, 0.0f, 0.0f);
// const Vec3f GREEN (0.0f, 1.0f, 0.0f);
// const Vec3f BLUE  (0.0f, 0.0f, 1.0f);
// const Vec3f BLACK (0.0f, 0.0f, 0.0f);

const TGAColor WHITE = TGAColor(255, 255, 255, 255);
const TGAColor RED   = TGAColor(255, 0,   0,   255);
const TGAColor GREEN = TGAColor(0,   255, 0,   255);
const int width  = 1000;
const int height = 1000;

TGAColor to_tgacolor(Vec3f color);
void draw_line(Vec2i p0, Vec2i p1, TGAImage& image, Vec3f color);
void draw_triangle(Vec2i a, Vec2i b, Vec2i c, TGAImage& image, Vec3f a_color, Vec3f b_color, Vec3f c_color, float a_depth, float b_depth, float c_depth, float** z_buffer);

int main()
{   
    TGAImage image(width, height, TGAImage::RGB);
    Mesh mesh = load_mesh("./obj/african_head.obj");

    float** z_buffer = new float*[height];
    for (int i = 0; i < height; i++)
    {
        z_buffer[i] = new float[width];
        for (int j = 0; j < width; j++)
        {
            z_buffer[i][j] = std::numeric_limits<float>::max();
        }
    }

    for (int i = 0; i < mesh.faces.size(); i++)
    {
        Vec3f v0 = mesh.vertices[mesh.faces[i][0].x];
        Vec3f v1 = mesh.vertices[mesh.faces[i][1].x];
        Vec3f v2 = mesh.vertices[mesh.faces[i][2].x];

        // Back facing triangle check
        if (dot(get_triangle_normal(v0, v1, v2), Vec3f(0.0f, 0.0f, 1.0f)) <= 0.0f)
        {
            continue;
        }

        Vec2i v0_trans ((v0.x + 1.0f) * ((width  - 1)/ 2.0f), (v0.y + 1.0f) * ((height - 1) / 2.0f));
        Vec2i v1_trans ((v1.x + 1.0f) * ((width  - 1)/ 2.0f), (v1.y + 1.0f) * ((height - 1) / 2.0f));
        Vec2i v2_trans ((v2.x + 1.0f) * ((width  - 1)/ 2.0f), (v2.y + 1.0f) * ((height - 1) / 2.0f));

        // Vec3f a_color;
        // a_color.x = (std::rand() % 101) / 100.0f;
        // a_color.y = (std::rand() % 101) / 100.0f;
        // a_color.z = (std::rand() % 101) / 100.0f;
        // Vec3f b_color;
        // b_color.x = (std::rand() % 101) / 100.0f;
        // b_color.y = (std::rand() % 101) / 100.0f;
        // b_color.z = (std::rand() % 101) / 100.0f;
        // Vec3f c_color;
        // c_color.x = (std::rand() % 101) / 100.0f;
        // c_color.y = (std::rand() % 101) / 100.0f;
        // c_color.z = (std::rand() % 101) / 100.0f;

        // Shading
        Vec3f normal = get_triangle_normal(v0, v1, v2);
        Vec3f light_dir = Vec3f(0.0f, 0.0f, 1.0f);
        float ambient = 0.00f;
        float diffuse = dot(normal, light_dir);
        Vec3f flat_gray_shading = clampedVec3f(Vec3f(diffuse + ambient), 0.0f, 1.0f);

        // a_color = clampedVec3f((a_color * diffuse) + Vec3f(ambient), 0.0f, 1.0f); 
        // b_color = clampedVec3f((b_color * diffuse) + Vec3f(ambient), 0.0f, 1.0f); 
        // c_color = clampedVec3f((c_color * diffuse) + Vec3f(ambient), 0.0f, 1.0f);        

        float cam_z = 3.0f;
        draw_triangle(v0_trans, v1_trans, v2_trans, image, flat_gray_shading, flat_gray_shading, flat_gray_shading, cam_z - v0.z, cam_z - v1.z, cam_z - v2.z, z_buffer);

        // for (int j = 0; j < 3; j++)
        // {
        //     Vec3f v0 = mesh.vertices[mesh.faces[i][j].x];
        //     Vec3f v1 = mesh.vertices[mesh.faces[i][(j + 1) % 3].x];

        //     int x0 = (v0.x + 1.0f) * ((image.width  - 1)/ 2.0f);
        //     int y0 = (v0.y + 1.0f) * ((image.height - 1) / 2.0f);
        //     int x1 = (v1.x + 1.0f) * ((image.width  - 1)/ 2.0f);
        //     int y1 = (v1.y + 1.0f) * ((image.height - 1) / 2.0f);

        //     draw_line(Vec2i(x0, y0), Vec2i(x1, y1), image, RED);
        // }
    }

    // save_image("./image.bmp", image);
    
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    return 0;
}


TGAColor to_tgacolor(Vec3f color)
{
    return TGAColor(int(color.x * 255.99), int(color.y * 255.99), int(color.z * 255.99), 255);
}


void draw_line(Vec2i p0, Vec2i p1, TGAImage& image, Vec3f color)
{
    assert(p0.x >= 0 && p0.x < width && p0.y >= 0 && p0.y < height);
    assert(p1.x >= 0 && p1.x < width && p1.y >= 0 && p1.y < height);

    // left-to-right (for symmetry)
    if (p0.x > p1.x) std::swap(p0, p1);

    Vec2i direction (p1.x - p0.x, p1.y - p0.y);
    Vec2f dt (1.0f / std::abs(direction.x), 1.0f / std::abs(direction.y)); // note: (float) 1.0 / (int) 0 == infinity
    Vec2f t_next (dt); // the next values at which a 'x' or 'y intersection' occur
    Vec2i delta_pixel;
    delta_pixel.x = direction.x > 0 ? 1 : -1;
    delta_pixel.y = direction.y > 0 ? 1 : -1;
    Vec2i current_pixel = p0;

    while (true)
    {
        image.set(current_pixel.x, current_pixel.y, to_tgacolor(color));
        if (current_pixel.x == p1.x && current_pixel.y == p1.y)
        {
            break;
        }

        if (t_next.x <= t_next.y)
        {
            current_pixel.x += delta_pixel.x;
            t_next.x += dt.x;
        }
        else // t_next.y > t_next.x
        {
            current_pixel.y += delta_pixel.y;
            t_next.y += dt.y;
        }
    }
}


void draw_triangle(Vec2i a, Vec2i b, Vec2i c, TGAImage& image, Vec3f a_color, Vec3f b_color, Vec3f c_color, float a_depth, float b_depth, float c_depth, float** z_buffer)
{
    // Bounding box for triangle
    Vec2i min;
    Vec2i max;
    min.x = std::min(a.x, std::min(b.x, c.x));
    min.y = std::min(a.y, std::min(b.y, c.y));
    max.x = std::max(a.x, std::max(b.x, c.x));
    max.y = std::max(a.y, std::max(b.y, c.y));

    // Clip box to be within image
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, width - 1);
    max.y = std::min(max.y, height - 1);

    // Center of pixels
    Vec2f pixel_half (0.5f, 0.5f);
    Vec2f a_center (a.x + pixel_half.x, a.y + pixel_half.y);
    Vec2f b_center (b.x + pixel_half.x, b.y + pixel_half.y);
    Vec2f c_center (c.x + pixel_half.x, c.y + pixel_half.y);

    Mat3x3 D (Vec3f(a_center, 1), Vec3f(b_center, 1), Vec3f(c_center, 1));
    float determinant_of_d = determinant(D);

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f pixel_center (col + pixel_half.x, row + pixel_half.y);

            // Cramer's rule
            Mat3x3 D_x (Vec3f(pixel_center, 1), Vec3f(b_center, 1), Vec3f(c_center, 1));
            Mat3x3 D_y (Vec3f(a_center, 1), Vec3f(pixel_center, 1), Vec3f(c_center, 1));
            Mat3x3 D_z (Vec3f(a_center, 1), Vec3f(b_center, 1), Vec3f(pixel_center, 1));
            float alpha = determinant(D_x) / determinant_of_d;
            float beta = determinant(D_y) / determinant_of_d;
            float gamma = determinant(D_z) / determinant_of_d;

            // Check if pixel center is inside triangle
            if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f)
            {
                float interpolated_depth = (a_depth * alpha) + (b_depth * beta) + (c_depth * gamma);

                // depth check
                if (z_buffer[row][col] > interpolated_depth)
                {
                    Vec3f interpolated_color = clampedVec3f((a_color * alpha) + (b_color * beta) + (c_color * gamma), 0.0f, 1.0f);
                    image.set(col, row, to_tgacolor(interpolated_color));
                    z_buffer[row][col] = interpolated_depth;
                }
            }
        }
    }
}