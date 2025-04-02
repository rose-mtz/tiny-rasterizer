#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include "Image.h"
#include "Mesh.h"


/**
 * Should I add:
 * struct Vertex { color, position, normal }
 * struct Triangle { Vertex v1, v2, v3 }
 * maybe later 
 */


bool back_facing_triangle(Vec3f a, Vec3f b, Vec3f c);
void draw_line(Vec2i p0, Vec2i p1, Image& image, Vec3f color);
void draw_triangle(Vec2i a, Vec2i b, Vec2i c, Image& image, Vec3f a_color, Vec3f b_color, Vec3f c_color);

int main()
{   
    Mesh mesh = load_mesh("./obj/african_head.obj");

    Image image;
    init_image(image, 100, 100);

    Vec3f white (1.0f, 1.0f, 1.0f);
    Vec3f red   (1.0f, 0.0f, 0.0f);
    Vec3f green (0.0f, 1.0f, 0.0f);
    Vec3f blue  (0.0f, 0.0f, 1.0f);

    // for (int i = 0; i < mesh.faces.size(); i++)
    // {
    //     Vec3f v0 = mesh.vertices[mesh.faces[i][0].x];
    //     Vec3f v1 = mesh.vertices[mesh.faces[i][1].x];
    //     Vec3f v2 = mesh.vertices[mesh.faces[i][2].x];

    //     if (back_facing_triangle(v0, v1, v2))
    //     {
    //         continue;
    //     }

    //     Vec2i v0_trans ((v0.x + 1.0f) * ((image.width  - 1)/ 2.0f), (v0.y + 1.0f) * ((image.height - 1) / 2.0f));
    //     Vec2i v1_trans ((v1.x + 1.0f) * ((image.width  - 1)/ 2.0f), (v1.y + 1.0f) * ((image.height - 1) / 2.0f));
    //     Vec2i v2_trans ((v2.x + 1.0f) * ((image.width  - 1)/ 2.0f), (v2.y + 1.0f) * ((image.height - 1) / 2.0f));

    //     Vec3f color_rnd;
    //     color_rnd.x = (std::rand() % 101) / 100.0f;
    //     color_rnd.y = (std::rand() % 101) / 100.0f;
    //     color_rnd.z = (std::rand() % 101) / 100.0f;

    //     draw_triangle(v0_trans, v1_trans, v2_trans, image, color_rnd);

    //     // for (int j = 0; j < 3; j++)
    //     // {
    //     //     Vec3f v0 = mesh.vertices[mesh.faces[i][j].x];
    //     //     Vec3f v1 = mesh.vertices[mesh.faces[i][(j + 1) % 3].x];

    //     //     int x0 = (v0.x + 1.0f) * ((image.width  - 1)/ 2.0f);
    //     //     int y0 = (v0.y + 1.0f) * ((image.height - 1) / 2.0f);
    //     //     int x1 = (v1.x + 1.0f) * ((image.width  - 1)/ 2.0f);
    //     //     int y1 = (v1.y + 1.0f) * ((image.height - 1) / 2.0f);

    //     //     draw_line(Vec2i(x0, y0), Vec2i(x1, y1), image, white);
    //     // }
    // }

    Vec2i a(0,0), b(82,34), c(24,99);
    draw_triangle(a, b, c, image, red, blue, green);
    // draw_line(a, b, image, red);
    // draw_line(a, c, image, red);
    // draw_line(b, c, image, red);

    save_image("./image.bmp", image);

    return 0;
}


void draw_line(Vec2i p0, Vec2i p1, Image& image, Vec3f color)
{
    assert(p0.x >= 0 && p0.x < image.width);
    assert(p0.y >= 0 && p0.y < image.height);
    assert(p1.x >= 0 && p1.x < image.width);
    assert(p1.y >= 0 && p1.y < image.height);

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
        image.image[current_pixel.y][current_pixel.x] = color;
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


bool back_facing_triangle(Vec3f a, Vec3f b, Vec3f c)
{
    Vec3f v1 = subtract(b, a);
    Vec3f v2 = subtract(c, a);
    Vec3f normal_ish = cross(v1,v2);
    
    // WARNING: assumes camera is facing down z-ed axis
    return dot(normal_ish, Vec3f(0.0f, 0.0f, 1.0f)) <= 0.0f;
}


void draw_triangle(Vec2i a, Vec2i b, Vec2i c, Image& image, Vec3f a_color, Vec3f b_color, Vec3f c_color)
{
    Vec2i min;
    Vec2i max;
    // bounding box for triangle
    min.x = std::min(a.x, std::min(b.x, c.x));
    min.y = std::min(a.y, std::min(b.y, c.y));
    max.x = std::max(a.x, std::max(b.x, c.x));
    max.y = std::max(a.y, std::max(b.y, c.y));
    // clip box to be within image
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, image.width - 1);
    max.y = std::min(max.y, image.height - 1);

    // center of pixels
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
            Mat3x3 D_x (Vec3f(pixel_center, 1), Vec3f(b_center, 1), Vec3f(c_center, 1));
            Mat3x3 D_y (Vec3f(a_center, 1), Vec3f(pixel_center, 1), Vec3f(c_center, 1));
            Mat3x3 D_z (Vec3f(a_center, 1), Vec3f(b_center, 1), Vec3f(pixel_center, 1));

            float determinant_of_dx = determinant(D_x);
            float determinant_of_dy = determinant(D_y);
            float determinant_of_dz = determinant(D_z);

            float alpha = determinant_of_dx / determinant_of_d;
            float beta = determinant_of_dy / determinant_of_d;
            float gamma = determinant_of_dz / determinant_of_d;

            // assert(alpha + beta + gamma == 1.0f); // FAILS sometimes!

            if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f)
            {
                Vec3f inter_color = (a_color * alpha) + (b_color * beta) + (c_color * gamma);
                image.image[row][col] = inter_color;
            }
        }
    }
}