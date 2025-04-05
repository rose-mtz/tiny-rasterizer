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
const int width  = 10;
const int height = 10;
const float virtual_screen_width = 1.0;
const float virtual_screen_height = 1.0;

TGAColor to_tgacolor(Vec3f color);
// void draw_line(Vec2i p0, Vec2i p1, TGAImage& image, Vec3f color);
void draw_line(Vec3f p0, Vec3f p1, TGAImage& image, Vec3f color);
void draw_triangle(Vec3f a, Vec3f b, Vec3f c, TGAImage& image, Vec3f a_color, Vec3f b_color, Vec3f c_color, float** z_buffer);

int main()
{   
    TGAImage image(width, height, TGAImage::RGB);
    // Model model ("./obj/african_head.obj");

    // float** z_buffer = new float*[height];
    // for (int i = 0; i < height; i++)
    // {
    //     z_buffer[i] = new float[width];
    //     for (int j = 0; j < width; j++)
    //     {
    //         z_buffer[i][j] = -std::numeric_limits<float>::max();
    //     }
    // }

    // for (int i = 0; i < model.nfaces(); i++)
    // {
    //     std::vector<int> face = model.face(i);
    //     Vec3f v0 = model.vert(face[0]);
    //     Vec3f v1 = model.vert(face[1]);
    //     Vec3f v2 = model.vert(face[2]);

    //     // Back facing triangle check
    //     if (get_triangle_normal(v0, v1, v2) * Vec3f(0.0f, 0.0f, 1.0f) <= 0.0f)
    //     {
    //         continue;
    //     }

    //     // // pretty much local --> global --> screen space --> pixel space
    //     // Vec2i v0_trans ((v0.x + 1.0f) * ((width - 1)/ 2.0f), (v0.y + 1.0f) * ((height - 1) / 2.0f));
    //     // Vec2i v1_trans ((v1.x + 1.0f) * ((width - 1)/ 2.0f), (v1.y + 1.0f) * ((height - 1) / 2.0f));
    //     // Vec2i v2_trans ((v2.x + 1.0f) * ((width - 1)/ 2.0f), (v2.y + 1.0f) * ((height - 1) / 2.0f));

    //     // transformation: local --> global
    //         // used to scale, rotate, postion object in world
    //     Vec3f v0_global = v0;
    //     Vec3f v1_global = v1;
    //     Vec3f v2_global = v2;
    //     // transformation: global --> camera
    //         // used to get camera view of world
    //     Vec3f v0_cam = v0_global + Vec3f(0.0f, 0.0f, -2.0f);
    //     Vec3f v1_cam = v1_global + Vec3f(0.0f, 0.0f, -2.0f);
    //     Vec3f v2_cam = v2_global + Vec3f(0.0f, 0.0f, -2.0f);
    //     // transformation: camera --> projection to virtual screen
    //         // project 3d world onto virtual screen
    //     Vec3f v0_virtual = v0_cam;
    //     Vec3f v1_virtual = v1_cam;
    //     Vec3f v2_virtual = v2_cam;
    //     // transformation: virtual screen --> device coordinates
    //         // map virtual screen to actual device screen
    //         // virtual_screen is 2x2, origin is center of it
    //     float x_scale = (width / virtual_screen_width) / 2.0f;
    //     float y_scale = (height / virtual_screen_height) / 2.0f;
    //     Vec3f v0_device = Vec3f((v0_virtual.x * x_scale) + (width / 2.0f), (v0_virtual.y * y_scale) + (height / 2.0f), v0_virtual.z);
    //     Vec3f v1_device = Vec3f((v1_virtual.x * x_scale) + (width / 2.0f), (v1_virtual.y * y_scale) + (height / 2.0f), v1_virtual.z);
    //     Vec3f v2_device = Vec3f((v2_virtual.x * x_scale) + (width / 2.0f), (v2_virtual.y * y_scale) + (height / 2.0f), v2_virtual.z);
        
    //     // send to rasterizer

    //     // Vec3f a_color;
    //     // a_color.x = (std::rand() % 101) / 100.0f;
    //     // a_color.y = (std::rand() % 101) / 100.0f;
    //     // a_color.z = (std::rand() % 101) / 100.0f;
    //     // Vec3f b_color;
    //     // b_color.x = (std::rand() % 101) / 100.0f;
    //     // b_color.y = (std::rand() % 101) / 100.0f;
    //     // b_color.z = (std::rand() % 101) / 100.0f;
    //     // Vec3f c_color;
    //     // c_color.x = (std::rand() % 101) / 100.0f;
    //     // c_color.y = (std::rand() % 101) / 100.0f;
    //     // c_color.z = (std::rand() % 101) / 100.0f;

    //     // Shading
    //     Vec3f normal = get_triangle_normal(v0, v1, v2);
    //     Vec3f light_dir = Vec3f(0.0f, 0.0f, 1.0f);
    //     float ambient = 0.0f;
    //     float diffuse = normal * light_dir;
    //     Vec3f flat_gray_shading = clampedVec3f(Vec3f(diffuse + ambient), 0.0f, 1.0f);

    //     // a_color = clampedVec3f((a_color * diffuse) + Vec3f(ambient), 0.0f, 1.0f); 
    //     // b_color = clampedVec3f((b_color * diffuse) + Vec3f(ambient), 0.0f, 1.0f); 
    //     // c_color = clampedVec3f((c_color * diffuse) + Vec3f(ambient), 0.0f, 1.0f);        

    //     float cam_z = 3.0f;
    //     // draw_triangle(v0_trans, v1_trans, v2_trans, image, flat_gray_shading, flat_gray_shading, flat_gray_shading, cam_z - v0.z, cam_z - v1.z, cam_z - v2.z, z_buffer);
    //     draw_triangle(v0_device, v1_device, v2_device, image, flat_gray_shading, flat_gray_shading, flat_gray_shading, z_buffer);

    //     draw_line(v0_device, v1_device, image, RED);
    //     draw_line(v1_device, v2_device, image, RED);
    //     draw_line(v2_device, v0_device, image, RED);

    //     // for (int j = 0; j < 3; j++)
    //     // {
    //     //     Vec3f v0 = mesh.vertices[mesh.faces[i][j].x];
    //     //     Vec3f v1 = mesh.vertices[mesh.faces[i][(j + 1) % 3].x];

    //     //     int x0 = (v0.x + 1.0f) * ((image.width  - 1)/ 2.0f);
    //     //     int y0 = (v0.y + 1.0f) * ((image.height - 1) / 2.0f);
    //     //     int x1 = (v1.x + 1.0f) * ((image.width  - 1)/ 2.0f);
    //     //     int y1 = (v1.y + 1.0f) * ((image.height - 1) / 2.0f);

    //     //     draw_line(Vec2i(x0, y0), Vec2i(x1, y1), image, RED);
    //     // }
    // }
    
    Vec3f a (0,0,0), b (10,10,0), c (0, 10, 0), d (10, 0, 0);
    // draw_line(a, b, image, RED);
    draw_line(c, d, image, WHITE);

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    return 0;
}


TGAColor to_tgacolor(Vec3f color)
{
    return TGAColor(int(color.x * 255.99), int(color.y * 255.99), int(color.z * 255.99), 255);
}


// Old way of drawing lines
// void draw_line(Vec2i p0, Vec2i p1, TGAImage& image, Vec3f color)
// {
//     assert(p0.x >= 0 && p0.x < width && p0.y >= 0 && p0.y < height);
//     assert(p1.x >= 0 && p1.x < width && p1.y >= 0 && p1.y < height);

//     // left-to-right (for symmetry)
//     if (p0.x > p1.x) std::swap(p0, p1);

//     Vec2i direction (p1.x - p0.x, p1.y - p0.y);
//     Vec2f dt (1.0f / std::abs(direction.x), 1.0f / std::abs(direction.y)); // note: (float) 1.0 / (int) 0 == infinity
//     Vec2f t_next (dt); // the next values at which a 'x' or 'y intersection' occur
//     Vec2i delta_pixel;
//     delta_pixel.x = direction.x > 0 ? 1 : -1; // shouldn't this always be 1, due to left-to-right symmetry?
//     delta_pixel.y = direction.y > 0 ? 1 : -1;
//     Vec2i current_pixel = p0;

//     while (true)
//     {
//         image.set(current_pixel.x, current_pixel.y, to_tgacolor(color));
//         if (current_pixel.x == p1.x && current_pixel.y == p1.y)
//         {
//             break;
//         }

//         if (t_next.x <= t_next.y)
//         {
//             current_pixel.x += delta_pixel.x;
//             t_next.x += dt.x;
//         }
//         else // t_next.y > t_next.x
//         {
//             current_pixel.y += delta_pixel.y;
//             t_next.y += dt.y;
//         }
//     }
// }


/**
 * Vertices should be in DEVICE COORDINATES.
 * Bottom left of device is origin.
 */
void draw_line(Vec3f p0, Vec3f p1, TGAImage& image, Vec3f color)
{
    // Right now it ignores depth.
    // Deal with symmetry later.

    // March from t = 0 to t = 1
    // March from p0 to p1

    Vec2f direction (p1.x - p0.x, p1.y - p0.y);
    Vec2f delta_t (1.0f / std::abs(direction.x), 1.0f / std::abs(direction.y));
    Vec2i current_pixel (int(p0.x), int(p0.y));
    // Edge case: p0 is on 'boarder/edge' of device coordinates
    if (current_pixel.x == width) current_pixel.x--;
    if (current_pixel.y == height) current_pixel.y--;
    Vec2i delta_pixel (direction.x > 0 ? 1 : -1, direction.y > 0 ? 1 : -1); // Inquiry: direction.x might always be 1 due to symmetry constraint (left-to-right), figure out
    float delta_x = (delta_pixel.x == 1) ? (current_pixel.x + 1) - p0.x : p0.x - current_pixel.x;
    float delta_y = (delta_pixel.y == 1) ? (current_pixel.y + 1) - p0.y : p0.y - current_pixel.y;
    Vec2f next_t (delta_x / direction.x, delta_y / direction.y);


    while (next_t.x <= 1.0f + 1e-5f || next_t.y <= 1.0f + 1e-5f)
    {
        assert(current_pixel.x >= 0 && current_pixel.x < width && current_pixel.y >= 0 && current_pixel.y < height); // quick test, remove later
        image.set(current_pixel.x, current_pixel.y, to_tgacolor(color));

        if (next_t.x < next_t.y) // move 1 pixel horizontally
        {
            next_t.x += delta_t.x;
            current_pixel.x += delta_pixel.x;
        }
        else if (next_t.y > next_t.x) // move 1 pixel vertically
        {
            next_t.y += delta_t.y;
            current_pixel.y += delta_pixel.y;
        }
        else // move 1 pixel horizontally & move 1 pixel vertically (essentially move diagonally)
        {
            next_t.x += delta_t.x;
            current_pixel.x += delta_pixel.x;
            next_t.y += delta_t.y;
            current_pixel.y += delta_pixel.y;
        }
    }
}


/**
 * Triangle vertices should be in DEVICE COORDINATES.
 * Bottom left of device is origin.
 */
void draw_triangle(Vec3f a, Vec3f b, Vec3f c, TGAImage& image, Vec3f a_color, Vec3f b_color, Vec3f c_color, float** z_buffer)
{
    // QUESTION: are positive depth values in front of screen? if so should they be culled?

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
    max.x = std::min(max.x, width - 1);
    max.y = std::min(max.y, height - 1);

    Mat3x3f D (Vec3f(a.x, a.y, 1.0f), Vec3f(b.x, b.y, 1.0f), Vec3f(c.x, c.y, 1.0f));
    float determinant = D.determinant();

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
                float interpolated_depth = (a.z * alpha) + (b.z * beta) + (c.z * gamma);

                // Depth check
                if (z_buffer[row][col] < interpolated_depth)
                {
                    Vec3f interpolated_color = clampedVec3f((a_color * alpha) + (b_color * beta) + (c_color * gamma), 0.0f, 1.0f);
                    image.set(col, row, to_tgacolor(interpolated_color));
                    z_buffer[row][col] = interpolated_depth;
                }
            }
        }
    }
}