#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include "Image.h"
#include "Mesh.h"


void draw_line(Vec2i p0, Vec2i p1, Image& image, Vec3f color);
void draw_triangle(Vec2i a, Vec2i b, Vec2i c, Image& image, Vec3f color);

int main()
{   
    // Mesh mesh = load_mesh("./obj/african_head.obj");

    Image image;
    init_image(image, 100, 100);

    Vec3f white (1.0f, 1.0f, 1.0f);
    Vec3f red   (1.0f, 0.0f, 0.0f);

    // for (int i = 0; i < mesh.faces.size(); i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //     {
    //         Vec3f v0 = mesh.vertices[mesh.faces[i][j].x];
    //         Vec3f v1 = mesh.vertices[mesh.faces[i][(j + 1) % 3].x];

    //         int x0 = (v0.x + 1.0f) * ((image.width  - 1)/ 2.0f);
    //         int y0 = (v0.y + 1.0f) * ((image.height - 1) / 2.0f);
    //         int x1 = (v1.x + 1.0f) * ((image.width  - 1)/ 2.0f);
    //         int y1 = (v1.y + 1.0f) * ((image.height - 1) / 2.0f);

    //         draw_line(Vec2i(x0, y0), Vec2i(x1, y1), image, white);
    //     }
    // }

    Vec2i a(0,0), b(82,34), c(24,99);
    draw_triangle(a, b, c, image, white);
    draw_line(a, b, image, red);
    draw_line(a, c, image, red);
    draw_line(b, c, image, red);

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


void draw_triangle(Vec2i a, Vec2i b, Vec2i c, Image& image, Vec3f color)
{
    // Lets do the BB later
        // BoundingBox get_clipped_bounding_box(a,b,c)
        // BoundingBox { Vec2i bottom_left, Vec2i top_right }
    // No AA, just dead center of pixel in or out
    // lets also just use lower left instead of 'center'

    Mat3x3 D (Vec3f(a, 1), Vec3f(b, 1), Vec3f(c, 1));
    float determinant_of_d = determinant(D);

    for (int row = 0; row < image.height; row++)
    {
        for (int col = 0; col < image.width; col++)
        {
            Mat3x3 D_x (Vec3f(col, row, 1), Vec3f(b, 1), Vec3f(c, 1));
            Mat3x3 D_y (Vec3f(a, 1), Vec3f(col, row, 1), Vec3f(c, 1));
            Mat3x3 D_z (Vec3f(a, 1), Vec3f(b, 1), Vec3f(col, row, 1));

            float determinant_of_dx = determinant(D_x);
            float determinant_of_dy = determinant(D_y);
            float determinant_of_dz = determinant(D_z);

            float alpha = determinant_of_dx / determinant_of_d;
            float beta = determinant_of_dy / determinant_of_d;
            float gamma = determinant_of_dz / determinant_of_d;

            // assert(alpha + beta + gamma == 1.0f); // FAILS sometimes!

            if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f)
            {
                image.image[row][col] = color;
            }
        }
    }
}