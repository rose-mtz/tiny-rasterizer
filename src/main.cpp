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


void draw_horizontal_line(Vec2i a, Vec2i b, Image& image, Vec3f color)
{
    assert(a.y == b.y);
    assert(a.x >= 0 && a.x < image.width && a.y >= 0 && a.y < image.height);
    assert(b.x >= 0 && b.x < image.width && b.y >= 0 && b.y < image.height);

    if (a.x > b.x) std::swap(a,b);

    Vec2i pixel (a);
    while (pixel.x <= b.x)
    {
        assert(pixel.x >= 0 && pixel.x < image.width && pixel.y >= 0 && pixel.y < image.height);
        image.image[pixel.y][pixel.x] = color;
        pixel.x++;
    }
}


void draw_triangle(Vec2i a, Vec2i b, Vec2i c, Image& image, Vec3f color)
{
    if (a.y > b.y) std::swap(a,b);
    if (a.y > c.y) std::swap(a,c);
    if (b.y > c.y) std::swap(b,c);

    if (a.y == b.y) // special case: 'flat base' triangle
    {
        // TODO
        assert(false);
    }
    else
    {
        // Set up march from a to b
        Vec2i direction_ab (b.x - a.x, b.y - a.y);
        Vec2f dt_ab (1.0f / std::abs(direction_ab.x), 1.0f / std::abs(direction_ab.y)); // note: (float) 1.0 / (int) 0 == infinity
        Vec2f t_next_ab (dt_ab); // the next values at which a 'x' or 'y intersection' occur
        Vec2i delta_pixel_ab;
        delta_pixel_ab.x = direction_ab.x > 0 ? 1 : -1;
        delta_pixel_ab.y = direction_ab.y > 0 ? 1 : -1;
        Vec2i current_pixel_ab = a;
        Vec2i target_pixel_ab = b;

        // Set up march from a to c
        Vec2i direction_ac (c.x - a.x, c.y - a.y);
        Vec2f dt_ac (1.0f / std::abs(direction_ac.x), 1.0f / std::abs(direction_ac.y)); // note: (float) 1.0 / (int) 0 == infinity
        Vec2f t_next_ac (dt_ac); // the next values at which a 'x' or 'y intersection' occur
        Vec2i delta_pixel_ac;
        delta_pixel_ac.x = direction_ac.x > 0 ? 1 : -1;
        delta_pixel_ac.y = direction_ac.y > 0 ? 1 : -1;
        Vec2i current_pixel_ac = a;
        Vec2i target_pixel_ac = c;

        // march a to b and a to c at the same time
        // synchronizing at y-values
        while (true)
        {
            draw_horizontal_line(current_pixel_ab, current_pixel_ac, image, color);
            if (
                (current_pixel_ab.x == target_pixel_ab.x && current_pixel_ab.y == target_pixel_ab.y) || 
                (current_pixel_ac.x == target_pixel_ac.x && current_pixel_ac.y == target_pixel_ac.y)
                )
            {
                break;
            }

            // march line ab to next y-value
            while (true)
            {
                if (t_next_ab.x <= t_next_ab.y)
                {
                    current_pixel_ab.x += delta_pixel_ab.x;
                    t_next_ab.x += dt_ab.x;
                }
                else
                {
                    current_pixel_ab.y += delta_pixel_ab.y;
                    t_next_ab.y += dt_ab.y;
                    break;
                }
            }

            // march line ac to next y-value
            while (true)
            {
                if (t_next_ac.x <= t_next_ac.y)
                {
                    current_pixel_ac.x += delta_pixel_ac.x;
                    t_next_ac.x += dt_ac.x;
                }
                else
                {
                    current_pixel_ac.y += delta_pixel_ac.y;
                    t_next_ac.y += dt_ac.y;
                    break;
                }
            }
        }

        // EDGE CASE: march ends for both at same time! --> V
        if (
            (current_pixel_ab.x == target_pixel_ab.x && current_pixel_ab.y == target_pixel_ab.y) &&
            (current_pixel_ac.x == target_pixel_ac.x && current_pixel_ac.y == target_pixel_ac.y)
            )
        {
            return;
        }
        else if (current_pixel_ab.x == target_pixel_ab.x && current_pixel_ab.y == target_pixel_ab.y)
        {
            // Set up march from b to c
            direction_ab  = Vec2i(c.x - b.x, c.y - b.y);
            dt_ab = Vec2f(1.0f / std::abs(direction_ab.x), 1.0f / std::abs(direction_ab.y)); // note: (float) 1.0 / (int) 0 == infinity
            t_next_ab = Vec2f(dt_ab); // the next values at which a 'x' or 'y intersection' occur
            delta_pixel_ab.x = direction_ab.x > 0 ? 1 : -1;
            delta_pixel_ab.y = direction_ab.y > 0 ? 1 : -1;
            current_pixel_ab = b;
            target_pixel_ab = c;
        }
        else
        {
            // Set up march from c to b
            direction_ac  = Vec2i(b.x - c.x, b.y - c.y);
            dt_ac = Vec2f(1.0f / std::abs(direction_ac.x), 1.0f / std::abs(direction_ac.y)); // note: (float) 1.0 / (int) 0 == infinity
            t_next_ac = Vec2f(dt_ac); // the next values at which a 'x' or 'y intersection' occur
            delta_pixel_ac.x = direction_ac.x > 0 ? 1 : -1;
            delta_pixel_ac.y = direction_ac.y > 0 ? 1 : -1;
            current_pixel_ac = c;
            target_pixel_ac = b;
        }

        // Set up march from (b to c OR c to b, depending on BEFORE shit)
        // I think easiest way to do this is just re-use old variables

        // continue marching
        while (true)
        {
            // In-Efficient: 'middle/perpendicular' line getting drawn twice
            draw_horizontal_line(current_pixel_ab, current_pixel_ac, image, color);
            if (
                (current_pixel_ab.x == target_pixel_ab.x && current_pixel_ab.y == target_pixel_ab.y) || 
                (current_pixel_ac.x == target_pixel_ac.x && current_pixel_ac.y == target_pixel_ac.y)
                )
            {
                break;
            }

            // march line ab to next y-value
            while (true)
            {
                if (t_next_ab.x <= t_next_ab.y)
                {
                    current_pixel_ab.x += delta_pixel_ab.x;
                    t_next_ab.x += dt_ab.x;
                }
                else
                {
                    current_pixel_ab.y += delta_pixel_ab.y;
                    t_next_ab.y += dt_ab.y;
                    break;
                }
            }

            // march line ac to next y-value
            while (true)
            {
                if (t_next_ac.x <= t_next_ac.y)
                {
                    current_pixel_ac.x += delta_pixel_ac.x;
                    t_next_ac.x += dt_ac.x;
                }
                else
                {
                    current_pixel_ac.y += delta_pixel_ac.y;
                    t_next_ac.y += dt_ac.y;
                    break;
                }
            }
        }
    }
}