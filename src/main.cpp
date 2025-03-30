#include <iostream>
#include <cassert>
#include "Image.h"
#include "Mesh.h"


void draw_line(int x0, int y0, int x1, int y1, Image& image, Vec3f color);


int main()
{   
    Mesh mesh = load_mesh("./obj/african_head.obj");

    Image image;
    init_image(image, 1000, 1000);

    Vec3f white (1.0f, 1.0f, 1.0f);

    for (int i = 0; i < mesh.faces.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Vec3f v0 = mesh.vertices[mesh.faces[i][j].x];
            Vec3f v1 = mesh.vertices[mesh.faces[i][(j + 1) % 3].x];

            int x0 = (v0.x + 1.0f) * ((image.width  - 1)/ 2.0f);
            int y0 = (v0.y + 1.0f) * ((image.height - 1) / 2.0f);
            int x1 = (v1.x + 1.0f) * ((image.width  - 1)/ 2.0f);
            int y1 = (v1.y + 1.0f) * ((image.height - 1) / 2.0f);

            draw_line(x0, y0, x1, y1, image, white);
        }
    }


    save_image("./image.bmp", image);

    return 0;
}


// pixel-coordinates
void draw_line(int x0, int y0, int x1, int y1, Image& image, Vec3f color)
{
    assert(x0 >= 0 && x0 < image.width);
    assert(y0 >= 0 && y0 < image.height);
    assert(x1 >= 0 && x1 < image.width);
    assert(y1 >= 0 && y1 < image.height);

    // left-to-right (for symmetry)
    if (x0 > x1)
    {
        int temp_x = x0;
        x0 = x1;
        x1 = temp_x;
        
        int temp_y = y0;
        y0 = y1;
        y1 = temp_y;
    }

    int dir_x = x1 - x0;
    int dir_y = y1 - y0;
    // Note: 1.0f / 0 --> inf
    float dt_x = 1.0f / dir_x;
    float dt_y = 1.0f / dir_y;
    // Note: Ray(t = 0) = (x0, y0), and Ray(t = 1) = (x1, y1), so dt_x, dt_y can't be negative where going from 0 to 1
    dt_x = dt_x < 0.0f ? dt_x * -1.0f : dt_x;
    dt_y = dt_y < 0.0f ? dt_y * -1.0f : dt_y;
    float t_x = dt_x;
    float t_y = dt_y;

    int dPix_x = x0 < x1 ? 1 : -1;
    int dPix_y = y0 < y1 ? 1 : -1;

    int pix_x = x0;
    int pix_y = y0;
    image.image[pix_y][pix_x] = color;
    while (true)
    {
        if (pix_x == x1 && pix_y == y1)
        {
            break;
        }

        if (t_x < t_y)
        {
            pix_x += dPix_x;
            t_x += dt_x;
            image.image[pix_y][pix_x] = color;
        }
        else if (t_y < t_x)
        {
            pix_y += dPix_y;
            t_y += dt_y;
            image.image[pix_y][pix_x] = color;   
        }
        else
        {
            pix_x += dPix_x;
            t_x += dt_x;
            image.image[pix_y][pix_x] = color;
            pix_y += dPix_y;
            t_y += dt_y;
            image.image[pix_y][pix_x] = color;  
        }
    }
}