#include <iostream>
#include "Image.h"


void draw_line(int x0, int y0, int x1, int y1, Image& image, Vec3f color);


int main()
{   
    Image image;
    init_image(image, 500, 500);

    Vec3f white (1.0f, 1.0f, 1.0f);
    Vec3f red (1.0f, 0.0f, 0.0f);
    Vec3f green (0.0f, 1.0f, 0.0f);    
    draw_line(0, 0, 0, 199, image, red);
    draw_line(99, 20, 0, 0, image, white);
    draw_line(145, 9, 12, 99, image, red);
    draw_line(100, 100, 499, 499, image, red);
    draw_line(25, 250, 499, 250, image, red);
    draw_line(10, 10, 10, 499, image, red);
    

    save_image("./image.bmp", image);

    return 0;
}


// pixel-coordinates
void draw_line(int x0, int y0, int x1, int y1, Image& image, Vec3f color)
{
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

        if (pix_x == x1 && pix_y == y1)
        {
            break;
        }
    }
}