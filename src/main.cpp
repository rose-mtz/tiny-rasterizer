#include <iostream>
#include "Image.h"

int main()
{   
    Image image;
    init_image(image, 10, 10);

    Vec3f white (1.0f, 1.0f, 1.0f);
    Vec3f red (1.0f, 0.0f, 0.0f);
    Vec3f green (0.0f, 1.0f, 0.0f);    
    image.image[0][0] = red;
    image.image[9][0] = green;
    image.image[9][9] = white;

    save_image("./image.bmp", image);

    return 0;
}