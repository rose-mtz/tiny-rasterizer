struct Vec3f
{
    float x, y, z;

    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
};


// [0,0] should be bottom left
// [height - 1, width - 1] should be top right
struct Image
{
    Vec3f** image = nullptr;
    int width;
    int height;
};


void save_image(const char* filename, Image& image);
void init_image(Image& image, int width, int height);
// void destruct_image(Image& image);