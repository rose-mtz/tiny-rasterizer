#ifndef __SCENE__
#define __SCENE__

#include "Math.h"
#include "model.h"
#include "tgaimage.h"
#include <string>
#include <vector>


struct ImageMetadata
{
    Vec3f background_color;
    Vec2i aspect_ratio;
    int width_pixels;
    int supersample_factor;
};


struct Material 
{
    float k_d; // diffuse
    float k_s; // specular
    float k_a; // ambient
    float shininess; // q
};


struct Object3D
{
    Model* model;
    TGAImage* texture;
    Material* mat;
    std::string shading;
    Vec3f pos;
    float scale;

    // Modes
    bool colored_triangle_normals_mode;
    bool colored_vertex_normals_mode;
    bool wireframe_mode;
    bool fill_mode;
};


struct Camera 
{
    std::string type; // refactor out later
    Vec3f pos;
    Vec3f look_at;
    float zoom;
};


struct Light
{
    std::string type; // refactor out later
    Vec3f color;

    // Directional lights

    Vec3f direction; // direction light moves towards

    // Point lights

    Vec3f pos;
    float intensity;
};


class Scene
{
public:
    std::vector<Light*> lights;
    std::vector<Object3D*> objects;
    Camera* camera;
    ImageMetadata* metadata;

    Scene(const char* filename);
    ~Scene();
};


#endif