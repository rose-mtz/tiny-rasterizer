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
    std::string save_location;
};


struct Material 
{
    float k_d; // diffuse
    float k_s; // specular
    float k_a; // ambient
    float shininess; // q
};

// Fill mode for filling in faces of object
enum class FILL_MODE { TEXTURE, VERTEX_COLORS, COLORED_FACE_NORMALS, COLORED_VERTEX_NORMALS };
enum class SHADING_TYPE { NONE, FLAT, GOURAUD, PHONG };

struct Object3D
{
    Model* model;
    TGAImage* texture;
    Material* mat;
    SHADING_TYPE shading;

    // Transformations
    Vec3f pos;
    float scale; // TODO: make it vec3
    float yaw;
    float pitch;
    float roll;

    // Modes
    bool wireframe;
    bool fill;
    FILL_MODE fill_mode;
};


struct Camera 
{
    std::string type; // refactor out later
    Vec3f pos;
    Vec3f look_at;
    Vec3f up;

    // orthographic
    float zoom;

    // perspective
    float fov;
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


enum FACE { FRONT = 0, RIGHT, BACK, LEFT, TOP, BOTTOM }; // Skybox faces

struct SkyBox
{
    TGAImage* textures[6];
    float fov;
    float yaw;
    float pitch;
    float roll;   
};


class Scene
{
public:
    std::vector<Light*> lights;
    std::vector<Object3D*> objects;
    Camera* camera;
    ImageMetadata* metadata;
    SkyBox* skybox = nullptr;

    Scene(const char* filename);
    ~Scene();
};


#endif