#ifndef __SCENE__
#define __SCENE__

#include "Math.h"
#include "model.h"
#include "tgaimage.h"
#include <string>
#include <vector>


/**
 * I really should just make this shit into a class.
 * And subclasses! KISS, later.
 */


struct Object3D
{
    Model* model;
    TGAImage* texture;
    bool gouraud_shading;
    Vec3f pos;
    float scale;
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
    Vec3f direction;
    Vec3f color;
    Vec3f pos;
};


class Scene
{
// private:
public:
    std::vector<Light*> lights;
    std::vector<Object3D*> objects;
    Camera* camera;

    Scene(const char* filename);
    ~Scene();
};


#endif