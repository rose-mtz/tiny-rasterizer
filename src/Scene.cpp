#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include "Scene.h"


std::vector<Model*> models;
std::vector<TGAImage*> textures;
std::vector<Material*> materials;


Vec3f parse_vec3f(std::istringstream& iss)
{
    Vec3f v;
    assert(iss >> v.x);
    assert(iss >> v.y);
    assert(iss >> v.z);
    return v;
}


bool parse_bool(std::istringstream& iss)
{
    std::string boolean; iss >> boolean;
    return boolean == "true" ? true : false;
}


Material* parse_material(std::ifstream& in)
{
    float k_d = 1.0f;
    float k_s = 1.0f;
    float k_a = 1.0f;
    int shininess = 10;

    std::string line; std::getline(in, line);
    while (line.size() != 0) // blank line indicates end of attributes
    {
        std::istringstream iss(line.c_str());
        std::string attribute; iss >> attribute;

        if (attribute == "k_d")
        {
            assert(iss >> k_d);
        }
        else if (attribute == "k_s")
        {
            assert(iss >> k_s);   
        }
        else if (attribute == "k_a")
        {
            assert(iss >> k_a);
        }
        else if (attribute == "shininess")
        {
            assert(iss >> shininess);
        }

        if (in.eof()) // was last line of file
        {
            break;
        }
        std::getline(in, line); // get next line
    }

    Material* mat = new Material();
    mat->k_d = k_d;
    mat->k_s = k_s;
    mat->k_a = k_a;
    mat->shininess = shininess;

    return mat;
}


Model* parse_model(std::ifstream& in)
{
    std::string line; std::getline(in, line);
    std::istringstream iss(line.c_str());
    std::string attribute; assert(iss >> attribute);
    std::string filename; assert(iss >> filename);

    Model* model = new Model(filename.c_str());

    return model;
}


TGAImage* parse_texture(std::ifstream& in)
{
    std::string line; std::getline(in, line);
    std::istringstream iss(line.c_str());
    std::string attribute; assert(iss >> attribute);
    std::string filename; assert(iss >> filename);

    TGAImage* texture = new TGAImage();
    texture->read_tga_file(filename.c_str());
    texture->flip_vertically(); // can later allow user to specify if they want it flipped in txt file

    return texture;
}


Object3D* parse_object3d(std::ifstream& in)
{
    // Default arguments
    Model* model = nullptr;
    TGAImage* texture = nullptr;
    Material* mat = nullptr;
    Vec3f position = Vec3f(0.0f, 0.0f, 0.0f);
    float scale = 1.0f;
    std::string shading = "flat";

    std::string line; std::getline(in, line);
    while (line.size() != 0) // blank line indicates end of attributes
    {
        std::istringstream iss(line.c_str());
        std::string attribute; iss >> attribute;

        if (attribute == "model")
        {
            int index; assert(iss >> index);
            assert(models.size() > index);
            model = models[index];
        }
        else if (attribute == "texture")
        {
            int index; assert(iss >> index);
            assert(textures.size() > index);
            texture = textures[index];
        }
        else if (attribute == "material")
        {
            int index; assert(iss >> index);
            assert(materials.size() > index);
            mat = materials[index];
        }
        else if (attribute == "position")
        {
            position = parse_vec3f(iss);
        }
        else if (attribute == "scale")
        {
            assert(iss >> scale);   
        }
        else if (attribute == "shading")
        {
            assert(iss >> shading);
            assert(shading == "flat" || shading == "gouraud" || shading == "phong");
        }
        else
        {
            std::cout << "Error:: unkown Object3D attribute " << attribute << '\n';
            assert(false);
        }

        if (in.eof()) // was last line of file
        {
            break;
        }
        std::getline(in, line); // get next line
    }

    Object3D* obj = new Object3D();
    obj->model = model;
    obj->texture = texture;
    obj->mat = mat;
    obj->shading = shading;
    obj->pos = position;
    obj->scale = scale;

    return obj;
}


Camera* parse_camera(std::ifstream& in)
{
    std::string type = "orthographic";
    Vec3f pos = Vec3f(0.0f, 0.0f, 10.0f);
    Vec3f look_at = Vec3f(0.0f, 0.0f, 0.0f);
    float zoom = 1.0f;

    std::string line; std::getline(in, line);
    while (line.size() != 0) // blank line indicates end of attributes
    {
        std::istringstream iss(line.c_str());
        std::string attribute; iss >> attribute;

        if (attribute == "type")
        {
            assert(iss >> type);
            assert(type == "perspective" || type == "orthographic");
        }
        else if (attribute == "position")
        {
            pos = parse_vec3f(iss);
        }
        else if (attribute == "lookAt")
        {
            look_at = parse_vec3f(iss);
        }
        else if (attribute == "zoom")
        {
            assert(iss >> zoom);   
        }
        else
        {
            std::cout << "Error:: unkown Object3D attribute " << attribute << '\n';
            assert(false);
        }

        if (in.eof()) // was last line of file
        {
            break;
        }
        std::getline(in, line); // get next line
    }

    Camera* cam = new Camera();
    cam->type = type;
    cam->pos = pos;
    cam->look_at = look_at;
    cam->zoom = zoom;

    return cam;
}


Light* parse_light(std::ifstream& in)
{
    std::string type = "directional";
    Vec3f direction = Vec3f(0.0f, 0.0f, -1.0f);
    Vec3f color = Vec3f(1.0f, 1.0f, 1.0f);
    Vec3f pos = Vec3f(10.0f, 10.0f, 10.0f);
    float intensity = 1.0f;

    std::string line; std::getline(in, line);
    while (line.size() != 0) // blank line indicates end of attributes
    {
        std::istringstream iss(line.c_str());
        std::string attribute; iss >> attribute;

        if (attribute == "type")
        {
            assert(iss >> type);
            assert(type == "directional" || type == "point");
        }
        else if (attribute == "position")
        {
            pos = parse_vec3f(iss);
        }
        else if (attribute == "direction")
        {
            direction = parse_vec3f(iss).normalize();
        }
        else if (attribute == "color")
        {
            color = parse_vec3f(iss);
        }
        else if (attribute == "intensity")
        {
            assert(iss >> intensity);
        }
        else
        {
            std::cout << "Error:: unkown Object3D attribute " << attribute << '\n';
            assert(false);
        }

        if (in.eof()) // was last line of file
        {
            break;
        }
        std::getline(in, line); // get next line
    }

    Light* light = new Light();
    light->type = type;
    light->direction = direction;
    light->color = color;
    light->pos = pos;
    light->intensity = intensity;

    return light;
}


Scene::Scene(const char* filename)
{
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail())
    {
        std::cerr << "Error:: could not open scene file " << filename << '\n';
    }
    std::string line;
    while (!in.eof())
    {
        std::getline(in, line);
        std::istringstream iss(line.c_str());

        std::string type; iss >> type;
        if (type == "Model")
        {
            models.push_back(parse_model(in));
        }
        else if (type == "Texture")
        {
            textures.push_back(parse_texture(in));
        }
        else if (type == "Material")
        {
            materials.push_back(parse_material(in));
        }
        else if (type == "Camera")
        {
            this->camera = parse_camera(in);
        }
        else if (type == "Object3D")
        {
            this->objects.push_back(parse_object3d(in));
        }
        else if (type == "Light")
        {
            this->lights.push_back(parse_light(in));
        }
        else if (line.size() != 0)
        {
            std::cout << "Error:: unkown object type " << line << '\n';
            assert(false); // Error: unkown type
        }
    }
}


Scene::~Scene() {}