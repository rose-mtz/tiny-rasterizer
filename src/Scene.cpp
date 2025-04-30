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


SkyBox* parse_skybox(std::ifstream& in)
{
    // Defaults
    TGAImage* front = nullptr;
    TGAImage* right = nullptr;
    TGAImage* back = nullptr;
    TGAImage* left = nullptr;
    TGAImage* top = nullptr;
    TGAImage* bottom = nullptr;
    float yaw = 0.0f;
    float pitch = 0.0f;
    float roll = 0.0f;
    float fov = 90.0f;

    std::string line; std::getline(in, line);
    while (line.size() != 0) // Empty line indicates end of attributes
    {
        std::istringstream iss(line.c_str());
        std::string attribute; iss >> attribute;

        if (attribute == "front")
        {
            std::string filename; assert(iss >> filename);
            front = new TGAImage();
            front->read_tga_file(filename.c_str());
            front->flip_vertically();
        }
        else if (attribute == "back")
        {
            std::string filename; assert(iss >> filename);
            back = new TGAImage();
            back->read_tga_file(filename.c_str());
            back->flip_vertically();
        }
        else if (attribute == "right")
        {
            std::string filename; assert(iss >> filename);
            right = new TGAImage();
            right->read_tga_file(filename.c_str());
            right->flip_vertically();
        }
        else if (attribute == "left")
        {
            std::string filename; assert(iss >> filename);
            left = new TGAImage();
            left->read_tga_file(filename.c_str());
            left->flip_vertically();
        }
        else if (attribute == "top")
        {
            std::string filename; assert(iss >> filename);
            top = new TGAImage();
            top->read_tga_file(filename.c_str());
            top->flip_vertically();
        }
        else if (attribute == "bottom")
        {
            std::string filename; assert(iss >> filename);
            bottom = new TGAImage();
            bottom->read_tga_file(filename.c_str());
            bottom->flip_vertically();
        }
        else if (attribute == "rotations")
        {
            // rotations <yaw> <pitch> <roll> 
            // degrees 
            assert(iss >> yaw);
            assert(iss >> pitch);
            assert(iss >> roll);
        }
        else if (attribute == "fov")
        {
            assert(iss >> fov);
        }
        else
        {
            std::cout << "Error:: unkown Image_Metadata attribute " << attribute << '\n';
            assert(false);
        }

        if (in.eof()) // was last line of file
        {
            break;
        }
        std::getline(in, line); // get next line
    }

    SkyBox* skybox = new SkyBox();
    skybox->textures[FRONT] = front;
    skybox->textures[BACK] = back;
    skybox->textures[RIGHT] = right;
    skybox->textures[LEFT] = left;
    skybox->textures[TOP] = top;
    skybox->textures[BOTTOM] = bottom;
    skybox->fov = radians(fov);
    skybox->yaw = radians(yaw);
    skybox->pitch = radians(pitch);
    skybox->roll = radians(roll);

    return skybox;
}


ImageMetadata* parse_image_metadata(std::ifstream& in)
{
    // Defaults
    Vec3f background_color = Vec3f(0.0f);
    Vec2i aspect_ratio = Vec2i(1, 1);
    int width_pixels = 500;
    int supersample_factor = 1;
    std::string save_location = "./output.tga";

    std::string line; std::getline(in, line);
    while (line.size() != 0) // Empty line indicates end of attributes
    {
        std::istringstream iss(line.c_str());
        std::string attribute; iss >> attribute;

        if (attribute == "background_color")
        {
            background_color = parse_vec3f(iss);
        }
        else if (attribute == "aspect_ratio")
        {
            assert(iss >> aspect_ratio.x);
            assert(iss >> aspect_ratio.y);
        }
        else if (attribute == "width")
        {
            assert(iss >> width_pixels);
        }
        else if (attribute == "supersample_factor")
        {
            assert(iss >> supersample_factor);
        }
        else if (attribute == "save_location")
        {
            assert(iss >> save_location);
        }
        else
        {
            std::cout << "Error:: unkown Image_Metadata attribute " << attribute << '\n';
            assert(false);
        }

        if (in.eof()) // was last line of file
        {
            break;
        }
        std::getline(in, line); // get next line
    }

    ImageMetadata* metadata = new ImageMetadata();
    metadata->background_color = background_color;
    metadata->aspect_ratio = aspect_ratio;
    metadata->width_pixels = width_pixels;
    metadata->supersample_factor = supersample_factor;
    metadata->save_location = save_location;

    return metadata;
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
    std::string filename; 
    bool clockwise_winding = true;

    std::string line; std::getline(in, line);
    while (line.size() != 0) // blank line indicates end of attributes
    {
        std::istringstream iss(line.c_str());
        std::string attribute; iss >> attribute;

        if (attribute == "filename")
        {
            assert(iss >> filename);
        }
        else if (attribute == "winding")
        {
            std::string winding; iss >> winding;

            if (winding == "clockwise")
            {
                clockwise_winding = true;
            }
            else if (winding == "counter_clockwise")
            {
                clockwise_winding = false;
            }
            else
            {
                assert(false);
            }
        }

        if (in.eof()) // was last line of file
        {
            break;
        }
        std::getline(in, line); // get next line   
    }

    Model* model = new Model(filename.c_str());
    model->clockwise_winding = clockwise_winding;

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
    SHADING_TYPE shading = SHADING_TYPE::FLAT;
    bool wireframe = false;
    bool fill = false;
    FILL_MODE fill_mode = FILL_MODE::TEXTURE;
    float yaw = 0;
    float pitch = 0;
    float roll = 0;

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
            std::string shading_type_str; iss >> shading_type_str;

            if (shading_type_str == "flat")
            {
                shading = SHADING_TYPE::FLAT;
            }
            else if (shading_type_str == "gouraud")
            {
                shading = SHADING_TYPE::GOURAUD;
            }
            else if (shading_type_str == "phong")
            {
                shading = SHADING_TYPE::PHONG;
            }
            else if (shading_type_str == "none")
            {
                shading = SHADING_TYPE::NONE;
            }
            else 
            {
                assert(false);
            }
        }
        else if (attribute == "modes")
        {
            // Flags turn on mode, else default is off
            std::string mode;
            while (iss >> mode)
            {
                if (mode == "fill")
                {
                    fill = true;
                }
                else if (mode == "wireframe")
                {
                    wireframe = true;
                }
                else if (mode == "colored_face_normals")
                {
                    fill_mode = FILL_MODE::COLORED_FACE_NORMALS;
                }
                else if (mode == "colored_vertex_normals")
                {
                    fill_mode = FILL_MODE::COLORED_VERTEX_NORMALS;
                }
                else if (mode == "texture")
                {
                    fill_mode = FILL_MODE::TEXTURE;
                }
                else if (mode == "vertex_colors")
                {
                    fill_mode = FILL_MODE::VERTEX_COLORS;
                }
                else
                {
                    std::cout << "Error:: unkown mode value " << mode << '\n';
                    assert(false);
                }
            }
        }
        else if (attribute == "rotations")
        {
            // rotations <yaw> <pitch> <roll> 
            // degrees 
            assert(iss >> yaw);
            assert(iss >> pitch);
            assert(iss >> roll);
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
    obj->wireframe = wireframe;
    obj->fill = fill;
    obj->fill_mode = fill_mode;
    obj->yaw = radians(yaw);
    obj->pitch = radians(pitch);
    obj->roll = radians(roll);

    return obj;
}


Camera* parse_camera(std::ifstream& in)
{
    std::string type = "orthographic";
    Vec3f pos = Vec3f(0.0f, 0.0f, 10.0f);
    Vec3f look_at = Vec3f(0.0f, 0.0f, 0.0f);
    Vec3f up = Vec3f(0.0f, 1.0f, 0.0f);
    float zoom = 1.0f;
    float fov = radians(45);

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
        else if (attribute == "fov")
        {
            assert(iss >> fov);
            fov = radians(fov);
        }
        else if (attribute == "up")
        {
            up = parse_vec3f(iss);
        }
        else
        {
            std::cout << "Error:: unkown Camera attribute " << attribute << '\n';
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
    cam->up = up.normalize();
    cam->zoom = zoom;
    cam->fov = fov;

    assert(!((up ^ (pos - look_at)) == Vec3f(0.0f, 0.0f, 0.0f))); // invalid up

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
            std::cout << "Error:: unkown Light attribute " << attribute << '\n';
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
        else if (type == "Image_Metadata")
        {
            this->metadata = parse_image_metadata(in);
        }
        else if (type == "Skybox")
        {
            this->skybox = parse_skybox(in);
        }
        else if (type == "#") // comment
        {
            continue;
        }
        else if (line.size() != 0)
        {
            std::cout << "Error:: unkown type " << line << '\n';
            assert(false); // Error: unkown type
        }
    }
}


Scene::~Scene() {}