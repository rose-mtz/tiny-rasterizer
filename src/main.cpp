#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <limits>
#include <chrono>

#include "tgaimage.h"
#include "model.h"
#include "Scene.h"
#include "Rasterize.h"


const Vec3f WHITE (1.0f, 1.0f, 1.0f);
const Vec3f RED   (1.0f, 0.0f, 0.0f);
const Vec3f GREEN (0.0f, 1.0f, 0.0f);
const Vec3f BLUE  (0.0f, 0.0f, 1.0f);
const Vec3f BLACK (0.0f, 0.0f, 0.0f);

int supersample_factor;
int device_width;
int device_height;
int supersample_device_width;
int supersample_device_height;
float aspect_ratio;

float** z_buffer = nullptr;
Vec3f** color_buffer = nullptr;


struct Vertex
{
    Vec3f pos_device;
    Vec3f pos_camera;
    Vec3f norm_camera;
    Vec2f uv;
    Vec3f shading;
    Vec3f color;
    Vec3f norm_world;
    Vec3f pos_world;
};

struct VirtualScreen
{
    float top, right, bottom, left, near, far;
};


Vec3f    get_colored_normal(Vec3f normal);
TGAColor to_tga_color(Vec3f color);
Vec3f    to_vec3f_color(TGAColor color);
Vec3f    sample_texture(Vec2f uv, TGAImage* texture);
float    calculate_shading(Vec3f p, Vec3f n, Material mat, Light* light, Vec3f camera);

void draw_line(Vec3f v0, Vec3f v1, float thickness, Vec3f color, bool depth_check = true);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, FILL_MODE fill, TGAImage* texture = nullptr);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, std::vector<Light*> lights, Material* mat, FILL_MODE fill, TGAImage* texture = nullptr);
void draw_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3, FILL_MODE fill, TGAImage* texture = nullptr);

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    if (argc != 2)
    {
        std::cout << "Error: expected only 1 command line argument (scene file). \n";
        return 0;
    }
    Scene scene (argv[1]);

    // Downsampled device dimensions
    device_width  = scene.metadata->width_pixels;
    device_height = (scene.metadata->aspect_ratio.y * device_width) / scene.metadata->aspect_ratio.x;
    aspect_ratio  = ((float) device_width) / ((float) device_height); // the one true aspect ratio

    // Supersampled device dimensions
    supersample_factor        = scene.metadata->supersample_factor;
    supersample_device_width  = device_width * supersample_factor;
    supersample_device_height = device_height * supersample_factor;

    // Virtual screen dimensions
    VirtualScreen virtual_screen;
    if (scene.camera->type == "perspective")
    {
        virtual_screen.top    = 0.5f;
        virtual_screen.right  = aspect_ratio/2.0f;
        virtual_screen.bottom = -0.5f;
        virtual_screen.left   = -aspect_ratio/2.0f;
        virtual_screen.near   = (aspect_ratio/2.0f)/tan(scene.camera->fov/2.0f);
        virtual_screen.far    = 100.0f;
    }
    else // orthographic
    {
        virtual_screen.top    = 0.5f * scene.camera->zoom;
        virtual_screen.right  = aspect_ratio/2.0f * scene.camera->zoom;
        virtual_screen.bottom = -0.5f * scene.camera->zoom;
        virtual_screen.left   = -aspect_ratio/2.0f * scene.camera->zoom;
        virtual_screen.near   = 0.0f;
        virtual_screen.far    = 100.0f;
    }

    // Initialize depth, and color buffer
    color_buffer = new Vec3f*[supersample_device_height];
    for (int i = 0; i < supersample_device_height; i++)
    {
        color_buffer[i] = new Vec3f[supersample_device_width];
        for (int j = 0; j < supersample_device_width; j++)
        {
            color_buffer[i][j] = scene.metadata->background_color;
        }
    }

    z_buffer = new float*[supersample_device_height];
    for (int i = 0; i < supersample_device_height; i++)
    {
        z_buffer[i] = new float[supersample_device_width];
        for (int j = 0; j < supersample_device_width; j++)
        {
            z_buffer[i][j] = std::numeric_limits<float>::max(); // positive infinity is far away
        }
    }

    // Camera
    Mat4x4f camera = look_at(scene.camera->pos, scene.camera->look_at, scene.camera->up);

    // Projection
    Mat4x4f projection;
    if (scene.camera->type == "perspective") projection = perspective(virtual_screen.top, virtual_screen.right, virtual_screen.bottom, virtual_screen.left, virtual_screen.near, virtual_screen.far);
    else projection = orthographic(virtual_screen.top, virtual_screen.right, virtual_screen.bottom, virtual_screen.left, virtual_screen.near, virtual_screen.far);

    // Device (from NDC)
    float ndc_width  = 2.0f;
    float ndc_height = 2.0f;
    Vec3f device_center    = Vec3f(supersample_device_width / 2.0f, supersample_device_height / 2.0f, 0.0f);
    Vec3f ndc_scale_factor = Vec3f((supersample_device_width / ndc_width) / 2.0f, (supersample_device_height / ndc_height) / 2.0f, 1.0f);
    Mat4x4f device = translation(device_center) * scale(ndc_scale_factor);

    // Transform lights

    Mat3x3f camera_inv_trans = camera.truncated().inv().transposed();
    for (int l = 0; l < scene.lights.size(); l++)
    {
        // world to camera
        Light* a_light = scene.lights[l];
        a_light->pos = (camera * Vec4f(a_light->pos, 1.0f)).xyz();
        a_light->direction = (camera_inv_trans * a_light->direction).normalize();
    }

    // Render objects

    for (int i = 0; i < scene.objects.size(); i++)
    {
        Object3D* obj = scene.objects[i];

        // Set up transformation matrices of this object

        // From local to world
        Mat4x4f world = translation(obj->pos) * rotation_y(obj->yaw) * rotation_x(obj->pitch) * rotation_z(obj->roll) * scale(Vec3f(obj->scale));
        Mat3x3f world_inv_trans = (scale(Vec3f(1.0f/obj->scale)) * rotation_z(-obj->roll) * rotation_x(-obj->pitch) * rotation_y(-obj->yaw)).truncated().transposed();
        // From local to world to camera
        Mat4x4f camera_world = camera * world;
        Mat3x3f camera_world_inv_trans = camera_world.truncated().inv().transposed();

        // Render mesh faces of object

        for (int f = 0; f < obj->model->nfaces(); f++)
        {
            std::vector<int> face = obj->model->face(f);
            // three indicies per vertex
            std::vector<Vertex> vertices (face.size() / 3);

            // Do transformations on each vertex
            for (int v = 0; v < vertices.size(); v++)
            {
                int pos_index   = (v * 3);
                int uv_index    = (v * 3) + 1;
                int color_index = (v * 3) + 1;
                int norm_index  = (v * 3) + 2;

                vertices[v].uv    = obj->model->uv(face[uv_index]);
                vertices[v].color = obj->model->color(face[color_index]);

                // local to world
                vertices[v].pos_world  = (world * Vec4f(obj->model->vert(face[pos_index]), 1.0f)).xyz();
                vertices[v].norm_world = (world_inv_trans * obj->model->norm(face[norm_index])).normalize();
                // local to world to camera
                vertices[v].pos_camera  = (camera_world * Vec4f(obj->model->vert(face[pos_index]), 1.0f)).xyz();
                vertices[v].norm_camera = (camera_world_inv_trans * obj->model->norm(face[norm_index])).normalize();
                // project to NDC
                Vec3f vertex_ndc = (projection * Vec4f(vertices[v].pos_camera, 1.0f)).homogenize();
                // NDC to device screen
                vertices[v].pos_device = (device * Vec4f(vertex_ndc, 1.0f)).xyz();
            }


            // Back face culling calculation
            Vec3f triangle_normal_device;
            if (obj->model->clockwise_winding)
            {
                triangle_normal_device = triangle_normal(Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, 0.0f), Vec3f(vertices[vertices.size() - 1].pos_device.x, vertices[vertices.size() - 1].pos_device.y, 0.0f), Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, 0.0f));
            }
            else
            {
                triangle_normal_device = triangle_normal(Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, 0.0f), Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, 0.0f), Vec3f(vertices[vertices.size() - 1].pos_device.x, vertices[vertices.size() - 1].pos_device.y, 0.0f));
            }
            bool back_facing = triangle_normal_device == Vec3f(0.0f, 0.0f, 1.0f);


            // Do shading for none, flat, or gouraud
            if (obj->shading == SHADING_TYPE::NONE && !back_facing)
            {
                for (int v = 0; v < vertices.size(); v++)
                { 
                    vertices[v].shading = Vec3f(1.0f, 1.0f, 1.0f);
                }
            }
            else if (obj->shading == SHADING_TYPE::FLAT && !back_facing)
            {
                Vec3f face_normal_camera = triangle_normal(vertices[0].pos_camera, vertices[1].pos_camera, vertices[2].pos_camera); // assumes 'flat' quads
                Vec3f face_center_camera (0.0f);
                Vec3f barycenter_weight = Vec3f(1.0f / vertices.size()); // 1/3 for triangles, and 1/4 for quads
                for (int v = 0; v < vertices.size(); v++) 
                {
                    face_center_camera += Vec3f::hadamard_product(vertices[v].pos_camera, barycenter_weight);
                }

                Vec3f total_shading (0.0f);
                for (int l = 0; l < scene.lights.size(); l++)
                {
                    Light* a_light = scene.lights[l];
                    float shading  = calculate_shading(face_center_camera, face_normal_camera, *obj->mat, a_light, Vec3f(0.0f)); // shading done in camera space
                    total_shading += Vec3f::hadamard_product(a_light->color, Vec3f(shading));
                }
                total_shading = clampedVec3f(total_shading, 0.0f, 1.0f);

                for (int v = 0; v < vertices.size(); v++)
                {
                    vertices[v].shading = total_shading;
                }
            }
            else if (obj->shading == SHADING_TYPE::GOURAUD && !back_facing)
            {
                for (int v = 0; v < vertices.size(); v++)
                {
                    Vec3f total_shading (0.0f);
                    for (int l = 0; l < scene.lights.size(); l++)
                    {
                        Light* a_light = scene.lights[l];
                        float shading  = calculate_shading(vertices[v].pos_camera, vertices[v].norm_camera, *obj->mat, a_light, Vec3f(0.0f)); // shading done in camera space
                        total_shading  = Vec3f::hadamard_product(a_light->color, Vec3f(shading));
                    }
                    vertices[v].shading = clampedVec3f(total_shading, 0.0f, 1.0f);
                }
            }


            // Special case fill
            if (obj->fill_mode == FILL_MODE::COLORED_FACE_NORMALS && !back_facing)
            {
                Vec3f triangle_normal_world;
                if (obj->model->clockwise_winding)
                {
                    triangle_normal_world = triangle_normal(vertices[0].pos_world, vertices[vertices.size() - 1].pos_world, vertices[1].pos_world);
                }
                else
                {
                    triangle_normal_world = triangle_normal(vertices[0].pos_world, vertices[1].pos_world, vertices[vertices.size() - 1].pos_world);
                }
                Vec3f colored_triangle_normal = get_colored_normal(triangle_normal_world);

                for (int v = 0; v < vertices.size(); v++)
                {
                    vertices[v].color = colored_triangle_normal;
                }

                obj->fill_mode = FILL_MODE::VERTEX_COLORS;
            }


            // Render face to device buffers 
            if (obj->fill && vertices.size() == 3 && !back_facing)
            {
                if (obj->shading == SHADING_TYPE::PHONG)
                {
                    draw_triangle(vertices[0], vertices[1], vertices[2], scene.lights, obj->mat, obj->fill_mode, obj->texture);
                }
                else // none, flat, or gouraud shading
                {
                    draw_triangle(vertices[0], vertices[1], vertices[2], obj->fill_mode, obj->texture);
                }
            }
            else if (obj->fill && vertices.size() == 4 && !back_facing)
            {
                if (obj->shading == SHADING_TYPE::PHONG)
                {
                    // TODO: phong shading for quads
                }
                else // none, flat, or gouraud shading
                {
                    draw_quad(vertices[0], vertices[1], vertices[2], vertices[3], obj->fill_mode, obj->texture);
                }
            }

            if (obj->wireframe)
            {
                Vec3f wireframe_color = RED;
                float wireframe_thickness = 4.0f;
                Vec3f wireframe_epsilon (0.0f, 0.0f, 0.01f); // for preventing z-fighting with face 
                for (int i = 0; i < vertices.size(); i++)
                {
                    draw_line(vertices[i].pos_device + wireframe_epsilon, vertices[(i + 1) % vertices.size()].pos_device + wireframe_epsilon, wireframe_thickness, wireframe_color);
                }
            }
        }

        // Local axis (for debugging)
        // Vec3f world_origin = (device * Vec4f((projection * (camera * world * Vec4f(0.0f, 0.0f, 0.0f, 1.0f))).homogenize(), 1.0f)).xyz();
        // Vec3f world_x_axis_endpoint = (device * Vec4f((projection * (camera * world * Vec4f(1.0f, 0.0f, 0.0f, 1.0f))).homogenize(), 1.0f)).xyz();
        // Vec3f world_y_axis_endpoint = (device * Vec4f((projection * (camera * world * Vec4f(0.0f, 1.0f, 0.0f, 1.0f))).homogenize(), 1.0f)).xyz();
        // Vec3f world_z_axis_endpoint = (device * Vec4f((projection * (camera * world * Vec4f(0.0f, 0.0f, 1.0f, 1.0f))).homogenize(), 1.0f)).xyz();
        // draw_line(world_origin, world_x_axis_endpoint, 10.0f, Vec3f(0.5f, 0.0f, 0.0f), false);
        // draw_line(world_origin, world_y_axis_endpoint, 10.0f, Vec3f(0.0f, 0.5f, 0.0f), false);
        // draw_line(world_origin, world_z_axis_endpoint, 10.0f, Vec3f(0.0f, 0.0f, 0.5f), false);
    }

    // World axis (for debugging)
    // Vec3f world_origin = (device * Vec4f((projection * (camera * Vec4f(0.0f, 0.0f, 0.0f, 1.0f))).homogenize(), 1.0f)).xyz();
    // Vec3f world_x_axis_endpoint = (device * Vec4f((projection * (camera * Vec4f(1.0f, 0.0f, 0.0f, 1.0f))).homogenize(), 1.0f)).xyz();
    // Vec3f world_y_axis_endpoint = (device * Vec4f((projection * (camera * Vec4f(0.0f, 1.0f, 0.0f, 1.0f))).homogenize(), 1.0f)).xyz();
    // Vec3f world_z_axis_endpoint = (device * Vec4f((projection * (camera * Vec4f(0.0f, 0.0f, 1.0f, 1.0f))).homogenize(), 1.0f)).xyz();
    // draw_line(world_origin, world_x_axis_endpoint, 10.0f, RED, false);
    // draw_line(world_origin, world_y_axis_endpoint, 10.0f, GREEN, false);
    // draw_line(world_origin, world_z_axis_endpoint, 10.0f, BLUE, false);

    if (scene.skybox != nullptr)
    {
        // For orthographic it won't be actual dimensions of virtual screen
        // But, this is better since skybox won't be tied to 'zoom' of orthographic camera
        float skybox_virtual_screen_height = 1.0f;
        float skybox_virtual_screen_width  = aspect_ratio;
        float skybox_virtual_screen_near   = (aspect_ratio/2.0f)/tan(scene.skybox->fov/2.0f);

        // From device to ndc to virtual screen to skybox space
        Mat4x4f device_inv = scale(Vec3f(1.0f/ndc_scale_factor.x, 1.0f/ndc_scale_factor.y, 1.0f)) * translation(-device_center);
        Mat4x4f ndc_inv = scale(Vec3f((skybox_virtual_screen_width / ndc_width) / 2.0f, (skybox_virtual_screen_height / ndc_height) / 2.0f, 1.0f));
        Mat4x4f camera_for_skybox = Mat4x4f(camera.truncated().transposed()) * translation(Vec3f(0.0f, 0.0f, -skybox_virtual_screen_near));
        Mat4x4f skybox_rotation_inv = rotation_z(-scene.skybox->roll) * rotation_x(-scene.skybox->pitch) * rotation_y(-scene.skybox->yaw);

        // Cube is 1x1x1
        // Face of cube is 1x1
        float half_length_of_cube = 1.0f / 2.0f;

        // Project points on virtual screen to points on the cube
        // Then create ray vector, use ray to figure out which cube face to sample and where to sample it
        
        for (int y = 0; y < supersample_device_height; y++)
        {
            for (int x = 0; x < supersample_device_width; x++)
            {
                if (z_buffer[y][x] != std::numeric_limits<float>::max()) continue;

                Vec3f pixel_center = Vec3f(0.5f + x, 0.5f + y, 0.0f);
                Vec3f ray = (skybox_rotation_inv * camera_for_skybox * ndc_inv * device_inv * Vec4f(pixel_center, 1.0f)).xyz().normalize();

                // Find the dominant axis (the largest absolute value in ray.x, ray.y, ray.z)
                float absX = std::abs(ray.x);
                float absY = std::abs(ray.y);
                float absZ = std::abs(ray.z);

                FACE face;
                Vec3f ray_relative_to_face;

                if (absZ >= absX && absZ >= absY && ray.z < 0.0f)
                {
                    ray_relative_to_face = Vec3f(ray.x, ray.y, std::abs(ray.z));
                    face = FRONT;
                }
                else if (absZ >= absX && absZ >= absY && ray.z > 0.0f)
                {
                    ray_relative_to_face = Vec3f(-ray.x, ray.y, std::abs(ray.z));
                    face = BACK;
                }
                else if (absY >= absX && absY >= absZ && ray.y > 0.0f)
                {
                    ray_relative_to_face = Vec3f(ray.x, ray.z, std::abs(ray.y));
                    face = TOP;
                }
                else if (absY >= absX && absY >= absZ && ray.y < 0.0f)
                {
                    ray_relative_to_face = Vec3f(ray.x, -ray.z, std::abs(ray.y));
                    face = BOTTOM;
                }
                else if (absX >= absY && absX >= absZ && ray.x > 0.0f)
                {
                    ray_relative_to_face = Vec3f(ray.z, ray.y, std::abs(ray.x));
                    face = RIGHT;
                }
                else
                {
                    ray_relative_to_face = Vec3f(-ray.z, ray.y, std::abs(ray.x));
                    face = LEFT;
                }

                // cube face is [0,1]x[0,1], so just translate origin
                Vec2f projected = Vec2f((ray_relative_to_face.x / ray_relative_to_face.z) * half_length_of_cube, (ray_relative_to_face.y / ray_relative_to_face.z) * half_length_of_cube);
                Vec2f uv = projected + Vec2f(0.5f, 0.5f);
                color_buffer[y][x] = sample_texture(uv, scene.skybox->textures[face]);
            }
        }
    }

    // Init. downsampled color buffer
    Vec3f** downsampled_color_buffer = new Vec3f*[device_height];
    for (int i = 0; i < device_height; i++)
    {
        downsampled_color_buffer[i] = new Vec3f[device_width];
    }

    // Downsample supersampled color buffer
    for (int y = 0; y < device_height; y++)
    {
        for (int x = 0; x < device_width; x++)
        {
            Vec3f sum (0.0f);
            int supersample_x = x * supersample_factor;
            int supersample_y = y * supersample_factor;
            for (int i = 0; i < supersample_factor * supersample_factor; i++)
            {
                sum = sum + color_buffer[supersample_y][supersample_x];

                supersample_x += 1;
                if (supersample_x % supersample_factor == 0)
                {
                    supersample_x = x * supersample_factor;
                    supersample_y++;
                }
            }

            Vec3f average = Vec3f::hadamard_product(sum, Vec3f(1.0f / (supersample_factor * supersample_factor)));
            downsampled_color_buffer[y][x] = average;
        }
    }

    // Save image
    TGAImage tga_image = TGAImage(device_width, device_height, TGAImage::RGB);
    for (int y = 0; y < device_height; y++)
    {
        for (int x = 0; x < device_width; x++)
        {
            tga_image.set(x, (device_height - 1 - y), to_tga_color(downsampled_color_buffer[y][x]));
        }
    }
    tga_image.write_tga_file(scene.metadata->save_location.c_str());

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken: " << (duration.count() / 1000.0f)  << " seconds" << std::endl;
    return 0;
}


Vec3f sample_texture(Vec2f uv, TGAImage* texture)
{
    assert(uv.x >= 0.0f && uv.x <= 1.0f);
    assert(uv.y >= 0.0f && uv.y <= 1.0f);

    Vec2i pixel = Vec2i((texture->get_width() - 1) * uv.x, (texture->get_height() - 1) * uv.y);
    return to_vec3f_color(texture->get(pixel.x, pixel.y));
}


Vec3f get_colored_normal(Vec3f normal)
{
    return Vec3f(std::abs(normal.x), std::abs(normal.y), std::abs(normal.z));
}


Vec3f to_vec3f_color(TGAColor color)
{
    return Vec3f::hadamard_product(Vec3f(color.r, color.g, color.b), Vec3f(1.0f / 255.0f));
}


TGAColor to_tga_color(Vec3f color)
{
    assert(color.x >= 0.0f && color.x <= 1.0f);
    assert(color.y >= 0.0f && color.y <= 1.0f);
    assert(color.z >= 0.0f && color.z <= 1.0f);
    
    return TGAColor(int(color.x * 255.99f), int(color.y * 255.99f), int(color.z * 255.99f), 255);
}


// Assumes light, point, and normal are all in SAME space
float calculate_shading(Vec3f p, Vec3f n, Material mat, Light* light, Vec3f camera)
{
    float ambient  = 0.0f;
    float diffuse  = 0.0f;
    float specular = 0.0f;
    
    Vec3f v = (camera - p).normalize(); // view vector (points towards camera)

    if (light->type == "directional")
    {
        Vec3f r = reflect(n, light->direction);
        bool back_face_lighting = light->direction * n >= 0.0f;

        // ambient  = (light.direction * n < 0.0f ? mat.k_a : 0.0f); // lambertian-ish ambient
        ambient  = mat.k_a;
        diffuse  = mat.k_d * max(-(light->direction * n), 0.0f);
        specular = mat.k_s * power(max(v * r, 0.0f), mat.shininess);
    }
    else // point
    {
        Vec3f light_dir = (p - light->pos).normalize(); // points towards direction light travels
        Vec3f r = reflect (n, light_dir);

        float radius = (p - light->pos).length(); // distance between point and light
        float inv_sqr = min(light->intensity / (radius * radius + 0.01f), 1.0f); // [0, 1]

        // ambient  = (light_dir * n < 0.0f ? mat.k_a * inv_sqr : 0.0f); // lambertian-ish ambient
        ambient  = mat.k_a * inv_sqr;
        diffuse  = mat.k_d * max(-(light_dir * n), 0.0f) * inv_sqr;
        specular = mat.k_s * power(max(v * r, 0.0f), mat.shininess) * inv_sqr;
    }

    float shading = clampedf(ambient + diffuse + specular, 0.0f, 1.0f);
    return shading;
}


// Device coordinates
// No shading, or textures for lines.
void draw_line(Vec3f v0, Vec3f v1, float thickness, Vec3f color, bool depth_check)
{
    std::vector<LinePixel> raster_line = rasterize_line(v0.xy(), v1.xy(), thickness, Vec2i(supersample_device_width, supersample_device_height));

    for (int i = 0; i < raster_line.size(); i++)
    {
        LinePixel pixel = raster_line[i];
        float interpolated_depth = v0.z * (1.0f - pixel.t) + v1.z * pixel.t;

        // Depth check
        if (depth_check && z_buffer[pixel.pixel.y][pixel.pixel.x] < interpolated_depth) continue; // TODO: cull fragments behind or very close to front of camera

        color_buffer[pixel.pixel.y][pixel.pixel.x] = color;
        z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
    }
}


// Gouraud shading
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, FILL_MODE fill, TGAImage* texture)
{
    std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y), Vec2i(supersample_device_width, supersample_device_height));

    // For each fragment
    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        // Depth check
        if (z_buffer[pixel.pixel.y][pixel.pixel.x] > interpolated_depth) // TODO: cull fragments behind or very close to front of camera
        {
            // Figure out fragment color
            Vec3f pixel_color;
            if (fill == FILL_MODE::TEXTURE)
            {
                Vec2f interpolated_uv = barycentric_vec2f(v0.uv, v1.uv, v2.uv, pixel.barycentric);
                Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
                pixel_color = to_vec3f_color(texture->get(texture_pixel.x, texture_pixel.y));
            }
            else if (fill == FILL_MODE::VERTEX_COLORS)
            {
                pixel_color = barycentric_vec3f(v0.color, v1.color, v2.color, pixel.barycentric);
            }
            else // COLORED_NORMALS
            {
                Vec3f interpolated_normal = barycentric_vec3f(v0.norm_world, v1.norm_world, v2.norm_world, pixel.barycentric).normalize();
                pixel_color = get_colored_normal(interpolated_normal);
            }

            // Shade the fragment's color
            Vec3f interpolated_light = barycentric_vec3f(v0.shading, v1.shading, v2.shading, pixel.barycentric);
            Vec3f shaded_pixel_color = Vec3f::hadamard_product(pixel_color, interpolated_light);

            color_buffer[pixel.pixel.y][pixel.pixel.x] = shaded_pixel_color;
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}


// Phong shading
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, std::vector<Light*> lights, Material* mat, FILL_MODE fill, TGAImage* texture)
{
    std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y), Vec2i(supersample_device_width, supersample_device_height));

    // For each fragment
    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        // Depth check
        if (z_buffer[pixel.pixel.y][pixel.pixel.x] > interpolated_depth) // TODO: cull fragments behind or very close to front of camera
        {
            Vec3f interpolated_point = barycentric_vec3f(v0.pos_camera, v1.pos_camera, v2.pos_camera, pixel.barycentric); // camera space
            Vec3f interpolated_normal = barycentric_vec3f(v0.norm_camera, v1.norm_camera, v2.norm_camera, pixel.barycentric).normalize(); // camera space

            // Calculate phong (per-pixel) shading

            Vec3f total_light (0.0f);
            for (int l = 0; l < lights.size(); l++)
            {   
                // Shading done in camera space
                float shading = calculate_shading(interpolated_point, interpolated_normal, *mat, lights[l], Vec3f(0.0f));
                total_light += Vec3f::hadamard_product(lights[l]->color, Vec3f(shading));
            }
            total_light = clampedVec3f(total_light, 0.0f, 1.0f);

            // Figure out pixel color

            Vec3f pixel_color;
            if (fill == FILL_MODE::TEXTURE)
            {
                Vec2f interpolated_uv = barycentric_vec2f(v0.uv, v1.uv, v2.uv, pixel.barycentric);
                Vec2i texture_pixel_coord = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
                pixel_color = to_vec3f_color(texture->get(texture_pixel_coord.x, texture_pixel_coord.y));
            }
            else if (fill == FILL_MODE::VERTEX_COLORS)
            {
                pixel_color = barycentric_vec3f(v0.color, v1.color, v2.color, pixel.barycentric);
            }
            else // COLORED_NORMALS
            {
                // World normals colored
                Vec3f interpolated_normal = barycentric_vec3f(v0.norm_world, v1.norm_world, v2.norm_world, pixel.barycentric).normalize();
                pixel_color = get_colored_normal(interpolated_normal);
            }

            // Shade pixel color
            Vec3f shaded_pixel_color = Vec3f::hadamard_product(pixel_color, total_light);

            color_buffer[pixel.pixel.y][pixel.pixel.x] = shaded_pixel_color;
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}


// Gouraud shading
void draw_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3, FILL_MODE fill, TGAImage* texture)
{
    std::vector<QuadPixel> raster_quad = rasterize_quad(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y), Vec2f(v3.pos_device.x, v3.pos_device.y), Vec2i(supersample_device_width, supersample_device_height));

    for (int i = 0; i < raster_quad.size(); i++)
    {
        QuadPixel pix = raster_quad[i];
        float interpolated_depth = v0.pos_device.z * pix.barycentric.x + v1.pos_device.z * pix.barycentric.y + v2.pos_device.z * pix.barycentric.z + v3.pos_device.z * pix.barycentric.w;

        // Depth check
        if (z_buffer[pix.pixel.y][pix.pixel.x] > interpolated_depth) // TODO: cull fragments behind or very close to front of camera
        {
            // Figure out pixel color

            Vec3f pixel_color;
            if (fill == FILL_MODE::TEXTURE)
            {
                Vec2f interpolated_uv = Vec2f::hadamard_product(v0.uv, Vec2f(pix.barycentric.x, pix.barycentric.x)) + Vec2f::hadamard_product(v1.uv, Vec2f(pix.barycentric.y, pix.barycentric.y)) + Vec2f::hadamard_product(v2.uv, Vec2f(pix.barycentric.z, pix.barycentric.z)) + Vec2f::hadamard_product(v3.uv, Vec2f(pix.barycentric.w, pix.barycentric.w));
                Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
                pixel_color = to_vec3f_color(texture->get(texture_pixel.x, texture_pixel.y));
            }
            else if (fill == FILL_MODE::VERTEX_COLORS)
            {
                pixel_color = Vec3f::hadamard_product(v0.color, Vec3f(pix.barycentric.x)) + Vec3f::hadamard_product(v1.color, Vec3f(pix.barycentric.y)) + Vec3f::hadamard_product(v2.color, Vec3f(pix.barycentric.z)) + Vec3f::hadamard_product(v3.color, Vec3f(pix.barycentric.w));
            }
            else // COLORED_NORMALS
            {
                Vec3f interpolated_normal = (Vec3f::hadamard_product(v0.norm_world, Vec3f(pix.barycentric.x)) + Vec3f::hadamard_product(v1.norm_world, Vec3f(pix.barycentric.y)) + Vec3f::hadamard_product(v2.norm_world, Vec3f(pix.barycentric.z)) + Vec3f::hadamard_product(v3.norm_world, Vec3f(pix.barycentric.w))).normalize();
                pixel_color = get_colored_normal(interpolated_normal);
            }

            // Shade pixel color
            Vec3f interpolated_light = Vec3f::hadamard_product(v0.shading, Vec3f(pix.barycentric.x)) + Vec3f::hadamard_product(v1.shading, Vec3f(pix.barycentric.y)) + Vec3f::hadamard_product(v2.shading, Vec3f(pix.barycentric.z)) + Vec3f::hadamard_product(v3.shading, Vec3f(pix.barycentric.w));
            Vec3f shaded_pixel_color = Vec3f::hadamard_product(pixel_color, interpolated_light);

            color_buffer[pix.pixel.y][pix.pixel.x] = shaded_pixel_color;
            z_buffer[pix.pixel.y][pix.pixel.x] = interpolated_depth;

            // float grid_thickness = 0.0010;
            // int number_of_grid_lines = 20;
            // for (int j = 1; j < number_of_grid_lines; j++)
            // {
            //     float grid_line_pos = j * (1.0f/number_of_grid_lines);
            //     if (std::abs(grid_line_pos - pix.u) < grid_thickness)
            //     {
            //         color_buffer[pix.pixel.y][pix.pixel.x] = WHITE;
            //     }
            //     else if (std::abs(grid_line_pos - pix.v) < grid_thickness)
            //     {
            //         color_buffer[pix.pixel.y][pix.pixel.x] = WHITE;
            //     }
            // }
        }
    }
}