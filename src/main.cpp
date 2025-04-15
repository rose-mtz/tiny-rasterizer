#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <limits>
#include "tgaimage.h"
#include "model.h"
#include "Scene.h"

const Vec3f WHITE (1.0f, 1.0f, 1.0f);
const Vec3f RED   (1.0f, 0.0f, 0.0f);
const Vec3f GREEN (0.0f, 1.0f, 0.0f);
const Vec3f BLUE  (0.0f, 0.0f, 1.0f);
const Vec3f BLACK (0.0f, 0.0f, 0.0f);

const TGAColor TGA_WHITE = TGAColor(255, 255, 255, 255);
const TGAColor TGA_RED   = TGAColor(255, 0,   0,   255);
const TGAColor TGA_GREEN = TGAColor(0,   255, 0,   255);

const int   device_width  = 1000;
const int   device_height = 1000;
const float virtual_screen_width = 0.50;
const float virtual_screen_height = 0.50;

float** z_buffer = nullptr;
TGAImage* image = nullptr;

struct TrianglePixel
{
    Vec2i pixel;
    Vec3f barycentric;
};

struct LinePixel
{
    Vec2i pixel;
    float t;
};

struct Vertex
{
    Vec3f pos_device;
    Vec3f pos_camera;
    Vec3f norm_camera;
    Vec2f uv;
    Vec3f shading;
    Vec3f color;
};


Vec3f normal_colored(Vec3f normal);
TGAColor to_tgacolor(Vec3f color);
Vec3f to_vec3fcolor(TGAColor color);
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness = 1.0f);
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c);
float calculate_shading(Vec3f p, Vec3f n, Material mat, Light light, Vec3f camera);
void draw_line(Vec3f v0, Vec3f v1, float thickness, Vec3f color);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, TGAImage* texture);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, TGAImage* texture, std::vector<Light*> lights, Mat4x4f camera, Mat3x3f camera_inv_trans, Material* mat);


int main()
{
    Scene scene ("scenes/scene.txt");
    image = new TGAImage(device_width, device_height, TGAImage::RGB);
    image->fill(to_tgacolor(scene.background_color));
    
    z_buffer = new float*[device_height];
    for (int i = 0; i < device_height; i++)
    {
        z_buffer[i] = new float[device_width];
        for (int j = 0; j < device_width; j++)
        {
            z_buffer[i][j] = -std::numeric_limits<float>::max();
        }
    }


    // Camera
    Vec3f up = Vec3f(0.0f, 1.0f, 0.0f); // should maybe refactor this out, move it ot the scene file
    Mat4x4f camera = look_at(scene.camera->pos, scene.camera->look_at, up);

    // Projection
    Mat4x4f projection;
    projection.cols[0].x = scene.camera->zoom;
    projection.cols[1].y = scene.camera->zoom;
    if (scene.camera->type == "perspective") // homogenize later
    {
        projection.cols[2] = Vec4f(0.0f, 0.0f, 1.0f, 1.0f);
        projection.cols[3] = Vec4f(0.0f, 0.0f, 0.0f, 0.0f);
    }

    // Device
    float virtual_scale_x = (device_width / virtual_screen_width) / 2.0f;
    float virtual_scale_y = (device_height / virtual_screen_height) / 2.0f;
    Mat3x3f device (
        Vec3f(    virtual_scale_x,                 0.0f, 0.0f),
        Vec3f(               0.0f,      virtual_scale_y, 0.0f),
        Vec3f(device_width / 2.0f, device_height / 2.0f, 1.0f)
    ); // takes in homogenies virtual screen coords


    for (int i = 0; i < scene.objects.size(); i++)
    {
        Object3D* obj = scene.objects[i];

        // World
        Mat4x4f world = get_transformation(obj->pos, obj->scale);

        Mat4x4f camera_world = camera * world;
        Mat3x3f camera_world_inv_trans = camera_world.truncated().inv().transposed();
        Mat3x3f camera_inv_trans = camera.truncated().inv().transposed();

        for (int f = 0; f < obj->model->nfaces(); f++)
        {
            std::vector<int> face = obj->model->face(f);
            Vec3f v0 = obj->model->vert(face[0]);
            Vec3f v1 = obj->model->vert(face[3]);
            Vec3f v2 = obj->model->vert(face[6]);
            Vec3f n0 = obj->model->norm(face[2]);
            Vec3f n1 = obj->model->norm(face[5]);
            Vec3f n2 = obj->model->norm(face[8]);
            Vec2f uv0 = obj->model->uv(face[1]);
            Vec2f uv1 = obj->model->uv(face[4]);
            Vec2f uv2 = obj->model->uv(face[7]);

            Vec4f v0_camera = camera_world * Vec4f(v0, 1.0f);
            Vec4f v1_camera = camera_world * Vec4f(v1, 1.0f);
            Vec4f v2_camera = camera_world * Vec4f(v2, 1.0f);

            Vec3f n0_camera = (camera_world_inv_trans * n0).normalize();
            Vec3f n1_camera = (camera_world_inv_trans * n1).normalize();
            Vec3f n2_camera = (camera_world_inv_trans * n2).normalize();

            // Culling: back facing triangle check
            // NOTE: some triangles that should be culled will not be culled
            //       but (hopefully) no triangle that can be seen (i.e not supposed to be culled)
            //       will not be culled
            Vec3f triangle_normal = get_triangle_normal(v0_camera.xyz(), v1_camera.xyz(), v2_camera.xyz()).normalize();
            if (triangle_normal * Vec3f(0.0f, 0.0f, 1.0f) <= -0.25f) continue;

            Vec4f v0_projected = projection * v0_camera;
            Vec4f v1_projected = projection * v1_camera;
            Vec4f v2_projected = projection * v2_camera;
            
            // Homogenize-ish
            Vec3f v0_virtual = Vec3f(v0_projected.x / std::abs(v0_projected.w), v0_projected.y / std::abs(v0_projected.w), v0_projected.z);
            Vec3f v1_virtual = Vec3f(v1_projected.x / std::abs(v1_projected.w), v1_projected.y / std::abs(v1_projected.w), v1_projected.z);
            Vec3f v2_virtual = Vec3f(v2_projected.x / std::abs(v2_projected.w), v2_projected.y / std::abs(v2_projected.w), v2_projected.z);

            Vec3f v0_device = device * Vec3f(v0_virtual.x, v0_virtual.y, 1.0f); v0_device.z = v0_virtual.z; // keep depth
            Vec3f v1_device = device * Vec3f(v1_virtual.x, v1_virtual.y, 1.0f); v1_device.z = v1_virtual.z;
            Vec3f v2_device = device * Vec3f(v2_virtual.x, v2_virtual.y, 1.0f); v2_device.z = v2_virtual.z;

            // Vec3f rnd_color;
            // rnd_color.x = (std::rand() % 101) / 100.0f;
            // rnd_color.y = (std::rand() % 101) / 100.0f;
            // rnd_color.z = (std::rand() % 101) / 100.0f;

            // Flat & gouraud shading

            Vec3f light_flat (0.0f);
            Vec3f light_gouraud_a (0.0f);
            Vec3f light_gouraud_b (0.0f);
            Vec3f light_gouraud_c (0.0f);
            if (obj->shading == "flat")
            {
                Vec3f triangle_center = interpolate_barycentric_vec3f(v0_camera.xyz(), v1_camera.xyz(), v2_camera.xyz(), Vec3f(1.0f / 3.0f));
                for (int l = 0; l < scene.lights.size(); l++)
                {
                    // World --> camera
                    Light light_camera = *scene.lights[l];
                    light_camera.pos = (camera * Vec4f(light_camera.pos, 1.0f)).xyz();
                    light_camera.direction = (camera_inv_trans * light_camera.direction).normalize();
                    
                    // Shading done in camera space
                    float shading = calculate_shading(triangle_center, triangle_normal, *obj->mat, light_camera, Vec3f(0.0f));
                    Vec3f light_color = light_camera.color;
                    light_flat = light_flat + Vec3f(light_color.x * shading, light_color.y * shading, light_color.z * shading);
                }
                
                light_flat = clampedVec3f(light_flat, 0.0f, 1.0f);
            }
            else if (obj->shading == "gouraud")
            {
                for (int l = 0; l < scene.lights.size(); l++)
                {
                    // World --> camera
                    Light light_camera = *scene.lights[l];
                    light_camera.pos = (camera * Vec4f(light_camera.pos, 1.0f)).xyz();
                    light_camera.direction = (camera_inv_trans * light_camera.direction).normalize();

                    // Shading done in camera space
                    float shading_a = calculate_shading(v0_camera.xyz(), n0_camera, *obj->mat, light_camera, Vec3f(0.0f));
                    float shading_b = calculate_shading(v1_camera.xyz(), n1_camera, *obj->mat, light_camera, Vec3f(0.0f));
                    float shading_c = calculate_shading(v2_camera.xyz(), n2_camera, *obj->mat, light_camera, Vec3f(0.0f));
                    Vec3f light_color = light_camera.color;
                    light_gouraud_a = light_gouraud_a + Vec3f(light_color.x * shading_a, light_color.y * shading_a, light_color.z * shading_a);
                    light_gouraud_b = light_gouraud_b + Vec3f(light_color.x * shading_b, light_color.y * shading_b, light_color.z * shading_b);
                    light_gouraud_c = light_gouraud_c + Vec3f(light_color.x * shading_c, light_color.y * shading_c, light_color.z * shading_c);
                }

                light_gouraud_a = clampedVec3f(light_gouraud_a, 0.0f, 1.0f);
                light_gouraud_b = clampedVec3f(light_gouraud_b, 0.0f, 1.0f);
                light_gouraud_c = clampedVec3f(light_gouraud_c, 0.0f, 1.0f);
            }

            // NOTE: Barycentric coordinates will be distorted by the non-affine transformation perspective projection
            //       So, unless they are re-adjusted for this distortion, they will be wrong
            //       I haven't noticed any problems caused by this, so for now I'll stick with this not correcting them

            if (scene.fill_mode)
            {
                if (obj->shading == "flat")
                {
                    Vertex vertex0 = {.pos_device = v0_device, .uv = uv0, .shading = light_flat};
                    Vertex vertex1 = {.pos_device = v1_device, .uv = uv1, .shading = light_flat};
                    Vertex vertex2 = {.pos_device = v2_device, .uv = uv2, .shading = light_flat};

                    draw_triangle(vertex0, vertex1, vertex2, obj->texture);
                }
                else if (obj->shading == "gouraud")
                {
                    Vertex vertex0 = {.pos_device = v0_device, .uv = uv0, .shading = light_gouraud_a};
                    Vertex vertex1 = {.pos_device = v1_device, .uv = uv1, .shading = light_gouraud_b};
                    Vertex vertex2 = {.pos_device = v2_device, .uv = uv2, .shading = light_gouraud_c};

                    draw_triangle(vertex0, vertex1, vertex2, obj->texture);
                }
                else if (obj->shading == "phong")
                {
                    Vertex vertex0 = {.pos_device = v0_device, .pos_camera = v0_camera.xyz(), .norm_camera = n0_camera, .uv = uv0, .shading = light_gouraud_a};
                    Vertex vertex1 = {.pos_device = v1_device, .pos_camera = v1_camera.xyz(), .norm_camera = n1_camera, .uv = uv1, .shading = light_gouraud_b};
                    Vertex vertex2 = {.pos_device = v2_device, .pos_camera = v2_camera.xyz(), .norm_camera = n2_camera, .uv = uv2, .shading = light_gouraud_c};

                    draw_triangle(vertex0, vertex1, vertex2, obj->texture, scene.lights, camera, camera_inv_trans, obj->mat);
                }
                else // shading none
                {
                    Vertex vertex0 = {.pos_device = v0_device, .uv = uv0, .shading = Vec3f(1.0f, 1.0f, 1.0f)};
                    Vertex vertex1 = {.pos_device = v1_device, .uv = uv1, .shading = Vec3f(1.0f, 1.0f, 1.0f)};
                    Vertex vertex2 = {.pos_device = v2_device, .uv = uv2, .shading = Vec3f(1.0f, 1.0f, 1.0f)};

                    draw_triangle(vertex0, vertex1, vertex2, obj->texture);
                }
            }

            if (scene.wireframe_mode)
            {
                // Wire frame render
                float line_thickness = 1.0f;
                float wireframe_epsilon = 0.01f; // for preventing z-fighting with triangle
                draw_line(Vec3f(v0_device.x, v0_device.y, v0_device.z + wireframe_epsilon), Vec3f(v1_device.x, v1_device.y, v1_device.z + wireframe_epsilon), line_thickness, RED);
                draw_line(Vec3f(v1_device.x, v1_device.y, v1_device.z + wireframe_epsilon), Vec3f(v2_device.x, v2_device.y, v2_device.z + wireframe_epsilon), line_thickness, RED);
                draw_line(Vec3f(v2_device.x, v2_device.y, v2_device.z + wireframe_epsilon), Vec3f(v0_device.x, v0_device.y, v0_device.z + wireframe_epsilon), line_thickness, RED);
            }
        }
    }

    image->flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image->write_tga_file("output.tga");

    return 0;
}


Vec3f normal_colored(Vec3f normal)
{
    return Vec3f(std::abs(normal.x), std::abs(normal.y), std::abs(normal.z));
}

Vec3f to_vec3fcolor(TGAColor color)
{
    return Vec3f(color.r, color.g, color.b) * (1.0f / 255.0f);
}


TGAColor to_tgacolor(Vec3f color)
{
    return TGAColor(int(color.x * 255.99), int(color.y * 255.99), int(color.z * 255.99), 255);
}


// Rasterize line
// Points must be in device coordinates
// Thickness is in pixels
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness)
{   
    // KNOWN BUG: bounding box clips thickness.

    std::vector<LinePixel> pixels;

    // Bounding box
    Vec2i min, max;
    min.x = std::min(p0.x, p1.x);
    min.y = std::min(p0.y, p1.y);
    max.x = std::max(p0.x, p1.x);
    max.y = std::max(p0.y, p1.y);
    // Clip bounding box
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, device_width - 1);
    max.y = std::min(max.y, device_height - 1);

    Vec2f dir = Vec2f(p1.x - p0.x, p1.y - p0.y).normalize();
    float line_length = (p1 - p0).length();
    float half_thickness = thickness / 2.0f;

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f p (col + 0.5f, row + 0.5f);
            float projected_distance = (p - p0) * dir;
            float perp_distance = ((p - p0) - (dir * projected_distance)).length();
            float t = projected_distance / line_length;

            if (perp_distance <= half_thickness)
            {
                pixels.push_back(LinePixel {Vec2i(col, row), t});
            }
        }
    }

    return pixels;
}


// Rasterize a triangle
// Points must be in device coordinates
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c)
{
    std::vector<TrianglePixel> pixels;

    // Bounding box for triangle (pixel grid coordinate system)
    Vec2i min;
    Vec2i max;
    min.x = std::min(a.x, std::min(b.x, c.x));
    min.y = std::min(a.y, std::min(b.y, c.y));
    max.x = std::max(a.x, std::max(b.x, c.x));
    max.y = std::max(a.y, std::max(b.y, c.y));
    // Clip box to be within pixel grid
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, device_width - 1);
    max.y = std::min(max.y, device_height - 1);

    Mat3x3f D (Vec3f(a.x, a.y, 1.0f), Vec3f(b.x, b.y, 1.0f), Vec3f(c.x, c.y, 1.0f));
    float determinant = D.determinant();

    if (std::abs(determinant) < 1e-5f)
    {
        std::cout << "degenerate triangle, you little bastard!\n"; // for testing, haven't seen it log once yet
        return pixels;
    };

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f pixel_center (col + 0.5f, row + 0.5f);

            // Cramer's rule
            Vec3f B (pixel_center.x, pixel_center.y, 1.0f);
            float alpha = Mat3x3f(B,        D.col(1), D.col(2)).determinant() / determinant;
            float beta  = Mat3x3f(D.col(0), B,        D.col(2)).determinant() / determinant;
            float gamma = Mat3x3f(D.col(0), D.col(1), B       ).determinant() / determinant;

            // Inside triangle check
            if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f)
            {
                pixels.push_back(TrianglePixel {Vec2i(col, row), Vec3f(alpha, beta, gamma)});
            }
        }
    }

    return pixels;
}


// Assumes light, point, and normal are all in SAME space
float calculate_shading(Vec3f p, Vec3f n, Material mat, Light light, Vec3f camera)
{
    float ambient  = 0.0f;
    float diffuse  = 0.0f;
    float specular = 0.0f;
    
    Vec3f v = (camera - p).normalize(); // view vector (points towards camera)

    if (light.type == "directional")
    {
        Vec3f r = reflect(n, light.direction);
        bool back_face_lighting = light.direction * n >= 0.0f;

        // ambient  = (light.direction * n < 0.0f ? mat.k_a : 0.0f); // lambertian-ish ambient
        ambient  = mat.k_a;
        diffuse  = mat.k_d * max(-(light.direction * n), 0.0f);
        specular = mat.k_s * power(max(v * r, 0.0f), mat.shininess);
    }
    else // point
    {
        Vec3f light_dir = (p - light.pos).normalize(); // points towards direction light travels
        Vec3f r = reflect (n, light_dir);

        float radius = (p - light.pos).length(); // distance between point and light
        float inv_sqr = min(light.intensity / (radius * radius + 0.01f), 1.0f); // [0, 1]

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
void draw_line(Vec3f v0, Vec3f v1, float thickness, Vec3f color)
{
    std::vector<LinePixel> raster_line = rasterize_line(v0.xy(), v1.xy(), thickness);

    for (int i = 0; i < raster_line.size(); i++)
    {
        LinePixel pixel = raster_line[i];
        float interpolated_depth = v0.z * (1.0f - pixel.t) + v1.z * pixel.t;

        if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
        {
            image->set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(color));
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}


// Phong shading
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, TGAImage* texture)
{
    // TODO: if texture == nullptr then use vertex color

    std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y));

    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = interpolate_barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
        {
            Vec2f interpolated_uv = clampedVec2f(interpolate_barycentric_vec2f(v0.uv, v1.uv, v2.uv, pixel.barycentric), 0.0f, 1.0f);
            Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
            Vec3f texture_color = to_vec3fcolor(texture->get(texture_pixel.x, texture_pixel.y));

            Vec3f interpolated_light = interpolate_barycentric_vec3f(v0.shading, v1.shading, v2.shading, pixel.barycentric);
            interpolated_light = clampedVec3f(interpolated_light, 0.0f, 1.0f); // just to make sure its [0,1]
            Vec3f shaded_texture = Vec3f(texture_color.x * interpolated_light.x, texture_color.y * interpolated_light.y, texture_color.z * interpolated_light.z);

            image->set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(shaded_texture));
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}


void draw_triangle(Vertex v0, Vertex v1, Vertex v2, TGAImage* texture, std::vector<Light*> lights, Mat4x4f camera, Mat3x3f camera_inv_trans, Material* mat)
{
    std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y));

    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = interpolate_barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
        {
            Vec2f interpolated_uv = clampedVec2f(interpolate_barycentric_vec2f(v0.uv, v1.uv, v2.uv, pixel.barycentric), 0.0f, 1.0f);
            Vec3f interpolated_point = interpolate_barycentric_vec3f(v0.pos_camera, v1.pos_camera, v2.pos_camera, pixel.barycentric); // camera
            Vec3f interpolated_normal = interpolate_barycentric_vec3f(v0.norm_camera, v1.norm_camera, v2.norm_camera, pixel.barycentric).normalize(); // camera space

            Vec3f light (0.0f);
            for (int l = 0; l < lights.size(); l++)
            {
                // World --> camera
                Light light_camera = *lights[l];
                light_camera.pos = (camera * Vec4f(light_camera.pos, 1.0f)).xyz();
                light_camera.direction = (camera_inv_trans * light_camera.direction).normalize();
                
                // Shading done in camera space
                float shading = calculate_shading(interpolated_point, interpolated_normal, *mat, light_camera, Vec3f(0.0f));
                Vec3f light_color = light_camera.color;
                light = light + Vec3f(light_color.x * shading, light_color.y * shading, light_color.z * shading);
            }
            light = clampedVec3f(light, 0.0f, 1.0f);

            Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
            Vec3f texture_color = to_vec3fcolor(texture->get(texture_pixel.x, texture_pixel.y));
            Vec3f shaded_texture = Vec3f(texture_color.x * light.x, texture_color.y * light.y, texture_color.z * light.z);

            image->set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(shaded_texture));
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}