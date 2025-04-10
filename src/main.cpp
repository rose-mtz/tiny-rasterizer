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


Vec3f normal_colored(Vec3f normal);
TGAColor to_tgacolor(Vec3f color);
Vec3f to_vec3fcolor(TGAColor color);
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness = 1.0f);
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c);


int main()
{
    TGAImage image(device_width, device_height, TGAImage::RGB);
    Scene scene ("scenes/scene.txt");
    
    float** z_buffer = new float*[device_height];
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

            // Back facing triangle check
            // NOTE: might need to use some epsilon, but that will have its own problems
            Vec3f triangle_normal = get_triangle_normal(v0_camera.xyz(), v1_camera.xyz(), v2_camera.xyz()).normalize();
            if (triangle_normal * Vec3f(0.0f, 0.0f, 1.0f) <= 0.0f)
            {
                continue;
            }

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

            // Flat shading

            Vec3f light_flat (0.0f);
            if (obj->shading == "flat")
            {
                for (int l = 0; l < scene.lights.size(); l++)
                {
                    float diffuse = 0.0f;
                    if (scene.lights[l]->type == "directional")
                    {
                        Vec3f dir_to_light_world = scene.lights[l]->direction * -1.0f; // world space
                        Vec3f dir_to_light_camera = (camera_inv_trans * dir_to_light_world).normalize(); // camera space
                        diffuse = clampedf(dir_to_light_camera * triangle_normal, 0.0f, 1.0f);
                    }
                    else // point
                    {
                        Vec3f light_pos_world = scene.lights[l]->pos;
                        Vec3f light_pos_camera = (camera * Vec4f(light_pos_world, 1.0f)).xyz();
                        Vec3f triangle_center_camera = ((v0_camera + v1_camera + v2_camera) * (1.0f / 3.0f)).xyz(); // barycenter of triangle in camera space
                        Vec3f dir_to_light = (light_pos_camera - triangle_center_camera).normalize(); // Optimize: redundant calculation
                        float r = (triangle_center_camera - light_pos_camera).length(); // Optimize: redundant calculation
                        diffuse = clampedf(dir_to_light * triangle_normal, 0.0f, 1.0f); // lambertian shading
                        diffuse = diffuse * (1.0f / (r * r + 0.01f)) * scene.lights[l]->intensity; // geometric attenuation (drop off)
                    }

                    Vec3f light_color = scene.lights[l]->color;
                    light_flat = light_flat + Vec3f(light_color.x * diffuse, light_color.y * diffuse, light_color.z * diffuse);
                }
                light_flat = clampedVec3f(light_flat, 0.0f, 1.0f);
            }

            // Gouraud shading

            Vec3f light_gouraud_a (0.0f);
            Vec3f light_gouraud_b (0.0f);
            Vec3f light_gouraud_c (0.0f);
            if (obj->shading == "gouraud")
            {
                for (int l = 0; l < scene.lights.size(); l++)
                {
                    float diffuse_a = 0.0f;
                    float diffuse_b = 0.0f;
                    float diffuse_c = 0.0f;
                    if (scene.lights[l]->type == "directional")
                    {
                        Vec3f dir_to_light_world = scene.lights[l]->direction * -1.0f; // world
                        Vec3f dir_to_light_camera = (camera_inv_trans * dir_to_light_world).normalize(); // camera

                        diffuse_a = clampedf(dir_to_light_camera * n0_camera, 0.0f, 1.0f);
                        diffuse_b = clampedf(dir_to_light_camera * n1_camera, 0.0f, 1.0f);
                        diffuse_c = clampedf(dir_to_light_camera * n2_camera, 0.0f, 1.0f);
                    }
                    else // point
                    {
                        Vec3f light_pos_world = scene.lights[l]->pos;
                        Vec3f light_pos_camera = (camera * Vec4f(light_pos_world, 1.0f)).xyz();
                        
                        // Optimize: redundant calculations
                        Vec3f dir_to_light_a = (light_pos_camera - v0_camera.xyz()).normalize();
                        Vec3f dir_to_light_b = (light_pos_camera - v1_camera.xyz()).normalize();
                        Vec3f dir_to_light_c = (light_pos_camera - v2_camera.xyz()).normalize();
                        float r_a = (light_pos_camera - v0_camera.xyz()).length();
                        float r_b = (light_pos_camera - v1_camera.xyz()).length();
                        float r_c = (light_pos_camera - v2_camera.xyz()).length();

                        diffuse_a = clampedf(dir_to_light_a * n0_camera, 0.0f, 1.0f) * (1.0f / (r_a * r_a + 0.01f)) * scene.lights[l]->intensity;
                        diffuse_b = clampedf(dir_to_light_b * n1_camera, 0.0f, 1.0f) * (1.0f / (r_b * r_b + 0.01f)) * scene.lights[l]->intensity;
                        diffuse_c = clampedf(dir_to_light_c * n2_camera, 0.0f, 1.0f) * (1.0f / (r_c * r_c + 0.01f)) * scene.lights[l]->intensity;
                    }

                    Vec3f light_color = scene.lights[l]->color;
                    light_gouraud_a = light_gouraud_a + Vec3f(light_color.x * diffuse_a, light_color.y * diffuse_a, light_color.z * diffuse_a); 
                    light_gouraud_b = light_gouraud_b + Vec3f(light_color.x * diffuse_b, light_color.y * diffuse_b, light_color.z * diffuse_b);
                    light_gouraud_c = light_gouraud_c + Vec3f(light_color.x * diffuse_c, light_color.y * diffuse_c, light_color.z * diffuse_c);
                }

                light_gouraud_a = clampedVec3f(light_gouraud_a, 0.0f, 1.0f);
                light_gouraud_b = clampedVec3f(light_gouraud_b, 0.0f, 1.0f);
                light_gouraud_c = clampedVec3f(light_gouraud_c, 0.0f, 1.0f);
            }


            std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0_device.x, v0_device.y), Vec2f(v1_device.x, v1_device.y), Vec2f(v2_device.x, v2_device.y));
            for (int i = 0; i < raster_triangle.size(); i++)
            {
                TrianglePixel pixel = raster_triangle[i];
                float interpolated_depth = (v0_device.z * pixel.barycentric.x) + (v1_device.z * pixel.barycentric.y) + (v2_device.z * pixel.barycentric.z);

                if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
                {
                    Vec2f interpolated_uv;
                    Vec3f interpolated_normal;
                    interpolated_uv = clampedVec2f((uv0 * pixel.barycentric.x) + (uv1 * pixel.barycentric.y) + (uv2 * pixel.barycentric.z), 0.0f, 1.0f);
                    interpolated_normal = ((n0_camera * pixel.barycentric.x) + (n1_camera * pixel.barycentric.y) + (n2_camera * pixel.barycentric.z)).normalize();

                    // Per-pixel shading

                    Vec3f light_per_pixel (0.0f);
                    if (obj->shading == "per-pixel")
                    {
                        for (int l = 0; l < scene.lights.size(); l++)
                        {
                            float diffuse = 0.0f;
                            if (scene.lights[l]->type == "directional")
                            {
                                Vec3f dir_to_light_world = scene.lights[l]->direction * -1.0f; // world
                                Vec3f dir_to_light_camera = (camera_inv_trans * dir_to_light_world).normalize(); // camera
                                diffuse = clampedf(dir_to_light_camera * interpolated_normal, 0.0f, 1.0f);
                            }
                            else // point
                            {
                                // NOTE: will be slightly wrong for perspective cameras due perspective projection not being a affine transformation
                                //       and not correct for its projective warping

                                Vec3f light_pos_world = scene.lights[l]->pos;
                                Vec3f light_pos_camera = (camera * Vec4f(light_pos_world, 1.0f)).xyz();
                                Vec3f interpolated_point_camera = ((v0_camera * pixel.barycentric.x) + (v1_camera * pixel.barycentric.y) + (v2_camera * pixel.barycentric.z)).xyz(); // this part is the wrong part
                                Vec3f dir_to_light = (light_pos_camera - interpolated_point_camera).normalize(); // Optimize: redundant calculation
                                float r = (light_pos_camera - interpolated_point_camera).length(); // Optimize: redundant calculation
                                diffuse = clampedf(dir_to_light * interpolated_normal, 0.0f, 1.0f) * (1.0f / (r * r + 0.01f)) * scene.lights[l]->intensity;
                            }

                            Vec3f light_color = scene.lights[l]->color;
                            light_per_pixel = light_per_pixel + Vec3f(light_color.x * diffuse, light_color.y * diffuse, light_color.z * diffuse);
                        }
                        light_per_pixel = clampedVec3f(light_per_pixel, 0.0f, 1.0f);
                    }

                    Vec2i texture_pixel = Vec2i((obj->texture->get_width() - 1) * interpolated_uv.x, (obj->texture->get_height() - 1) * interpolated_uv.y);
                    Vec3f texture_color = to_vec3fcolor(obj->texture->get(texture_pixel.x, texture_pixel.y));

                    Vec3f shaded_texture;
                    if (obj->shading == "flat")
                    {
                        shaded_texture = Vec3f(texture_color.x * light_flat.x, texture_color.y * light_flat.y, texture_color.z * light_flat.z);
                    }
                    else if (obj->shading == "gouraud")
                    {
                        Vec3f interpolated_light = (light_gouraud_a * pixel.barycentric.x) + (light_gouraud_b * pixel.barycentric.y) + (light_gouraud_c * pixel.barycentric.z);
                        interpolated_light = clampedVec3f(interpolated_light, 0.0f, 1.0f); // just to make sure its [0,1]
                        shaded_texture = Vec3f(texture_color.x * interpolated_light.x, texture_color.y * interpolated_light.y, texture_color.z * interpolated_light.z);
                    }
                    else // per-pixel
                    {
                        shaded_texture = Vec3f(texture_color.x * light_per_pixel.x, texture_color.y * light_per_pixel.y, texture_color.z * light_per_pixel.z);
                    }

                    image.set(pixel.pixel.x, pixel.pixel.y, to_tgacolor(shaded_texture));
                    z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
                }
            }
        }
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

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


// Rastersizes line
// Points must be in device coordinates
// Thickness is in pixels
// KNOWN BUG: bounding box clips thickness.
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness)
{
    // TODO: test this function!!!
    
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

            if (perp_distance <= half_thickness)
            {
                pixels.push_back(LinePixel {Vec2i(col, row), (1.0f - (projected_distance / line_length))});
            }
        }
    }

    return pixels;
}


// Rastersizes a triangle
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