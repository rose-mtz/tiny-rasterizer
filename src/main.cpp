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

int supersample_factor;
int device_width;
int device_height;
int supersample_device_width;
int supersample_device_height;
float aspect_ratio;

float virtual_screen_width;
float virtual_screen_height;

float** z_buffer = nullptr;
Vec3f** color_buffer = nullptr;

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

struct QuadPixel
{
    Vec2i pixel;
    float alpha, beta, gamma, lambda;
};

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
}; // TODO: reorder struct members

enum FILL_MODE { VERTEX_UV, VERTEX_COLOR, COLORED_NORMALS };


Vec3f normal_colored(Vec3f normal);
TGAColor to_tgacolor(Vec3f color);
Vec3f to_vec3fcolor(TGAColor color);
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness = 1.0f);
std::vector<TrianglePixel> rasterize_triangle(Vec2f a, Vec2f b, Vec2f c);
float calculate_shading(Vec3f p, Vec3f n, Material mat, Light light, Vec3f camera);
void draw_line(Vec3f v0, Vec3f v1, float thickness, Vec3f color);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, FILL_MODE fill, TGAImage* texture = nullptr);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, std::vector<Light*> lights, Mat4x4f camera, Mat3x3f camera_inv_trans, Material* mat, FILL_MODE fill, TGAImage* texture = nullptr);
void draw_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3, FILL_MODE fill, TGAImage* texture);

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Error: Need to pass in scene as command-line argument.\n";
        return 0;
    }

    Scene scene (argv[1]);
    supersample_factor = scene.metadata->supersample_factor;
    device_width = scene.metadata->width_pixels;
    device_height = (scene.metadata->aspect_ratio.y * device_width) / scene.metadata->aspect_ratio.x;
    aspect_ratio = ((float) device_width) / ((float) device_height); // the one true aspect ratio
    supersample_device_width = device_width * supersample_factor;
    supersample_device_height = device_height * supersample_factor;
    virtual_screen_height = 1.0f;
    virtual_screen_width = aspect_ratio;
    
    z_buffer = new float*[supersample_device_height];
    for (int i = 0; i < supersample_device_height; i++)
    {
        z_buffer[i] = new float[supersample_device_width];
        for (int j = 0; j < supersample_device_width; j++)
        {
            z_buffer[i][j] = -std::numeric_limits<float>::max();
        }
    }

    color_buffer = new Vec3f*[supersample_device_height];
    for (int i = 0; i < supersample_device_height; i++)
    {
        color_buffer[i] = new Vec3f[supersample_device_width];
        for (int j = 0; j < supersample_device_width; j++)
        {
            color_buffer[i][j] = scene.metadata->background_color;
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
    float virtual_scale_x = (supersample_device_width / virtual_screen_width) / 2.0f;
    float virtual_scale_y = (supersample_device_height / virtual_screen_height) / 2.0f;
    Mat3x3f device (
        Vec3f(    virtual_scale_x,                 0.0f, 0.0f),
        Vec3f(               0.0f,      virtual_scale_y, 0.0f),
        Vec3f(supersample_device_width / 2.0f, supersample_device_height / 2.0f, 1.0f)
    ); // takes in homogenies virtual screen coords


    for (int i = 0; i < scene.objects.size(); i++)
    {
        Object3D* obj = scene.objects[i];

        // Transformation matrices
        Mat4x4f world = get_transformation(obj->pos, obj->scale);
        Mat3x3f world_inv_trans = world.truncated().inv().transposed();
        Mat4x4f camera_world = camera * world;
        Mat3x3f camera_world_inv_trans = camera_world.truncated().inv().transposed();
        Mat3x3f camera_inv_trans = camera.truncated().inv().transposed();

        for (int f = 0; f < obj->model->nfaces(); f++)
        {
            std::vector<int> face = obj->model->face(f);
            std::vector<Vertex> vertices (face.size() / 3); // 3 indices per vertex (pos, tex, norm)

            // Set up vertices
            for (int v = 0; v < vertices.size(); v++)
            {
                int pos_index  = (v * 3);
                int uv_index   = (v * 3) + 1;
                int norm_index = (v * 3) + 2;

                vertices[v].uv = obj->model->uv(face[uv_index]);
                vertices[v].color = obj->model->color(face[uv_index]);
                vertices[v].pos_camera = (camera_world * Vec4f(obj->model->vert(face[pos_index]), 1.0f)).xyz();
                vertices[v].norm_camera = (camera_world_inv_trans * obj->model->norm(face[norm_index])).normalize();
                vertices[v].pos_world = (world * Vec4f(obj->model->vert(face[pos_index]), 1.0f)).xyz();
                vertices[v].norm_world = (world_inv_trans * obj->model->norm(face[norm_index])).normalize();

                Vec4f vertex_projected = projection * Vec4f(vertices[v].pos_camera, 1.0f);
                Vec3f vertex_virtual = Vec3f(vertex_projected.x / std::abs(vertex_projected.w), vertex_projected.y / std::abs(vertex_projected.w), vertex_projected.z);
                vertices[v].pos_device = device * Vec3f(vertex_virtual.x, vertex_virtual.y, 1.0f); vertices[v].pos_device.z = vertex_virtual.z;
            }

            // Do shading
            if (obj->shading == "flat")
            {
                Vec3f total_shading (0.0f);

                // This won't work for quads
                Vec3f triangle_center = interpolate_barycentric_vec3f(vertices[0].pos_camera, vertices[1].pos_camera, vertices[2].pos_camera, Vec3f(1.0f / 3.0f)); // hardcoded weight (1/3 for triangles, 1/4 for quads)
                Vec3f triangle_normal_camera = get_triangle_normal(vertices[0].pos_camera, vertices[1].pos_camera, vertices[2].pos_camera).normalize();

                for (int l = 0; l < scene.lights.size(); l++)
                {
                    // World --> camera
                    Light light_camera = *scene.lights[l];
                    light_camera.pos = (camera * Vec4f(light_camera.pos, 1.0f)).xyz();
                    light_camera.direction = (camera_inv_trans * light_camera.direction).normalize();
                    
                    // Shading done in camera space
                    float shading = calculate_shading(triangle_center, triangle_normal_camera, *obj->mat, light_camera, Vec3f(0.0f));
                    Vec3f light_color = light_camera.color;
                    total_shading = total_shading + Vec3f(light_color.x * shading, light_color.y * shading, light_color.z * shading);
                }
                total_shading = clampedVec3f(total_shading, 0.0f, 1.0f);

                for (int v = 0; v < vertices.size(); v++)
                {
                    vertices[v].shading = total_shading;
                }
            }
            else if (obj->shading == "gouraud")
            {
                for (int v = 0; v < vertices.size(); v++)
                {
                    Vec3f total_shading (0.0f);
                    for (int l = 0; l < scene.lights.size(); l++)
                    {
                        // World --> camera
                        Light light_camera = *scene.lights[l];
                        light_camera.pos = (camera * Vec4f(light_camera.pos, 1.0f)).xyz();
                        light_camera.direction = (camera_inv_trans * light_camera.direction).normalize();

                        // Shading done in camera space
                        float shading = calculate_shading(vertices[v].pos_camera, vertices[v].norm_camera, *obj->mat, light_camera, Vec3f(0.0f));
                        Vec3f light_color = light_camera.color;
                        total_shading = total_shading + Vec3f(light_color.x * shading, light_color.y * shading, light_color.z * shading);
                    }
                    vertices[v].shading = clampedVec3f(total_shading, 0.0f, 1.0f);
                }
            }
            else if (obj->shading == "none")
            {
                for (int v = 0; v < vertices.size(); v++)
                {
                    vertices[v].shading = Vec3f(1.0f, 1.0f, 1.0f);
                }
            }

            // Triangle culling in device space
            // Vec3f triangle_normal_device = get_triangle_normal(Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, 0.0f), Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, 0.0f), Vec3f(vertices[2].pos_device.x, vertices[2].pos_device.y, 0.0f)).normalize();
            // if (triangle_normal_device == Vec3f(0.0f, 0.0f, -1.0f)) continue;

            // NOTE: Barycentric coordinates will be distorted by the non-affine transformation perspective projection
            //       So, unless they are re-adjusted for this distortion, they will be wrong
            //       I haven't noticed any problems caused by this, so for now I'll stick with this (not correcting them)

            if (vertices.size() == 3)
            {
                if (obj->fill_mode)
                {
                    if (obj->shading == "phong")
                    {
                        if (obj->colored_vertex_normals_mode) // colored vertex normals
                        {
                            draw_triangle(vertices[0], vertices[1], vertices[2], scene.lights, camera, camera_inv_trans, obj->mat, COLORED_NORMALS);
                        }
                        else if (obj->colored_triangle_normals_mode) // colored triangle normals
                        {
                            Vec3f triangle_normal_world = get_triangle_normal(vertices[0].pos_world, vertices[1].pos_world, vertices[2].pos_world);
                            Vec3f colored_triangle_normal = normal_colored(triangle_normal_world);
                            vertices[0].color = colored_triangle_normal;
                            vertices[1].color = colored_triangle_normal;
                            vertices[2].color = colored_triangle_normal;
                            draw_triangle(vertices[0], vertices[1], vertices[2], scene.lights, camera, camera_inv_trans, obj->mat, VERTEX_COLOR);
                        }
                        else // texture
                        {
                            draw_triangle(vertices[0], vertices[1], vertices[2], scene.lights, camera, camera_inv_trans, obj->mat, VERTEX_UV, obj->texture);
                        }
                    }
                    else // none, flat, or gouraud shading
                    {
                        if (obj->colored_vertex_normals_mode) // colored vertex normals
                        {
                            draw_triangle(vertices[0], vertices[1], vertices[2], COLORED_NORMALS);
                        }
                        else if (obj->colored_triangle_normals_mode) // colored triangle normals
                        {
                            Vec3f triangle_normal_world = get_triangle_normal(vertices[0].pos_world, vertices[1].pos_world, vertices[2].pos_world);
                            Vec3f colored_triangle_normal = normal_colored(triangle_normal_world);
                            vertices[0].color = colored_triangle_normal;
                            vertices[1].color = colored_triangle_normal;
                            vertices[2].color = colored_triangle_normal;
                            draw_triangle(vertices[0], vertices[1], vertices[2], VERTEX_COLOR);
                        }
                        else // texture
                        {
                            draw_triangle(vertices[0], vertices[1], vertices[2], VERTEX_UV, obj->texture);
                        }
                    }
                }

                if (obj->wireframe_mode)
                {
                    // Wire frame render
                    float line_thickness = 4.0f;
                    float wireframe_epsilon = 0.01f; // for preventing z-fighting with triangle
                    draw_line(Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, vertices[0].pos_device.z + wireframe_epsilon), Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, vertices[1].pos_device.z + wireframe_epsilon), line_thickness, RED);
                    draw_line(Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, vertices[1].pos_device.z + wireframe_epsilon), Vec3f(vertices[2].pos_device.x, vertices[2].pos_device.y, vertices[2].pos_device.z + wireframe_epsilon), line_thickness, RED);
                    draw_line(Vec3f(vertices[2].pos_device.x, vertices[2].pos_device.y, vertices[2].pos_device.z + wireframe_epsilon), Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, vertices[0].pos_device.z + wireframe_epsilon), line_thickness, RED);
                }
            }
            else // 4 vertex
            {
                draw_quad(vertices[0], vertices[1], vertices[2], vertices[3], VERTEX_COLOR, obj->texture);

                // // // Wire frame render
                // float line_thickness = 1.0f;
                // float wireframe_epsilon = 0.01f; // for preventing z-fighting with triangle
                // draw_line(Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, vertices[0].pos_device.z + wireframe_epsilon), Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, vertices[1].pos_device.z + wireframe_epsilon), line_thickness, RED);
                // draw_line(Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, vertices[1].pos_device.z + wireframe_epsilon), Vec3f(vertices[2].pos_device.x, vertices[2].pos_device.y, vertices[2].pos_device.z + wireframe_epsilon), line_thickness, RED);
                // draw_line(Vec3f(vertices[2].pos_device.x, vertices[2].pos_device.y, vertices[2].pos_device.z + wireframe_epsilon), Vec3f(vertices[3].pos_device.x, vertices[3].pos_device.y, vertices[3].pos_device.z + wireframe_epsilon), line_thickness, RED);
                // draw_line(Vec3f(vertices[3].pos_device.x, vertices[3].pos_device.y, vertices[3].pos_device.z + wireframe_epsilon), Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, vertices[0].pos_device.z + wireframe_epsilon), line_thickness, RED);
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
            tga_image.set(x, (device_height - 1 - y), to_tgacolor(downsampled_color_buffer[y][x]));
        }
    }
    tga_image.write_tga_file(scene.metadata->save_location.c_str());

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


std::vector<QuadPixel> rasterize_quad(Vec2f A, Vec2f B, Vec2f C, Vec2f D)
{
    std::vector<QuadPixel> pixels;

    for (int row = 0; row < supersample_device_height; row++)
    {
        for (int col = 0; col < supersample_device_width; col++)
        {
            Vec2f P (col + 0.5f, row + 0.5f);
            float u, v;

            float a = cross(A - B, A - B - D + C);
            float b = cross(P - A, A - B - D + C) + cross(A - B, D - A);
            float c = cross(P - A, D - A);

            if (a == 0) // degenerates to linear equation
            {
                u = -c / b;
                
                float numerator_x = (P.x - A.x + (u * A.x) - (u * B.x));
                float numerator_y = (P.y - A.y + (u * A.y) - (u * B.y));
                float denominator_x = (-A.x + (u * A.x) - (u * B.x) + D.x - (u * D.x) + (u * C.x));
                float denominator_y = (-A.y + (u * A.y) - (u * B.y) + D.y - (u * D.y) + (u * C.y));
                v = (denominator_x != 0.0f) ? numerator_x / denominator_x : numerator_y / denominator_y;

                if (u < 0.0f || u > 1.0f) continue;
                if (v < 0.0f || v > 1.0f) continue;
            }
            else // quadratic equation
            {
                float discriminant = (b * b) - (4 * a * c);
                if (discriminant < 0) continue;

                float u_plus = (-b + std::sqrt(discriminant)) / (2 * a);
                float u_minus = (-b - std::sqrt(discriminant)) / (2 * a);

                bool valid_u_plus = (u_plus >= 0.0f && u_plus <= 1.0f);
                bool valid_u_minus = (u_minus >= 0.0f && u_minus <= 1.0f);

                u = valid_u_plus ? u_plus : u_minus;

                float numerator_x = (P.x - A.x + (u * A.x) - (u * B.x));
                float numerator_y = (P.y - A.y + (u * A.y) - (u * B.y));
                float denominator_x = (-A.x + (u * A.x) - (u * B.x) + D.x - (u * D.x) + (u * C.x));
                float denominator_y = (-A.y + (u * A.y) - (u * B.y) + D.y - (u * D.y) + (u * C.y));

                v = (denominator_x != 0.0f) ? numerator_x / denominator_x : numerator_y / denominator_y;
                bool valid_v = (v >= 0.0f && v <= 1.0f);

                if (!valid_u_plus && !valid_u_minus) continue;
                if (!valid_v) continue;
            }

            QuadPixel pix;
            pix.pixel  = Vec2i(col, row);
            pix.alpha  = (1.0f - u) * (1.0f - v);
            pix.beta   = u * (1.0f - v);
            pix.gamma  = u * v;
            pix.lambda = v * (1.0f - u);
            pixels.push_back(pix);
        }
    }

    return pixels;
}


// Rasterize line
// Points must be in device coordinates
// Thickness is in pixels
std::vector<LinePixel> rasterize_line(Vec2f p0, Vec2f p1, float thickness)
{   
    std::vector<LinePixel> pixels;

    // Bounding box (includes thickness)
    Vec2i min, max;
    min.x = std::min(p0.x, p1.x) - thickness;
    min.y = std::min(p0.y, p1.y) - thickness;
    max.x = std::max(p0.x, p1.x) + thickness;
    max.y = std::max(p0.y, p1.y) + thickness;

    // Clip bounding box to device screen
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, supersample_device_width - 1);
    max.y = std::min(max.y, supersample_device_height - 1);

    Vec2f dir = Vec2f(p1.x - p0.x, p1.y - p0.y).normalize();
    float line_length = (p1 - p0).length();
    float half_thickness = thickness / 2.0f;

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f p (col + 0.5f, row + 0.5f);
            float projected_distance = (p - p0) * dir;

            if (projected_distance < 0.0f || projected_distance > line_length)
            {
                continue;
            }

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
    max.x = std::min(max.x, supersample_device_width - 1);
    max.y = std::min(max.y, supersample_device_height - 1);

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
            color_buffer[pixel.pixel.y][pixel.pixel.x] = color;
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}


// Gouraud shading
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, FILL_MODE fill, TGAImage* texture)
{
    std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y));

    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = interpolate_barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
        {
            Vec3f interpolated_light = interpolate_barycentric_vec3f(v0.shading, v1.shading, v2.shading, pixel.barycentric);

            Vec3f pixel_color;
            if (fill == VERTEX_UV)
            {
                Vec2f interpolated_uv = interpolate_barycentric_vec2f(v0.uv, v1.uv, v2.uv, pixel.barycentric);
                Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
                pixel_color = to_vec3fcolor(texture->get(texture_pixel.x, texture_pixel.y));
            }
            else if (fill == VERTEX_COLOR)
            {
                pixel_color = interpolate_barycentric_vec3f(v0.color, v1.color, v2.color, pixel.barycentric);
            }
            else // COLORED_NORMALS
            {
                Vec3f interpolated_normal = interpolate_barycentric_vec3f(v0.norm_world, v1.norm_world, v2.norm_world, pixel.barycentric).normalize();
                pixel_color = normal_colored(interpolated_normal);
            }

            Vec3f shaded_pixel_color = Vec3f::hadamard_product(pixel_color, interpolated_light);
            color_buffer[pixel.pixel.y][pixel.pixel.x] = shaded_pixel_color;
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}


// Phong shading
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, std::vector<Light*> lights, Mat4x4f camera, Mat3x3f camera_inv_trans, Material* mat, FILL_MODE fill, TGAImage* texture)
{
    std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y));

    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = interpolate_barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
        {
            Vec2f interpolated_uv = interpolate_barycentric_vec2f(v0.uv, v1.uv, v2.uv, pixel.barycentric);
            Vec3f interpolated_point = interpolate_barycentric_vec3f(v0.pos_camera, v1.pos_camera, v2.pos_camera, pixel.barycentric); // camera
            Vec3f interpolated_normal = interpolate_barycentric_vec3f(v0.norm_camera, v1.norm_camera, v2.norm_camera, pixel.barycentric).normalize(); // camera space

            Vec3f total_light (0.0f);
            for (int l = 0; l < lights.size(); l++)
            {
                // TODO: do something better, maybe let user handle these transformations

                // World --> camera
                Light light_camera = *lights[l];
                light_camera.pos = (camera * Vec4f(light_camera.pos, 1.0f)).xyz();
                light_camera.direction = (camera_inv_trans * light_camera.direction).normalize();
                
                // Shading done in camera space
                float shading = calculate_shading(interpolated_point, interpolated_normal, *mat, light_camera, Vec3f(0.0f));
                Vec3f light_color = light_camera.color;
                total_light = total_light + Vec3f::hadamard_product(light_color, Vec3f(shading));
            }
            total_light = clampedVec3f(total_light, 0.0f, 1.0f);

            Vec3f pixel_color;
            if (fill == VERTEX_UV)
            {
                Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
                pixel_color = to_vec3fcolor(texture->get(texture_pixel.x, texture_pixel.y));

            }
            else if (fill == VERTEX_COLOR)
            {
                pixel_color = interpolate_barycentric_vec3f(v0.color, v1.color, v2.color, pixel.barycentric);
            }
            else // COLORED_NORMALS
            {
                Vec3f interpolated_normal = interpolate_barycentric_vec3f(v0.norm_world, v1.norm_world, v2.norm_world, pixel.barycentric).normalize();
                pixel_color = normal_colored(interpolated_normal);
            }

            Vec3f shaded_pixel_color = Vec3f::hadamard_product(pixel_color, total_light);
            color_buffer[pixel.pixel.y][pixel.pixel.x] = shaded_pixel_color;
            z_buffer[pixel.pixel.y][pixel.pixel.x] = interpolated_depth;
        }
    }
}


// Gouraud shading
void draw_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3, FILL_MODE fill, TGAImage* texture)
{
    std::vector<QuadPixel> raster_quad = rasterize_quad(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y), Vec2f(v3.pos_device.x, v3.pos_device.y));

    for (int i = 0; i < raster_quad.size(); i++)
    {
        QuadPixel pix = raster_quad[i];
        float interpolated_depth = v0.pos_device.z * pix.alpha + v1.pos_device.z * pix.beta + v2.pos_device.z * pix.gamma + v3.pos_device.z * pix.lambda;

        if (z_buffer[pix.pixel.y][pix.pixel.x] <= interpolated_depth)
        {
            // Vec3f interpolated_light = interpolate_barycentric_vec3f(v0.shading, v1.shading, v2.shading, pixel.barycentric);
            // Vec3f interpolated_light = component_wise_product_vec3f(v0.shading, Vec3f(pix.alpha)) + component_wise_product_vec3f(v1.shading, Vec3f(pix.beta)) + component_wise_product_vec3f(v2.shading, Vec3f(pix.gamma)) + component_wise_product_vec3f(v3.shading, Vec3f(pix.lambda));

            Vec3f pixel_color = RED;
            // if (fill == TEXTURE)
            // {
            //     Vec2f interpolated_uv = component_wise_product_vec2f(v0.uv, Vec2f(pix.alpha, pix.alpha)) + component_wise_product_vec2f(v1.uv, Vec2f(pix.beta, pix.beta)) + component_wise_product_vec2f(v2.uv, Vec2f(pix.gamma, pix.gamma)) + component_wise_product_vec2f(v3.uv, Vec2f(pix.lambda, pix.lambda));


            //     Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
            //     pixel_color = to_vec3fcolor(texture->get(texture_pixel.x, texture_pixel.y));
            // }
            // else if (fill == VERTEX_COLOR)
            // {
            //     pixel_color = component_wise_product_vec3f(v0.color, Vec3f(pix.alpha)) + component_wise_product_vec3f(v1.color, Vec3f(pix.beta)) + component_wise_product_vec3f(v2.color, Vec3f(pix.gamma)) + component_wise_product_vec3f(v3.color, Vec3f(pix.lambda));
            // }
            // else // COLORED_NORMALS
            // {
            //     Vec3f interpolated_normal = (component_wise_product_vec3f(v0.norm_world, Vec3f(pix.alpha)) + component_wise_product_vec3f(v1.norm_world, Vec3f(pix.beta)) + component_wise_product_vec3f(v2.norm_world, Vec3f(pix.gamma)) + component_wise_product_vec3f(v3.norm_world, Vec3f(pix.lambda))).normalize();
            //     pixel_color = normal_colored(interpolated_normal);
            // }

            // Vec3f shaded_pixel_color = component_wise_product_vec3f(pixel_color, interpolated_light);
            pixel_color = Vec3f::hadamard_product(v0.color, Vec3f(pix.alpha)) + Vec3f::hadamard_product(v1.color, Vec3f(pix.beta)) + Vec3f::hadamard_product(v2.color, Vec3f(pix.gamma)) + Vec3f::hadamard_product(v3.color, Vec3f(pix.lambda));
            color_buffer[pix.pixel.y][pix.pixel.x] = pixel_color;
            z_buffer[pix.pixel.y][pix.pixel.x] = interpolated_depth;
        }
    }
}