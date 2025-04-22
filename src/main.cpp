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
    float u, v;
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
};


Vec3f normal_colored(Vec3f normal);
TGAColor to_tgacolor(Vec3f color);
Vec3f to_vec3fcolor(TGAColor color);
float calculate_shading(Vec3f p, Vec3f n, Material mat, Light* light, Vec3f camera);
void draw_line(Vec3f v0, Vec3f v1, float thickness, Vec3f color);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, FILL_MODE fill, TGAImage* texture = nullptr);
void draw_triangle(Vertex v0, Vertex v1, Vertex v2, std::vector<Light*> lights, Material* mat, FILL_MODE fill, TGAImage* texture = nullptr);
void draw_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3, FILL_MODE fill, TGAImage* texture = nullptr);

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Error: Need to pass in scene as command-line argument.\n";
        return 0;
    }

    // Calculate dimensions of downsampled device, supersampled device, and virtual screen

    Scene scene (argv[1]);
    supersample_factor = scene.metadata->supersample_factor;
    device_width = scene.metadata->width_pixels;
    device_height = (scene.metadata->aspect_ratio.y * device_width) / scene.metadata->aspect_ratio.x;
    aspect_ratio = ((float) device_width) / ((float) device_height); // the one true aspect ratio
    supersample_device_width = device_width * supersample_factor;
    supersample_device_height = device_height * supersample_factor;
    virtual_screen_height = 1.0f;
    virtual_screen_width = aspect_ratio;

    // Initialize depth, and color buffer

    z_buffer = new float*[supersample_device_height];
    color_buffer = new Vec3f*[supersample_device_height];
    for (int i = 0; i < supersample_device_height; i++)
    {
        z_buffer[i] = new float[supersample_device_width];
        color_buffer[i] = new Vec3f[supersample_device_width];
        for (int j = 0; j < supersample_device_width; j++)
        {
            z_buffer[i][j] = -std::numeric_limits<float>::max();
            color_buffer[i][j] = scene.metadata->background_color;
        }
    }

    // Set up camera, projection, and device transformation matrices

    Vec3f up = Vec3f(0.0f, 1.0f, 0.0f); // TODO: use should control 'up' vector
    Mat4x4f camera = look_at(scene.camera->pos, scene.camera->look_at, up);
    Mat3x3f camera_inv_trans = camera.truncated().inv().transposed();

    Mat4x4f projection;
    projection.cols[0].x = scene.camera->zoom;
    projection.cols[1].y = scene.camera->zoom;
    if (scene.camera->type == "perspective")
    {
        projection.cols[2] = Vec4f(0.0f, 0.0f, 1.0f, 1.0f);
        projection.cols[3] = Vec4f(0.0f, 0.0f, 0.0f, 0.0f);
    }

    float virtual_scale_x = (supersample_device_width / virtual_screen_width) / 2.0f;
    float virtual_scale_y = (supersample_device_height / virtual_screen_height) / 2.0f;
    Mat3x3f device (
        Vec3f(    virtual_scale_x,                 0.0f, 0.0f),
        Vec3f(               0.0f,      virtual_scale_y, 0.0f),
        Vec3f(supersample_device_width / 2.0f, supersample_device_height / 2.0f, 1.0f)
    );

    // Transform lights

    for (int l = 0; l < scene.lights.size(); l++)
    {
        // Transfrom: world to camera
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
        Mat4x4f world = get_transformation(obj->pos, obj->scale);
        Mat3x3f world_inv_trans = world.truncated().inv().transposed();
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
                int pos_index  = (v * 3);
                int uv_or_color_index   = (v * 3) + 1;
                int norm_index = (v * 3) + 2;

                vertices[v].uv = obj->model->uv(face[uv_or_color_index]);
                vertices[v].color = obj->model->color(face[uv_or_color_index]);

                // Transform: local to world
                vertices[v].pos_world = (world * Vec4f(obj->model->vert(face[pos_index]), 1.0f)).xyz();
                vertices[v].norm_world = (world_inv_trans * obj->model->norm(face[norm_index])).normalize();
                // Transform: local to world to camera
                vertices[v].pos_camera = (camera_world * Vec4f(obj->model->vert(face[pos_index]), 1.0f)).xyz();
                vertices[v].norm_camera = (camera_world_inv_trans * obj->model->norm(face[norm_index])).normalize();
                // Transform: project to virtual screen
                Vec4f vertex_projected = projection * Vec4f(vertices[v].pos_camera, 1.0f);
                Vec3f vertex_virtual = Vec3f(vertex_projected.x / std::abs(vertex_projected.w), vertex_projected.y / std::abs(vertex_projected.w), vertex_projected.z);
                // Transform: virtual screen to device screen
                vertices[v].pos_device = device * Vec3f(vertex_virtual.x, vertex_virtual.y, 1.0f); vertices[v].pos_device.z = vertex_virtual.z;
            }

            // Do shading calculations

            if (obj->shading == "none")
            {
                for (int v = 0; v < vertices.size(); v++) vertices[v].shading = Vec3f(1.0f, 1.0f, 1.0f);
            }
            else if (obj->shading == "flat")
            {
                // Calculate point and normal at center of face
                Vec3f face_normal_camera = triangle_normal(vertices[0].pos_camera, vertices[1].pos_camera, vertices[2].pos_camera);
                Vec3f face_center_camera (0.0f);
                for (int v = 0; v < vertices.size(); v++) face_center_camera += Vec3f::hadamard_product(vertices[v].pos_camera, Vec3f(1.0f / vertices.size())); // weight is either 1/3 for triangles or 1/4 for quads

                // Calculate shading contribution of each light source
                Vec3f total_shading (0.0f);
                for (int l = 0; l < scene.lights.size(); l++)
                {
                    // Shading done in camera space
                    Light* a_light = scene.lights[l];
                    float shading = calculate_shading(face_center_camera, face_normal_camera, *obj->mat, a_light, Vec3f(0.0f));
                    total_shading += Vec3f::hadamard_product(a_light->color, Vec3f(shading));
                }
                total_shading = clampedVec3f(total_shading, 0.0f, 1.0f);

                for (int v = 0; v < vertices.size(); v++) vertices[v].shading = total_shading;
            }
            else if (obj->shading == "gouraud")
            {
                // Calculate total shading for each vertex
                for (int v = 0; v < vertices.size(); v++)
                {
                    // Calculate shading contribution of each light source
                    Vec3f total_shading (0.0f);
                    for (int l = 0; l < scene.lights.size(); l++)
                    {
                        // Shading done in camera space
                        Light* a_light = scene.lights[l];
                        float shading = calculate_shading(vertices[v].pos_camera, vertices[v].norm_camera, *obj->mat, a_light, Vec3f(0.0f));
                        total_shading = Vec3f::hadamard_product(a_light->color, Vec3f(shading));
                    }
                    vertices[v].shading = clampedVec3f(total_shading, 0.0f, 1.0f);
                }
            }

            // Back face culling calculation

            Vec3f triangle_normal_device = triangle_normal(Vec3f(vertices[0].pos_device.x, vertices[0].pos_device.y, 0.0f), Vec3f(vertices[1].pos_device.x, vertices[1].pos_device.y, 0.0f), Vec3f(vertices[2].pos_device.x, vertices[2].pos_device.y, 0.0f));
            bool back_facing = triangle_normal_device == Vec3f(0.0f, 0.0f, -1.0f); // TODO: use more numerically stable calculation

            // Render face to device buffer

            if (obj->fill && vertices.size() == 3 && !back_facing)
            {
                // Fill triangle

                if (obj->shading == "phong")
                {
                    if (obj->fill_mode == FILL_MODE::COLORED_FACE_NORMALS) // Colored normals fill
                    {
                        // Set colored face normal as vertex colors
                        Vec3f triangle_normal_world = triangle_normal(vertices[0].pos_world, vertices[1].pos_world, vertices[2].pos_world);
                        Vec3f colored_triangle_normal = normal_colored(triangle_normal_world);
                        vertices[0].color = colored_triangle_normal;
                        vertices[1].color = colored_triangle_normal;
                        vertices[2].color = colored_triangle_normal;
                        draw_triangle(vertices[0], vertices[1], vertices[2], scene.lights, obj->mat, FILL_MODE::VERTEX_COLORS);
                    }
                    else // Texture or vertex color fill
                    {
                        draw_triangle(vertices[0], vertices[1], vertices[2], scene.lights, obj->mat, obj->fill_mode, obj->texture);
                    }
                }
                else // none, flat, or gouraud shading
                {
                    if (obj->fill_mode == FILL_MODE::COLORED_FACE_NORMALS) // Colored normals fill
                    {
                        // Set colored face normal as vertex colors
                        Vec3f triangle_normal_world = triangle_normal(vertices[0].pos_world, vertices[1].pos_world, vertices[2].pos_world);
                        Vec3f colored_triangle_normal = normal_colored(triangle_normal_world);
                        vertices[0].color = colored_triangle_normal;
                        vertices[1].color = colored_triangle_normal;
                        vertices[2].color = colored_triangle_normal;
                        draw_triangle(vertices[0], vertices[1], vertices[2], FILL_MODE::COLORED_VERTEX_NORMALS);
                    }
                    else // Texture or vertex color fill
                    {
                        draw_triangle(vertices[0], vertices[1], vertices[2], obj->fill_mode, obj->texture);
                    }
                }
            }
            else if (obj->fill) // && !back_facing)
            {
                // Fill quad

                if (obj->shading == "phong")
                {
                    // TODO
                }
                else // none, flat, or gouraud shading
                {
                    if (obj->fill_mode == FILL_MODE::COLORED_FACE_NORMALS) // Colored normals fill
                    {
                        // Set colored face normal as vertex colors
                        Vec3f quad_normal_world = triangle_normal(vertices[0].pos_world, vertices[1].pos_world, vertices[2].pos_world);
                        Vec3f colored_quad_normal = normal_colored(quad_normal_world);
                        vertices[0].color = colored_quad_normal;
                        vertices[1].color = colored_quad_normal;
                        vertices[2].color = colored_quad_normal;
                        draw_quad(vertices[0], vertices[1], vertices[2], vertices[3], FILL_MODE::COLORED_VERTEX_NORMALS);
                    }
                    else // Texture or vertex color fill
                    {
                        draw_quad(vertices[0], vertices[1], vertices[2], vertices[3], obj->fill_mode, obj->texture);
                    }
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
    return Vec3f::hadamard_product(Vec3f(color.r, color.g, color.b), Vec3f(1.0f / 255.0f));
}


TGAColor to_tgacolor(Vec3f color)
{
    return TGAColor(int(color.x * 255.99f), int(color.y * 255.99f), int(color.z * 255.99f), 255);
}


std::vector<QuadPixel> rasterize_quad(Vec2f A, Vec2f B, Vec2f C, Vec2f D)
{
    std::vector<QuadPixel> pixels;

    // Bounding box (pixel grid coordinate system)
    Vec2i min;
    Vec2i max;
    min.x = std::min(A.x, std::min(B.x, std::min(C.x, D.x)));
    min.y = std::min(A.y, std::min(B.y, std::min(C.y, D.y)));
    max.x = std::max(A.x, std::max(B.x, std::max(C.x, D.x)));
    max.y = std::max(A.y, std::max(B.y, std::max(C.y, D.y)));
    // Clip box to be within pixel grid
    min.x = std::max(min.x, 0);
    min.y = std::max(min.y, 0);
    max.x = std::min(max.x, supersample_device_width - 1);
    max.y = std::min(max.y, supersample_device_height - 1);

    for (int row = min.y; row <= max.y; row++)
    {
        for (int col = min.x; col <= max.x; col++)
        {
            Vec2f P (col + 0.5f, row + 0.5f);
            float u = -1.0f;
            float v = -1.0f;

            float a = cross(A - B, A - B - D + C);
            float b = cross(P - A, A - B - D + C) + cross(A - B, D - A);
            float c = cross(P - A, D - A);

            if (a == 0) 
            {
                // Solve for u (degenerates to linear equation)
                u = -c / b;
            }
            else
            {
                // Solve for u (quadratic equation)
                float discriminant = (b * b) - (4 * a * c);
                if (discriminant < 0) continue;

                float u_plus = (-b + std::sqrt(discriminant)) / (2 * a);
                float u_minus = (-b - std::sqrt(discriminant)) / (2 * a);

                if (u_plus  >= 0.0f && u_plus  <= 1.0f) u = u_plus;
                if (u_minus >= 0.0f && u_minus <= 1.0f) u = u_minus;
            }

            // Solve for v
            float numerator_x = (P.x - A.x + (u * A.x) - (u * B.x));
            float numerator_y = (P.y - A.y + (u * A.y) - (u * B.y));
            float denominator_x = (-A.x + (u * A.x) - (u * B.x) + D.x - (u * D.x) + (u * C.x));
            float denominator_y = (-A.y + (u * A.y) - (u * B.y) + D.y - (u * D.y) + (u * C.y));
            v = (denominator_x != 0.0f) ? numerator_x / denominator_x : numerator_y / denominator_y;

            // Outside quad check
            if (u < 0.0f || u > 1.0f || v < 0.0f || v > 1.0f) continue;

            QuadPixel pix;
            pix.pixel  = Vec2i(col, row);
            pix.alpha  = (1.0f - u) * (1.0f - v);
            pix.beta   = u * (1.0f - v);
            pix.gamma  = u * v;
            pix.lambda = v * (1.0f - u);
            pix.u      = u;
            pix.v      = v;
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

    // Degenerate triangle check
    if (std::abs(determinant) < 1e-5f) return pixels;

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

    // For each fragment
    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        // Depth check
        if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
        {
            // Figure out fragment color
            Vec3f pixel_color;
            if (fill == FILL_MODE::TEXTURE)
            {
                Vec2f interpolated_uv = barycentric_vec2f(v0.uv, v1.uv, v2.uv, pixel.barycentric);
                Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
                pixel_color = to_vec3fcolor(texture->get(texture_pixel.x, texture_pixel.y));
            }
            else if (fill == FILL_MODE::VERTEX_COLORS)
            {
                pixel_color = barycentric_vec3f(v0.color, v1.color, v2.color, pixel.barycentric);
            }
            else // COLORED_NORMALS
            {
                Vec3f interpolated_normal = barycentric_vec3f(v0.norm_world, v1.norm_world, v2.norm_world, pixel.barycentric).normalize();
                pixel_color = normal_colored(interpolated_normal);
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
    std::vector<TrianglePixel> raster_triangle = rasterize_triangle(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y));

    // For each fragment
    for (int i = 0; i < raster_triangle.size(); i++)
    {
        TrianglePixel pixel = raster_triangle[i];
        float interpolated_depth = barycentric_f(v0.pos_device.z, v1.pos_device.z, v2.pos_device.z, pixel.barycentric);

        // Depth check
        if (z_buffer[pixel.pixel.y][pixel.pixel.x] <= interpolated_depth)
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
                pixel_color = to_vec3fcolor(texture->get(texture_pixel_coord.x, texture_pixel_coord.y));
            }
            else if (fill == FILL_MODE::VERTEX_COLORS)
            {
                pixel_color = barycentric_vec3f(v0.color, v1.color, v2.color, pixel.barycentric);
            }
            else // COLORED_NORMALS
            {
                // World normals colored
                Vec3f interpolated_normal = barycentric_vec3f(v0.norm_world, v1.norm_world, v2.norm_world, pixel.barycentric).normalize();
                pixel_color = normal_colored(interpolated_normal);
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
    std::vector<QuadPixel> raster_quad = rasterize_quad(Vec2f(v0.pos_device.x, v0.pos_device.y), Vec2f(v1.pos_device.x, v1.pos_device.y), Vec2f(v2.pos_device.x, v2.pos_device.y), Vec2f(v3.pos_device.x, v3.pos_device.y));

    for (int i = 0; i < raster_quad.size(); i++)
    {
        QuadPixel pix = raster_quad[i];
        float interpolated_depth = v0.pos_device.z * pix.alpha + v1.pos_device.z * pix.beta + v2.pos_device.z * pix.gamma + v3.pos_device.z * pix.lambda;

        // Depth check
        if (z_buffer[pix.pixel.y][pix.pixel.x] <= interpolated_depth)
        {
            // Calculate shading
            Vec3f interpolated_light = Vec3f::hadamard_product(v0.shading, Vec3f(pix.alpha)) + Vec3f::hadamard_product(v1.shading, Vec3f(pix.beta)) + Vec3f::hadamard_product(v2.shading, Vec3f(pix.gamma)) + Vec3f::hadamard_product(v3.shading, Vec3f(pix.lambda));

            // Figure out pixel color

            Vec3f pixel_color;
            if (fill == FILL_MODE::TEXTURE)
            {
                Vec2f interpolated_uv = Vec2f::hadamard_product(v0.uv, Vec2f(pix.alpha, pix.alpha)) + Vec2f::hadamard_product(v1.uv, Vec2f(pix.beta, pix.beta)) + Vec2f::hadamard_product(v2.uv, Vec2f(pix.gamma, pix.gamma)) + Vec2f::hadamard_product(v3.uv, Vec2f(pix.lambda, pix.lambda));
                Vec2i texture_pixel = Vec2i((texture->get_width() - 1) * interpolated_uv.x, (texture->get_height() - 1) * interpolated_uv.y);
                pixel_color = to_vec3fcolor(texture->get(texture_pixel.x, texture_pixel.y));
            }
            else if (fill == FILL_MODE::VERTEX_COLORS)
            {
                pixel_color = Vec3f::hadamard_product(v0.color, Vec3f(pix.alpha)) + Vec3f::hadamard_product(v1.color, Vec3f(pix.beta)) + Vec3f::hadamard_product(v2.color, Vec3f(pix.gamma)) + Vec3f::hadamard_product(v3.color, Vec3f(pix.lambda));
            }
            else // COLORED_NORMALS
            {
                Vec3f interpolated_normal = (Vec3f::hadamard_product(v0.norm_world, Vec3f(pix.alpha)) + Vec3f::hadamard_product(v1.norm_world, Vec3f(pix.beta)) + Vec3f::hadamard_product(v2.norm_world, Vec3f(pix.gamma)) + Vec3f::hadamard_product(v3.norm_world, Vec3f(pix.lambda))).normalize();
                pixel_color = normal_colored(interpolated_normal);
            }

            // Shade pixel color
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