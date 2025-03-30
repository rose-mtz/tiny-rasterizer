#include "Vec.h"
#include <vector>

struct Mesh
{
    // TODO: later remove vector use pointes instead, why? fuck vectors
    std::vector<Vec3f> vertices;
    std::vector<Vec3f> vertex_textures;
    std::vector<Vec3f> vertex_normals;
    std::vector<std::vector<Vec3i>> faces;
};

Mesh load_mesh(const char* filename);