#include "Mesh.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


// for debugging
void print_mesh(Mesh& mesh)
{
    std::cout << "Mesh:\n";
    
    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        std::cout << "v "<< mesh.vertices[i].x << ' ' << mesh.vertices[i].y << ' ' << mesh.vertices[i].z << '\n'; 
    }

    for (int i = 0; i < mesh.vertex_textures.size(); i++)
    {
        std::cout << "vt "<< mesh.vertex_textures[i].x << ' ' << mesh.vertex_textures[i].y << ' ' << mesh.vertex_textures[i].z << '\n'; 
    }

    for (int i = 0; i < mesh.vertex_normals.size(); i++)
    {
        std::cout << "vn "<< mesh.vertex_normals[i].x << ' ' << mesh.vertex_normals[i].y << ' ' << mesh.vertex_normals[i].z << '\n'; 
    }

    for (int i = 0; i < mesh.faces.size(); i++)
    {

        std::cout << "f "; 
        for (int j = 0; j < 3; j++)
        {
            std::cout << mesh.faces[i][j].x << '/' << mesh.faces[i][j].y << '/' << mesh.faces[i][j].z << ' '; 
        }
        std::cout << '\n';
    }
}


Mesh load_mesh(const char* filename)
{
    Mesh mesh;
    
    std::ifstream in;
    in.open(filename);

    if (in.fail())
    {
        std::cout << "Failed to open: " << filename << '\n';
        return mesh;
    }

    std::string line;
    while(!in.eof())
    {
        std::getline(in, line);
        std::istringstream iss(line.c_str());

        if (line.size() == 0)
        {
            continue;
        }

        std::string line_type;
        iss >> line_type;
        if (line_type == "#")
        {
            continue;
        }
        else if (line_type == "v")
        {
            Vec3f vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            mesh.vertices.push_back(vertex);
        }
        else if (line_type == "vt")
        {
            Vec3f vertex_texture;
            iss >> vertex_texture.x >> vertex_texture.y >> vertex_texture.z;
            mesh.vertex_textures.push_back(vertex_texture);
        }
        else if (line_type == "vn")
        {
            Vec3f vertex_normal;
            iss >> vertex_normal.x >> vertex_normal.y >> vertex_normal.z;
            mesh.vertex_normals.push_back(vertex_normal);
        }
        else if (line_type == "f")
        {
            std::vector<Vec3i> face (3);

            for (int i = 0; i < 3; i++)
            {
                Vec3i f;
                char discard;

                iss >> f.x >> discard >> f.y >> discard >> f.z;
                f.x--; // zero index
                f.y--;
                f.z--;
                face[i] = f;
            }

            mesh.faces.push_back(face);
        }
    }

    return mesh;
}