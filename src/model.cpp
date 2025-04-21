#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"


Model::Model(const char *filename) : verts_(), faces_()
{
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) 
    {
        std::cerr << "Error:: could not open model file " << filename << '\n';
        return;
    }
    std::string line;
    while (!in.eof())
    {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v "))
        {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            verts_.push_back(v);
        } 
        else if (!line.compare(0, 2, "vt") || !line.compare(0, 2, "vc"))
        {
            iss >> trash >> trash;
            Vec3f tex;
            for (int i=0;i<3;i++) iss >> tex.raw[i];
            textures_.push_back(tex);
        }
        else if (!line.compare(0, 2, "vn"))
        {
            iss >> trash >> trash;
            Vec3f n;
            for (int i=0;i<3;i++) iss >> n.raw[i];
            norms_.push_back(n);
        }
        else if (!line.compare(0, 2, "f "))
        {
            std::vector<int> f;
            int v_idx, t_idx, n_idx;
            iss >> trash;
            while (iss >> v_idx >> trash >> t_idx >> trash >> n_idx)
            {
                v_idx--; t_idx--; n_idx--; // zero index
                f.push_back(v_idx); f.push_back(t_idx); f.push_back(n_idx);
            }
            faces_.push_back(f);
        }
    }
}


Model::~Model()
{
}


int Model::nverts()
{
    return (int)verts_.size();
}


int Model::nfaces()
{
    return (int)faces_.size();
}


std::vector<int> Model::face(int idx)
{
    return faces_[idx];
}


Vec3f Model::vert(int i)
{
    return verts_[i];
}

Vec3f Model::norm(int i)
{
    return norms_[i];
}

Vec2f Model::uv(int i)
{
    return Vec2f(textures_[i].x, textures_[i].y);
}

Vec3f Model::color(int i)
{
    return textures_[i];
}