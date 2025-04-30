#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "Math.h"


/**
 * CREDIT: Dmitry V. Sokolov
 * REPO: https://github.com/ssloy/tinyrenderer/tree/a175be75a8a9a773bdfae7543a372e3bc859e02f
 */


class Model 
{
private:
	std::vector<Vec3f> verts_;
	std::vector<Vec3f> norms_;
	std::vector<Vec3f> textures_;
	std::vector<std::vector<int>> faces_;

public:
	Model(const char *filename);
	~Model();

	int nverts();
	int nfaces();
	Vec3f vert  (int i);
	Vec2f uv    (int i);
	Vec3f color (int i);
	Vec3f norm  (int i);

	// v/t/n, ...
	std::vector<int> face(int idx);
	bool clockwise_winding;
};


#endif //__MODEL_H__