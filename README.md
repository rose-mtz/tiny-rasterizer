# Tiny Rasterizer

A tiny rasterizer I am working on.

How it works. You specify the scene to render in a file. That scene will contain
models (meshes), lights, and meta data (supersampling, save location). The scene will
be rendered to TGA file.

**Current Features**
- Skybox rendering
- Triangle rasterization
- Quad rasterization
- Line rasterization
- Flat (per-privative) shading
- Gouraud (per-vertex) shading
- Phong (per-pixel) shading
- Specify scene to render from a file
- Supersampling
- Orthographic & perspective cameras
- Point & directional lights


## Some Rendered Images

Three heads each with different shading (flat, phong, gouraud). Taken with perspective camera. <br>
![rendered image](readme-images/three_african_heads_shaded_perspective.png)

Three heads each with different shading (flat, phong, gouraud). Taken with orthographic camera.<br>
![rendered image](readme-images/three_african_heads_shaded_orthographic.png)

Skybox. <br>
![rendered image](readme-images/skybox_front.png)

Cat made from a quad mesh. No shading. <br>
![rendered image](readme-images/cat_made_of_quads.png)

A quad made from two rasterized triangles. <br>
![rendered image](readme-images/square_vertex_colored_made_of_triangles.png)

A quad made from one rasterized quad. <br>
![rendered image](readme-images/vertex_colored_square_quad.png)

## Credit

I used Dmitry V. Sokolov's [tinyrenderer](https://github.com/ssloy/tinyrenderer) project as a guide for making this project. I did not follow 1-for-1 his project, instead used it as a guide. Although, I did make use of some of his code, like the TGA file import/save code, his vector class code, and the model loader code.