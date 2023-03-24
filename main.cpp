#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <vector>
#include <stdlib.h>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green   = TGAColor(0, 255, 0,   255);
int width = 100;
int height = 100;

void line(
    int x1, int y1,
    int x2, int y2,
    TGAImage &image,
    const TGAColor &color
)
{
    if (std::abs(x2-x1) > std::abs(y2-y1))
    {
        if (x2 - x1 < 0) {
            std::swap(x1, x2);
            std::swap(y1, y2);
        }
        bool neg = y1 > y2;
        int dx = x2 - x1;
        int dy = neg ? y1 - y2 : y2 - y1;
        
        int P = 2 * dy - dx;

        int x = x1;
        int y = y1;

        while (x <= x2)
        {
            image.set(x, y, color);
            ++x;
            if (P < 0)
            {
                P += 2 * dy;
            }
            else
            {
                P += 2 * dy - 2 * dx;
                y += neg? -1 : 1;
            }
        }
    }

    else
    {
        if (y2 - y1 < 0) {
            std::swap(x1, x2);
            std::swap(y1, y2);
        }
        bool neg = x1 > x2;
        int dx = neg ? x1 - x2 : x2 - x1;
        int dy = y2 - y1;
        int P = 2 * dx - dy;

        int x = x1;
        int y = y1;

        while (y <= y2)
        {
            image.set(x, y, color);
            ++y;
            if (P < 0)
            {
                P += 2 * dx;
            }
            else
            {
                P += 2 * dx - 2 * dy;
                x += neg ? -1 : 1;
            }
        }
    }
}

inline void line(Vec2i v1, Vec2i v2, TGAImage &image, const TGAColor &color)
{
    line(v1.x, v1.y, v2.x, v2.y, image, color);
}

Vec3f barycentric(Vec2i* pts, Vec2i P) { 
    Vec3f u = Vec3f(
        pts[2].x-pts[0].x,
        pts[1].x-pts[0].x,
        pts[0].x-P.x
    )^Vec3f(
        pts[2].y-pts[0].y,
        pts[1].y-pts[0].y,
        pts[0].y-P.y
    );
    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */
    if (std::abs(u.z)<1) return Vec3f(-1,1,1);
    return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z); 
} 

void hline(int x1, int x2, int y, TGAImage &image, const TGAColor &color)
{
    if (x1 > x2) std::swap(x1, x2);
    for (int x = x1; x < x2; ++x) image.set(x, y, color); 
}

void triangle_lines(Vec2i v1, Vec2i v2, Vec2i v3, TGAImage &image, TGAColor color)
{
    line(v1, v2, image, color);
    line(v2, v3, image, color);
    line(v3, v1, image, color);
}

void filled_triangle_float(Vec2i v1, Vec2i v2, Vec2i v3, TGAImage &image, TGAColor color)
{
    if (v1.y > v2.y) std::swap(v1, v2);
    if (v1.y > v3.y) std::swap(v1, v3);
    if (v2.y > v3.y) std::swap(v2, v3);

    float m1 = (float)(v2.x - v1.x) / (v2.y - v1.y);
    float m2 = (float)(v3.x - v1.x) / (v3.y - v1.y);

    int x1, x2;
    
    int x10 = v1.x;
    int y10 = v1.y;

    int x20 = v1.x;
    int y20 = v1.y;

    for (int y = v1.y; y <= v3.y; ++y)
    {
        if (y == v2.y)
        {
            m1 = (float)(v3.x - v2.x) / (v3.y - v2.y);
            x10 = v2.x;
            y10 = v2.y;
        }

        x1 = (y - y10) * m1 + x10;
        x2 = (y - y20) * m2 + x20;

        hline(x1, x2, y, image, color);
    }
}

void filled_triangle_bar(Vec2i *pts, TGAImage &image, TGAColor color) { 
    Vec2i bboxmin(image.get_width()-1,  image.get_height()-1); 
    Vec2i bboxmax(0, 0); 
    Vec2i clamp(image.get_width()-1, image.get_height()-1); 
    for (int i=0; i<3; i++) { 
        bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
	bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));

	bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
	bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
    } 
    Vec2i P; 
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) { 
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) { 
            Vec3f bc_screen  = barycentric(pts, P); 
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue; 
            image.set(P.x, P.y, color); 
        } 
    } 
}

void draw_mesh(const char* filename, TGAImage &image, TGAColor color)
{
    Model* model = new Model(filename);

    for (int i=0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i); 
        for (int j=0; j<3; j++) { 
            Vec3f v0 = model->vert(face[j]); 
            Vec3f v1 = model->vert(face[(j+1)%3]); 
            int x0 = (v0.x+1.)*width/2.; 
            int y0 = (v0.y+1.)*height/2.; 
            int x1 = (v1.x+1.)*width/2.; 
            int y1 = (v1.y+1.)*height/2.; 
            line(x0, y0, x1, y1, image, color); 
        } 
    }
}

void draw_object_rc(const char* filename, TGAImage &image)
{
    Model* model = new Model(filename);

    for (int i=0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i); 
        Vec2i screen_coords[3];
        for (int j=0; j<3; j++) { 
            Vec3f world_coords = model->vert(face[j]); 
            screen_coords[j] = Vec2i((world_coords.x+1.)*width/2., (world_coords.y+1.)*height/2.); 
        } 
        filled_triangle_bar(screen_coords, image, TGAColor(rand()%255, rand()%255, rand()%255, 255)); 
    }
}

void draw_object_wl(const char* filename, TGAImage &image, Vec3f &light_dir)
{
    Model* model = new Model(filename);

    for (int i=0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i); 
        Vec2i screen_coords[3]; 
        Vec3f world_coords[3]; 
        for (int j=0; j<3; j++) { 
            Vec3f v = model->vert(face[j]); 
            screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.); 
            world_coords[j]  = v; 
        } 
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
        n.normalize(); 
        float intensity = n*light_dir; 
        if (intensity>0) { 
            filled_triangle_bar(screen_coords, image, TGAColor(intensity*255, intensity*255, intensity*255, 255)); 
        } 
    }
}


int main(int argc, char** argv) {
    if (argc == 3) {
        width = std::atoi(argv[1]);
        height = std::atoi(argv[2]);

        if (width == 0) width = 1000;
        if (height == 0) height = 1000;
    }

	TGAImage image(width, height, TGAImage::RGB);
	Vec2i v1(0, 0);
	Vec2i v2(width - 1, 0);
	Vec2i v3(width / 2, height - 1);

    const char* filename = "head_wireframe.obj";
    Vec3f light_dir(1, 0, -1);
    // draw_object_rc(filename, image);
    draw_object_wl(filename, image, light_dir);
    
    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
