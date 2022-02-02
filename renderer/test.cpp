#include "tgaimage.h"
#include <cmath>
#include <vector>
#include "model.h"
#include "geometry.h"
#include <iostream>



const TGAColor white = TGAColor(255,255,255,255);
const TGAColor red = TGAColor(255,0,0,255);
const TGAColor green = TGAColor(0,255,0,255);
const int width = 800;
const int height = 800;
Model *model = NULL;


void draw_line(TGAImage& image, int x0, int y0, int x1, int y1, TGAColor color){
	bool steep = false;
	if(std::abs(x0-x1) < std::abs(y0-y1)){
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if(x0 > x1){
		std::swap(x0, x1);
		std::swap(y0, y1);
	}
	
	int dx = x1-x0;
	int dy = y1-y0;
	float derror= std::abs(dy/float(dx));
	float error = 0;
	int y = y0;
	for(int x=x0; x<=x1; x++){
		if(steep){
			image.set(y, x, color);
		}
		else{
			image.set(x, y, color);
		}
		error += derror;
		if(error > 0.5){
			y+= (y1>y0?1:-1);
			error -= 1.;
		}
	}
}

/*
void triangle(TGAImage& image, Vec2i t0, Vec2i t1, Vec2i t2, TGAColor color){
	if(t0.y > t1.y) std::swap(t0, t1);
	if(t1.y > t2.y) std::swap(t1, t2);
	if(t0.y > t2.y) std::swap(t0, t2);
	int height = t2.y - t0.y;
	
	for(int y = t0.y; y <=t1.y; y++){
		int segment = (t1.y-t0.y+1);
		float alpha = (float)(y-t0.y)/(float)height; //interpolate the largest side from y0 to y1
		float beta = (float)(y-t0.y)/(float)segment; // other side from y0 to y1
		
		Vec2i A = t0 + (t2-t0)*alpha;
		Vec2i B = t0 + (t1-t0)*beta;
		if(A.x > B.x) std::swap(A, B); // x-coord on the largest side to x-coord at the other side
		for(int j = A.x; j<=B.x; j++){ //draw from a.x to b.x per y
			image.set(j, y, color); 
		}
		
	}
	
	for(int y = t1.y; y <=t2.y; y++){ //same shit for upper triangle
		int segment = (t2.y-t1.y+1);
		float alpha = (float)(y-t0.y)/(float)height; //interpolate the largest side from y1 to y2
		float beta = (float)(y-t1.y)/(float)segment; //The other sides from y1 to y2
		
		Vec2i A = t0 + (t2-t0)*alpha; //you still go from each vertex of A from 0
		Vec2i B = t1 + (t2-t1)*beta; // Already drew up until t1 so start from t1.
		if(A.x > B.x) std::swap(A, B);
		for(int j = A.x; j<=B.x; j++){
			image.set(j, y, color);
		}
		
	}
	
}
*/


	
Vec3f cross(Vec3f v1, Vec3f v2) {
    return Vec3f(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

// Multi Core method for drawing lines
// cross product (ABx ACx PA) and (ABy ACy PA) barry center is 1-u-v u v
Vec3f barrycentric(Vec2i *pts, Vec2i P){
	Vec3f u = cross(Vec3f(pts[2].x -pts[0].x, pts[1].x - pts[0].x, pts[0].x - P.x), Vec3f(pts[2].y -pts[0].y, pts[1].y- pts[0].y, pts[0].y-P.y));
	if (std::abs(u.z) < 1) return Vec3f(-1,1,1);
	return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
}


void triangle(TGAImage& image, Vec2i *pts, TGAColor color){
	Vec2i bboxmin(image.width()-1, image.height()-1);
	Vec2i bboxmax(0, 0);
	Vec2i clamp(image.width()-1, image.height()-1);
 	for(int i=0; i<3; i++){
		for(int j=0; j<2; j++){
			if(j==0){
				bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
				bboxmax.x  = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
			}else{
				bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));
				bboxmax.y  = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
			}
			
		}
	}
	
	//std::cout << bboxmin.x << " " << bboxmin.y << std::endl;
	//std::cout << bboxmax.x << " " << bboxmax.y << std::endl;
	
	Vec2i P;
	for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
		for(P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
			Vec3f lol = barrycentric(pts, P);
			if(lol.x < 0 || lol.y <0 || lol.z < 0) continue;
			image.set(P.x, P.y, color);
		}
		
	}
}



void draw_wire_frame(Model* model, TGAImage& image){
	for(int i=0; i< model->nfaces(); i++){
		std::vector<int> face = model -> face(i);
		for(int j=0; j<3; j++){
			Vec3f v0 = model->vert(face[j]);
			Vec3f v1 = model->vert(face[(j+1)%3]);
			int x0 = (v0.x+1.)*width/2.;
			int y0 = (v0.y+1.)*height/2.;
			int x1 = (v1.x+1.)*width/2.;
			int y1 = (v1.y+1.)*height/2.;
			draw_line(image, x0, y0, x1, y1, white);
		}
	}
}

void flate_shade_render(Model *model, TGAImage& image, Vec3f light_dir){
	for(int i=0; i<model->nfaces(); i++){
		std::vector<int> face = model -> face(i);
		Vec2i screen_coords[3];
		Vec3f world_coords[3];
		for(int j=0; j<3; j++){
			Vec3f v = model->vert(face[j]);
			screen_coords[j] = Vec2i((v.x+1)*image.width()/2., (v.y+1)*image.height()/2.);
			world_coords[j] = v; 
		}
		//triangle(image, screen_coords, TGAColor(rand()%255, rand()%255, rand()%255, 255));
		
		Vec3f n = (world_coords[2] - world_coords[0])^(world_coords[1]-world_coords[0]);
		n.normalize();
		float intensity = n * light_dir;
		//std::cout << in << std::endl;
		if(intensity > 0){
			triangle(image, screen_coords, TGAColor(intensity*255, intensity*255, intensity*255, 255));
		}
	}
}

int main(int argc, char** argv){
	
	if(2==argc){
		model = new Model(argv[1]);
	}
	else {
		model = new Model("head.obj");
	}
	
	TGAImage image(width, height, TGAImage::RGB);
	Vec3f light_dir(0,0,-1);
	flate_shade_render(model, image, Vec3f(0,0,-1));
	
	//image.flip_vertically();
	image.write_tga_file("out.tga");
	//delete model;
	return 0;
}