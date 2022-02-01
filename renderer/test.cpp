#include "tgaimage.h"
#include <cmath>
#include <vector>
#include "model.h"
#include "geometry.h"



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

int main(int argc, char** argv){
	/*
	if(2==argc){
		model = new Model(argv[1]);
	}
	else {
		model = new Model("head.obj");
	}*/
	
	TGAImage image(width, height, TGAImage::RGB);
	
	Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)}; 
	Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
	Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 
	triangle(image, t0[0], t0[1], t0[2], red);
	triangle(image, t1[0], t1[1], t1[2], white);
	triangle(image, t2[0], t2[1], t2[2], green);
	
	//image.flip_vertically();
	image.write_tga_file("out.tga");
	//delete model;
	return 0;
}