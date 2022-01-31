#include "tgaimage.h"
#include <cmath>
#include <vector>
#include "model.h"
#include "geometry.h"



const TGAColor white = TGAColor(255,255,255,255);
const TGAColor red = TGAColor(255,0,0,255);
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

int main(int argc, char** argv){
	if(2==argc){
		model = new Model(argv[1]);
	}
	else {
		model = new Model("head.obj");
	}
	
	TGAImage image(width, height, TGAImage::RGB);
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
	//image.flip_vertically();
	image.write_tga_file("out.tga");
	delete model;
	return 0;
}