#include "tgaimage.h"
#include <cmath>
#include <vector>
#include "model.h"
#include "geometry.h"
#include <iostream>



const TGAColor white = TGAColor(255,255,255,255);
const TGAColor red = TGAColor(255,0,0,255);
const TGAColor green = TGAColor(0,255,0,255);
const TGAColor blue = TGAColor(0,0,255,255);
const int width = 800;
const int height = 800;
const int depth  = 255;
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
Vec3f barrycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}


void triangle(TGAImage& image, Vec3f *pts, float* zbuffer, TGAColor color ){
	Vec2f bboxmin(image.get_width()-1, image.get_height()-1);
	Vec2f bboxmax(0, 0);
	Vec2f clamp(image.get_width()-1, image.get_height()-1);
 	for(int i=0; i<3; i++){
		for(int j=0; j<2; j++){
			if(j==0){
				bboxmin.x = std::max((float)0, std::min(bboxmin.x, pts[i].x));
				bboxmax.x  = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
			}else{
				bboxmin.y = std::max((float)0, std::min(bboxmin.y, pts[i].y));
				bboxmax.y  = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
			}
			
		}
	}
	
	//std::cout << bboxmin.x << " " << bboxmin.y << std::endl;
	//std::cout << bboxmax.x << " " << bboxmax.y << std::endl;
	
	Vec3f P;
	for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
		for(P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
			Vec3f lol = barrycentric(pts[0],pts[1], pts[2], P);
			if(lol.x < 0 || lol.y <0 || lol.z < 0) continue;
			P.z =0;
			for(int i=0; i<3; i++) P.z += lol[i]*pts[i][2];
			if(zbuffer[(int)(P.x + image.get_width()*P.y)] < P.z){
				zbuffer[(int)(P.x + image.get_width()*P.y)] = P.z;
				image.set(P.x, P.y, color);
			}
		}
		
	}
}



void tex_triangle(TGAImage& image, Vec3f *pts, Vec2i *uvs, float* zbuffer, Model* model, float intensity){
	Vec2f bboxmin(image.get_width()-1, image.get_height()-1);
	Vec2f bboxmax(0,0);
	Vec2f clamp(image.get_width()-1, image.get_height()-1);
	
	for(int i=0; i<3; i++){
		for(int j=0; j<2; j++){
			if(j==0){
				//std::cout << pts[i].x << " "; 
				bboxmin.x = std::max((float)0, std::min(bboxmin.x, pts[i].x));
				bboxmax.x  = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
			}else{
				//std::cout << pts[i].y << " ";
				bboxmin.y = std::max((float)0, std::min(bboxmin.y, pts[i].y));
				bboxmax.y  = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
			}
			
		}
		//std::cout << "\n";
	}
	
	std::cout << "x: " << uvs[0].x << " " << uvs[1].x << "\n";
	std::cout << "y: " << uvs[0].y << " " << uvs[1].y << "\n";
	
	Vec3f P;
	for(P.x = bboxmin.x; P.x<=bboxmax.x; P.x++){
		//std::cout << "Hello\n";
		for(P.y = bboxmin.y; P.y<=bboxmax.y; P.y++){
			Vec3f lol = barrycentric(pts[0], pts[1], pts[2], P);

			if(lol.x<0 || lol.y< 0 || lol.z <0) continue;
			P.z =0;
			for(int i=0; i<3; i++) P.z += pts[i][2]*lol[i];
			if(zbuffer[(int)(P.x + image.get_width()*P.y)] < P.z){
				zbuffer[(int)(P.x + image.get_width()*P.y)] = P.z;
				Vec2i uv = {0,0};
				for(int j=0; j<3; j++){
					uv.x += uvs[j].x*lol[j];
					uv.y += uvs[j].y*lol[j];
				}
				
				TGAColor color = model->diffuse(uv);
				image.set(P.x, P.y, TGAColor(color.r*intensity, color.g*intensity, color.b*intensity, 255));
			}
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

void rasterize1D(TGAImage& image, Vec2i p0, Vec2i p1, TGAColor color, int *ybuffer){
	if(p0.x > p1.x) std::swap(p0, p1);
	
	//int c = (p1.x*p0.y - p1.x*p0.y);
	//float m = (p1.y-p0.y)/(float)(p1.x - p0.x);
	for(int x = p0.x; x<=p1.x; x++){
		float t = (x-p0.x)/(float)(p1.x-p0.x);
		int y = p0.y*(1.-t) + p1.y*t;
		if(ybuffer[x] < y){
			ybuffer[x] = y;
			image.set(x, 0, color);
		}
	}
}

/*
void flate_shade_render(Model *model, TGAImage& image, float *zbuffer, Vec3f light_dir){
	for(int i=0; i<model->nfaces(); i++){
		std::vector<int> face = model -> face(i);
		Vec3f screen_coords[3];
		Vec3f world_coords[3];
		for(int j=0; j<3; j++){
			Vec3f v = model->vert(face[j]);
			screen_coords[j] = Vec2i((v.x+1)*image.get_width()/2., (v.y+1)*image.height()/2.);
			world_coords[j] = v; 
		}
		//triangle(image, screen_coords, TGAColor(rand()%255, rand()%255, rand()%255, 255));
		
		Vec3f n = cross((world_coords[2] - world_coords[0]),(world_coords[1]-world_coords[0]));
		n.normalize();
		float intensity = n * light_dir;
		//std::cout << in << std::endl;
		if(intensity > 0){
			triangle(image, screen_coords, zbuffer, TGAColor(intensity*255, intensity*255, intensity*255, 255));
		}
	}
}*/

Vec3f world2screen(Vec3f v){
	return Vec3f(int((v.x+1)*width/2. + .5),int((v.y+1.)*height/2. +.5), (v.z +1.)*depth/2);
}

int main(int argc, char** argv){
	
	freopen("out0.txt","w",stdout);
	
	
	if(2==argc){
		model = new Model(argv[1]);
	}
	else {
		model = new Model("head.obj");
	}
	
	TGAImage image(width, height, TGAImage::RGB);
	Vec3f light_dir(0,0,-1);
	//flate_shade_render(model, image, Vec3f(0,0,-1));
	
	float *zbuffer = new float[width*height];
	
	for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

	
	for(int i=0; i< model->nfaces(); i++){
		std::vector<int> face = model->face(i);
        Vec3f pts[3];
		Vec3f world_coords[3];
        for (int i=0; i<3; i++) pts[i] = world2screen(model->vert(face[i]));
		for(int j=0; j<3; j++){
			Vec3f v = model->vert(face[j]);
			world_coords[j] = v; 
		}
		Vec3f n = cross((world_coords[2] - world_coords[0]),(world_coords[1]-world_coords[0]));
		n.normalize();
		float intensity = n * light_dir;
		if(intensity > 0){
			Vec2i uv[3];
            for (int k=0; k<3; k++) {
                uv[k] = model->uv(i,k);
            }
			tex_triangle(image, pts, uv, zbuffer, model, intensity);
		}
        
	}
	
	
	
	image.flip_vertically();
	image.write_tga_file("out.tga");
	delete model;
	return 0;
}