//Author: Ben Loisch. Date completed: 07/08/2016.
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;
#define imagewidth 1920
#define imageheight 1080
#define imageSamples 9
#define ambientSamples 9
#define lightSamples 25
#define dotproductmax(a, b) ((a > b) ? a : b)
double specularReflectionCoefficient(450), phongExponent(100), plr(5), epsilon(0.00001), lightIntensity(1), PI(3.1415926), ambCoef(0.2);
class Vector {public: double x; double y; double z; Vector(double xi, double yi, double zi) : x(xi), y(yi), z(zi) {}
inline Vector Vector::normalize() { return Vector(x / sqrt((x * x) + (y * y) + (z * z)), y / sqrt((x * x) + (y * y) + (z * z)), z / sqrt((x * x) + (y * y) + (z * z))); } //return normalized vector
inline Vector Vector::operator-(const Vector &rhs) const { return Vector(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z); }
inline Vector Vector::operator+(const Vector &rhs) const { return Vector(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z); }
inline Vector Vector::operator*(const double rhs) const { return Vector(this->x * rhs, this->y * rhs, this->z * rhs); }
inline double Vector::operator*(const Vector &rhs) const { return ((this->x * rhs.x) + (this->y * rhs.y) + (this->z * rhs.z)); } //dot product
inline Vector Vector::operator^(const Vector &rhs) const { return Vector((this->y * rhs.z) - (this->z * rhs.y), (this->z * rhs.x) - (this->x * rhs.z), (this->x * rhs.y) - (this->y * rhs.x)); } }; //cross product
class RGB { public: double r; double g; double b; RGB(double ri, double gi, double bi) : r(ri), g(gi), b(bi) {} 
inline RGB& RGB::operator/=(const double rhs) { r = r / rhs; g = g / rhs; b = b / rhs; return *this; }
inline RGB&	RGB::operator+=(RGB &rhs) { r = r + rhs.r; g = g + rhs.g; b = b + rhs.b; return *this; }
inline RGB RGB::operator*(const double rhs) { return RGB(r * rhs, g * rhs, b * rhs); }
inline RGB RGB::operator+(RGB &rhs) { return RGB(r + rhs.r, g + rhs.g, b + rhs.b); }
inline RGB& RGB::operator*=(double rhs) { r = r * rhs; g = g * rhs; b = b * rhs; return *this; }};
class Object { public: virtual double testHit(Vector &rayOrigin, Vector &ray, Vector &normal) = 0; virtual RGB getColor() = 0; };
class Plane : public Object { public: Vector normal; Vector pointOnPlane; RGB color;  Plane(int r, int g, int b, double nx, double ny, double nz, double px, double py, double pz) : color(r, g, b), normal(nx, ny, nz), pointOnPlane(px, py, pz) {}
RGB getColor() { return color; };
double testHit(Vector &rayOrigin, Vector &ray, Vector &normal) { //for hit function test depth <= 0 and also the abs(depth) against epsilon value to check for precision errors
	normal = this->normal; return (abs((((pointOnPlane - rayOrigin) * normal) / (ray * normal) > 0) ? ((pointOnPlane - rayOrigin) * normal) / (ray * normal) : 0) > epsilon) ? ((pointOnPlane - rayOrigin) * normal) / (ray * normal) : DBL_MAX;}};
class Sphere : public Object { public: RGB color; Vector origin; double radius; Sphere(int r, int g, int b, double ox, double oy, double oz, double rad) : color(r, g, b), origin(ox, oy, oz), radius(rad) {}
RGB getColor() { return color; };
double testHit(Vector &rayOrigin, Vector &ray, Vector &normal) {
	double a(ray * ray), b(((rayOrigin - origin) * ray) * 2), c(((rayOrigin - origin) * (rayOrigin - origin)) - (radius * radius));
	double discriminant = (b * b) - (4 * a * c);
	double depth = ((-b + sqrt(((b * b) - (4 * a * c)))) < (-b - sqrt(((b * b) - (4 * a * c)))) / (2 * a)) ? (-b + sqrt(((b * b) - (4 * a * c)))) / (2 * a) : (-b - sqrt(((b * b) - (4 * a * c)))) / (2 * a);
	normal = ((rayOrigin + (ray * depth)) - origin).normalize();
	return (abs((depth <= 0 || ((b * b) - (4 * a * c)) < 0) ? 0 : depth) > epsilon) ? depth : DBL_MAX;}}; //test discriminant then test abs(depth) against epsilon value to account for precision errors
int main() { vector<Object *> objects;
	objects.push_back(new Sphere(10, 120, 10, 0, 0, 10, 4));
	objects.push_back(new Plane(200, 200, 200, 0, 1, 0, 0, -4, 3));
	Vector pointLight(15, 15, 0), objectNormal(0, 0, 0);
	srand(NULL); //seed random number generator for random samples
	ofstream file("test.ppm", ios::out); //create ppm file for output image
	file << "P3\n" << imagewidth << " " << imageheight << "\n" << 255 << "\n"; //P3 for ppm ascii representation of rgb
	for (int h = imageheight - 1; h > 0; h--) {	//main loop
		if (!(h % 50)) { cout << (imageheight - double(h)) / imageheight * 100 << " percent done" << endl; }
		for (int w = 0; w < imagewidth; w++) { //calculate primary ray with respect to aspect ratio and field of view
			vector<Vector> samples; // used stochastic stratified sampling (jittered sampling)
			for (double i = 0; i < pow(int(sqrt(imageSamples)), 2); i++) {//16 jittered samples on 1x1 square centered at pixel origin [-0.5, 0.5] ^ 2
				samples.push_back(Vector(((((int(i) % int(sqrt(imageSamples)) / sqrt(imageSamples)) + (((double)rand() / (double)(RAND_MAX)) / sqrt(imageSamples))) * 2) - 1) * 0.5,
					(((((i - (int(i) % int(sqrt(imageSamples)))) / int(sqrt(imageSamples)) / sqrt(imageSamples)) + (((double)rand() / (double)(RAND_MAX)) / sqrt(imageSamples))) * 2) - 1) * 0.5, 0)); }
			RGB output(0, 0, 0), diffuse(0, 0, 0), specular(0, 0, 0), ambient(0, 0, 0);
			for (unsigned int s = 0; s < samples.size(); s++) {
				double px = (tan((90.0 / 2) * 3.141592 / 180) * (double(imagewidth) / imageheight) * ((w + samples[s].x) - (imagewidth / 2.0) + 0.5)) / (imagewidth / 2.0);
				double py = (tan((90.0 / 2) * 3.141592 / 180) * ((h + samples[s].y) - (imageheight / 2.0) + 0.5)) / (imageheight / 2.0); //calculate px and py and add in the sample[i] amount
				double d(DBL_MAX), minDepth(DBL_MAX);
				Object *tempObj = NULL;
				for (unsigned int o = 0; o < objects.size(); o++) { //perform depth testing, find closest object
					Vector tempNormal(0, 0, 0);
					d = objects[o]->testHit(Vector(0, 0, 0), Vector(px, py, 1).normalize(), tempNormal);
					if (d < minDepth) { minDepth = d;
						objectNormal = tempNormal; tempObj = objects[o]; }}
				if (tempObj != NULL) { //now that we have the closest intersected object, shade it
					for (int a = 0; a < ambientSamples; a++) { //generate ambient samples and apply ambient occlusion
						bool hit = false;
						for (unsigned int ainter = 0; ainter < objects.size(); ainter++) {
							Vector sampleV((double)rand() / (double)(RAND_MAX), (double)rand() / (double)(RAND_MAX),0); //random vector in range [0, 1]^2
							double z = sqrt(1.0 - pow(pow((1.0 - sampleV.y), 1.0 / (1 + 1.0)), 2)) * cos(2.0 * PI * sampleV.x);
							double x = sqrt(1.0 - pow(pow((1.0 - sampleV.y), 1.0 / (1 + 1.0)), 2)) * sin(2.0 * PI * sampleV.x);
							double y = pow((1.0 - sampleV.y), 1.0 / (1 + 1.0));
							Vector up(objectNormal), right(1, 0, 0), forward(0, 0, 1);
							if (!(up.x == 0 && up.y == 1 && up.z == 0)) { //create orthonormal base at point of intersection to shoot ambient rays off of
								right = Vector(0, 1, 0) ^ up;
								forward = up ^ right;}
							Vector ray((x*right.x)+(y*up.x)+(z*forward.x), (x*right.y)+(y*up.y)+(z*forward.y), (x*right.z)+(y*up.z)+(z*forward.z));//transform ambient ray
							if (objects[ainter]->testHit((Vector(px, py, 1).normalize() * minDepth), ray, Vector(0, 0, 0)) < DBL_MAX) { hit = true; }}
						if (!hit) {ambient += tempObj->getColor() * ambCoef;}}
					ambient /= ambientSamples;
					for (int shadow = 0; shadow < lightSamples; shadow++) { //create soft shadow from point light, //fake soft shadow by moving around origin of light for every light sample
					Vector random(((((double)rand() / (double)(RAND_MAX))*2)-1)*plr, ((((double)rand() / (double)(RAND_MAX)) * 2) - 1)*plr, ((((double)rand() / (double)(RAND_MAX)) * 2) - 1)*plr);
					bool hit = false;
					for (unsigned int objshd = 0; objshd < objects.size(); objshd++) { //hit calculations to gather how much light from point light is hitting pixel
						if (objects[objshd]->testHit((Vector(px, py, 1).normalize() * minDepth), (((pointLight + random)-(Vector(px, py, 1).normalize() * minDepth)).normalize()), Vector(0, 0, 0)) < DBL_MAX) { hit = true; }}
					if (!hit) {
						Vector wi(((pointLight + random)-(Vector(px, py, 1).normalize() * minDepth)).normalize()), wo((Vector(0, 0, 0) - (Vector(px, py, 1).normalize() * minDepth)).normalize());
						diffuse += tempObj->getColor() * dotproductmax((wi * objectNormal), 0) * lightIntensity; //diffuse and specular shading (phong model)
						specular += RGB(1, 1, 1) * specularReflectionCoefficient * (pow(abs(dotproductmax(((((wi * -1) + (objectNormal * ((objectNormal * wi) * 2)))) * wo), 0)), phongExponent) * dotproductmax((wi * objectNormal), 0)) * lightIntensity;}}
					diffuse /= lightSamples; specular /= lightSamples;}
				output += ambient + diffuse + specular;//combine ambient diffuse and specular shading for phong illumination model
				diffuse = specular = RGB(0, 0, 0);}
			output /= samples.size();
			if (fmax(fmax(output.r, output.g), output.b) > 255) { output /= fmax(fmax(output.r, output.g), output.b); output *= 255;}//perform crude clamping on out of range rgb values
			file << int(output.r) << " " << int(output.g) << " " << int(output.b) << " ";
		} file << "\n"; }
	for (unsigned int v = 0; v < objects.size(); v++) { delete objects[v]; }
	file.close();
	return 0;}