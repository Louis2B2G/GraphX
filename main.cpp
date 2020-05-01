/* ------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------- 

CSE 306 - Computer Graphics : Louis de Benoist
Rendering using Ray-Tracing

------------------------------------------------------------------------------------------------------------------------ 
------------------------------------------------------------------------------------------------------------------------ */

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "math.h"
#include <string>
#include <iostream>
#include <vector> 
#include <stdlib.h> 
#include <time.h>
#include <algorithm> 
#include <chrono> 
#include <fstream>
using namespace std::chrono; 

/* ------------------------------------------------------------------------------------------------------------------------

Vector Class Definition

------------------------------------------------------------------------------------------------------------------------ */

class Vector { 
    private:
        double coords[3];
    public:
        ~Vector() {}
        explicit Vector(double x = 0. , double y = 0. , double z = 0.){ 
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        const double& operator[](int i) const{
            return coords[i];
        };
        double& operator[](int i){
            return coords[i];
        };
        Vector& operator+=(const Vector& v) {
            coords[0] += v[0]; 
            coords[1] += v[1]; 
            coords[2] += v[2]; 
            return *this ;
        };
        Vector& operator-=(const Vector& v) {
            coords[0] -= v[0]; 
            coords[1] -= v[1]; 
            coords[2] -= v[2]; 
            return *this ;
        };
        Vector& normalize() {
            double norm = sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]);
            coords[0] = coords[0] / norm;
            coords[1] = coords[1] / norm;
            coords[2] = coords[2] / norm;
            return *this ;
        };
        std::string to_string(){
           return "Vector(" + std::to_string(coords[0]) + "," + std::to_string(coords[1]) + "," + std::to_string(coords[2]) + ")";
        }
};


/* ------------------------------------------------------------------------------------------------------------------------

Vector Operators
 - Addition
 - Substraction
 - Scalar Product
 - Cross Product
 - Pointwise Multiplication
 - Norm
 - Multiplication by a Scalar

------------------------------------------------------------------------------------------------------------------------ */

Vector operator+(const Vector& v, const Vector& w){
    return Vector(v[0] + w[0], v[1] + w[1], v[2] + w[2]);
};
Vector operator-(const Vector& v, const Vector& w){
    return Vector(v[0] - w[0], v[1] - w[1], v[2] - w[2]);
};
double dot(const Vector& v, const Vector& w){
    return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
};
Vector cross(const Vector& v, const Vector& w){
    return Vector(v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0]);
};
Vector mult(const Vector& v, const Vector& w){
    return Vector(v[0]*w[0], v[1]*w[1], v[2]*w[2]);
};
double norm(const Vector& v){
    return sqrt(dot(v, v));
};
Vector operator*(const Vector& v, double lambda) {
    return Vector(v[0]*lambda, v[1]*lambda, v[2]*lambda);
};
Vector operator*(double lambda, const Vector& v) {
    return Vector(v[0]*lambda, v[1]*lambda, v[2]*lambda);
};


/* ------------------------------------------------------------------------------------------------------------------------

Ray Class Definition

------------------------------------------------------------------------------------------------------------------------ */

class Ray { 
    public:
        Vector O; // Origin
        Vector u; // Direction (unit vector)

        explicit Ray(Vector origin, Vector direction){
            //direction.normalize(); //for speed, assume that direction is already normalized
            O = origin;
            u = direction;
        };

        ~Ray() {}

        // default constructor if no arguments are given
        Ray(){
            O = Vector();
            u = Vector();
        }

        std::string to_string(){
            return "Ray(" + O.to_string() + "," + u.to_string() + ")";
        }
};

/* ------------------------------------------------------------------------------------------------------------------------

Object Intersection Structure

------------------------------------------------------------------------------------------------------------------------ */

struct Intersection{
    bool exists; // bool saying if there's an intersection
    Vector P; // Closest point of intersection
    Vector N; // vector normal to point of intersection
    double t; // scaling to get to the intersection
};


/* ------------------------------------------------------------------------------------------------------------------------

Plane Class

------------------------------------------------------------------------------------------------------------------------ */

class Plane{
public:
    Vector A;
    Vector N; 

    Plane(Vector A_, Vector N_){
        A = A_;
        N = N_;
    }

    Intersection intersection(Ray ray){
        double t = dot(A-ray.O, N) / dot(ray.u, N);
        Intersection inter; inter.exists = false; 
        if (t > 0){
            inter.exists = true; 
            inter.t = t; 
        }
        return inter; 
    }
};


/* ------------------------------------------------------------------------------------------------------------------------

Geometry Class Definition

------------------------------------------------------------------------------------------------------------------------ */

class Geometry{
    public:
        virtual Intersection intersection(Ray ray) = 0;
        virtual std::string to_string() = 0; 
        bool is_mirror; 
        bool is_mesh;
        Vector albedo; 
        int ID; 
};  


/* ------------------------------------------------------------------------------------------------------------------------

Sphere Class Definition

------------------------------------------------------------------------------------------------------------------------ */

class Sphere : public Geometry { 
    public:
        Vector C; // center
        double R; // radius

        explicit Sphere(Vector center, double r, Vector color, std::string material, int id){ 
            C = center;
            R = r;
            albedo = color;
            if (material == "mirror"){is_mirror = true;} else {is_mirror = false;}
            is_mesh = false;
            ID = id; 
        };
        
        ~Sphere() {}
        // default constructor if no arguments are offered
        Sphere(){
            C = Vector(0,0,0);
            R = 0; 
            albedo = Vector(0,0,0); 
            is_mirror = false; 
            ID = -1;
        }

        // Where does the sphere intersect with a given ray?
        Intersection intersection(Ray ray){
            struct Intersection inter;
            double discr = pow(dot(ray.u, ray.O - C), 2) - (dot(ray.O - C, ray.O -C) - pow(R,2));
            inter.exists = false;

            // if the discriminant is negative, there are no intersections
            if (discr < 0) return inter;
            
            double t1 = dot(ray.u, C - ray.O) + sqrt(discr);
            double t2 = dot(ray.u, C - ray.O) - sqrt(discr);

            // we only consider the root with t > 0
            if (std::max(t1,t2) < 0.) return inter;
        
            double t = std::min(t1,t2); // keep the smallest positive root (we want the closest one)
            inter.exists = true;
            Vector P = ray.O + t*ray.u;

            // Find the normal vector to the sphere at P
            Vector N = P - C ; N.normalize(); 

            inter.P = P;
            inter.N = N; 

            return inter;
        };

        std::string to_string(){
            return "Sphere(ID=" + std::to_string(ID) + ")"; 
        };
};

/* ------------------------------------------------------------------------------------------------------------------------

Bounding Box Class Definition

------------------------------------------------------------------------------------------------------------------------ */

class BoundingBox {
public:
    Vector Bmin; 
    Vector Bmax;

    BoundingBox(Vector B_min, Vector B_max){
        Bmin = B_min;
        Bmax = B_max;
    };

    // returns whether the ray intersects the bounding box
    bool intersection(Ray ray){
        // Planes where x is constant
        Intersection inter1 = Plane(Bmin, Vector(1,0,0)).intersection(ray); 
        Intersection inter2 = Plane(Bmax, Vector(1,0,0)).intersection(ray);

        // Planes where y is constant
        Intersection inter3 = Plane(Bmin, Vector(0,1,0)).intersection(ray);
        Intersection inter4 = Plane(Bmax, Vector(0,1,0)).intersection(ray);

        // Planes where z is constant
        Intersection inter5 = Plane(Bmin, Vector(0,0,1)).intersection(ray);
        Intersection inter6 = Plane(Bmax, Vector(0,0,1)).intersection(ray);

        double tx_0 = - std::numeric_limits<double>::infinity(); 
        double tx_1 = std::numeric_limits<double>::infinity(); 
        double ty_0 = - std::numeric_limits<double>::infinity(); 
        double ty_1 = std::numeric_limits<double>::infinity(); 
        double tz_0 = - std::numeric_limits<double>::infinity(); 
        double tz_1 = std::numeric_limits<double>::infinity(); 
        
        if (inter1.exists && inter2.exists){
            tx_0 = std::min(inter1.t, inter2.t);
            tx_1 = std::max(inter1.t, inter2.t); 
        };
        if (inter3.exists && inter4.exists){
            ty_0 = std::min(inter3.t, inter4.t);
            ty_1 = std::max(inter3.t, inter4.t); 
        };
        if (inter5.exists && inter6.exists){
            tz_0 = std::min(inter5.t, inter6.t);
            tz_1 = std::max(inter5.t, inter6.t); 
        };

        double l0[] = {tx_0, ty_0, tz_0};
        double l1[] = {tx_1, ty_1, tz_1};
        if (*std::min_element(l1,l1+3) > *std::max_element(l0,l0+3)){return true;}
        return false;
    };
};


/* ------------------------------------------------------------------------------------------------------------------------

TriangleIndices and TriangleMesh Class Definitions

------------------------------------------------------------------------------------------------------------------------ */

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class TriangleMesh : public Geometry {
public:
    ~TriangleMesh() {}
	TriangleMesh() {
        is_mesh = true; 
    };

    void resize(double scaling, Vector translation){
        for (int i = 0; i < vertices.size(); i++){
            vertices[i] = scaling*vertices[i] + translation;
        }
    }

    std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

    std::string to_string(){
        // print the vertices
        return "TriangleMesh with " + std::to_string(vertices.size()) + " vertices\n";
    };

    // return bounding box for the mesh
    BoundingBox getBoundingBox(){
       double min_x = std::numeric_limits<double>::infinity(); 
       double min_y = std::numeric_limits<double>::infinity(); 
       double min_z = std::numeric_limits<double>::infinity(); 

       double max_x = -std::numeric_limits<double>::infinity();
       double max_y = -std::numeric_limits<double>::infinity();
       double max_z = -std::numeric_limits<double>::infinity();

       for (auto v : vertices){
           if (v[0] < min_x) min_x = v[0];
           if (v[0] > max_x) max_x = v[0];
           if (v[1] < min_y) min_y = v[1];
           if (v[1] > max_y) max_y = v[1];
           if (v[2] < min_z) min_z = v[2];
           if (v[2] > max_z) max_z = v[2];
       }
       Vector Bmin = Vector(min_x, min_y, min_z);
       Vector Bmax = Vector(max_x, max_y, max_z);
       return BoundingBox(Bmin, Bmax);
    }

    Intersection intersection(Ray ray){
        Intersection inter; inter.exists = false; 

        // first, check if the ray intersects the bounding box
        BoundingBox box = getBoundingBox();
        if (box.intersection(ray) == false){
            return inter; 
        }

        double smallest_t = std::numeric_limits<double>::infinity(); //initalize as infinity
        TriangleIndices intersected_triangle; 
        double alpha_intersected; double beta_intersected; double gamma_intersected; 

        for (auto triangle : indices){
            // vertices of the triangle associated to the index
            Vector A = vertices[triangle.vtxi]; 
            Vector B = vertices[triangle.vtxj]; 
            Vector C = vertices[triangle.vtxk]; 

            Vector e1 = B-A ; Vector e2 = C-A; 
            Vector non_normalized_normal = cross(e1, e2);

            // compute t
            double t = dot(A-ray.O, non_normalized_normal) / dot(ray.u, non_normalized_normal);
            if (t < smallest_t && t > 0){
                // compute beta, gamma, and alpha
                double beta = dot(e2, cross(A-ray.O, ray.u)) / dot(ray.u, non_normalized_normal);
                double gamma = - dot(e1, cross(A-ray.O, ray.u)) / dot(ray.u, non_normalized_normal);
                double alpha = 1 - beta - gamma; 

                // compute P
                inter.P = A + beta*e1 + gamma*e2;

                // check if P is inside the triangle
                if ((alpha + beta + gamma == 1)  && alpha > 0 && beta > 0 && gamma > 0){
                    inter.exists = true; 
                    alpha_intersected = alpha ; beta_intersected = beta; gamma_intersected = gamma; 
                    smallest_t = t; 
                    intersected_triangle = triangle;
                }
            }
        }

        // compute the shading normal (Phong interpolation)
        Vector N_A = normals[intersected_triangle.ni]; 
        Vector N_B = normals[intersected_triangle.nj]; 
        Vector N_C = normals[intersected_triangle.nk]; 

        inter.N = alpha_intersected*N_A + beta_intersected*N_B + gamma_intersected*N_C; 
        inter.N = inter.N.normalize();

        return inter; 
    }
	
	void readOBJ(const char* filename) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(filename, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}
			}
		}
		fclose(f);
	}
};


/* ------------------------------------------------------------------------------------------------------------------------

Scene Intersection Structure

------------------------------------------------------------------------------------------------------------------------ */

struct Scene_Intersection : Intersection{
    Geometry* object; // the object being intersected
};


/* ------------------------------------------------------------------------------------------------------------------------

Scene Class Definition

------------------------------------------------------------------------------------------------------------------------ */

class Scene { 
    public:
        std::vector<Geometry*> objects;

        explicit Scene(std::vector<Geometry*>& objects_list){
            objects = objects_list;
        };

        Scene_Intersection intersection(Ray ray){
            struct Scene_Intersection scene_intersection;
            scene_intersection.exists = false; // default
            double smallest_distance = std::numeric_limits<double>::infinity(); //initalize as infinity
            Vector P; //closest point
            Vector N; 
            Geometry* closest_object;

            for (int i = 0; i < objects.size(); i++){
                Intersection inter = objects[i] -> intersection(ray);
                if (inter.exists){
                    scene_intersection.exists = true;
                    double distance = norm(ray.O - inter.P); 
                    if (distance < smallest_distance){
                        smallest_distance = distance;
                        P = inter.P;
                        N = inter.N;
                        closest_object = objects[i]; 
                    }
                }
            }
            scene_intersection.P = P;
            scene_intersection.N = N;
            scene_intersection.object = closest_object; 
            return scene_intersection; 
        }

        std::string to_string(){
            std::string res = "";
            for (int i = 0; i < objects.size(); i++){
                res += objects[i] -> to_string();
            }
            res += "}\n"; return res;
        };
};


/* ------------------------------------------------------------------------------------------------------------------------

Get the coordinates of a pixel

------------------------------------------------------------------------------------------------------------------------ */

Vector pixel_to_coordinate(Vector& camera, double x, double y, double dist_from_screen, int H, int W){
    return Vector(camera[0] + x + 0.5 - W/2, camera[1] + y + 0.5 - H/2, camera[2] - dist_from_screen);
};


/* ------------------------------------------------------------------------------------------------------------------------

Get a random vector shooting off a point on a surface with a normal vector N

------------------------------------------------------------------------------------------------------------------------ */

Vector random_cos(Vector N){
    // get r1,r2 ~ U([0,1])
    double r1 = ((double) rand() / (RAND_MAX));
    double r2 = ((double) rand() / (RAND_MAX));

    double x = cos(2*M_PI*r1)*sqrt(1-r2);
    double y = sin(2*M_PI*r1)*sqrt(1-r2);
    double z = sqrt(r2);

    Vector T1;
    if ((pow(N[0], 2) < pow(N[1], 2)) && (pow(N[0], 2) < pow(N[2], 2))){
        T1[0] = 0;
        // swap the values of the other 2 and negate one of them
        T1[1] = N[2];
        T1[2] = -1*N[1];  
    }
    else if ((pow(N[1], 2) < pow(N[0], 2)) && (pow(N[1], 2) < pow(N[2], 2))){
        T1[1] = 0;
        // swap the values of the other 2 and negate one of them
        T1[0] = N[2];
        T1[2] = -1*N[0];  
    }
    else{
        T1[2] = 0;
        // swap the values of the other 2 and negate one of them
        T1[0] = N[1];
        T1[1] = -1*N[0];  
    }

    T1.normalize();
    Vector T2 = cross(N, T1);
    return x*T1 + y*T2 + z*N; 
}


/* ------------------------------------------------------------------------------------------------------------------------

Sample x and y from N(0, standard_deviation) - to avoid only shooting rays through the center of the pixel

------------------------------------------------------------------------------------------------------------------------ */

void boxMuller(double standard_deviation, double &x, double &y){
    // get r1,r2 ~ U([0,1])
    double r1 = ((double) rand() / (RAND_MAX));
    double r2 = ((double) rand() / (RAND_MAX));

    x = sqrt(-2 * log(r1))*cos(2 * M_PI*r2)*standard_deviation; 
    y = sqrt(-2 * log(r1))*sin(2 * M_PI*r2)*standard_deviation;
}


/* ------------------------------------------------------------------------------------------------------------------------

Gives the color value of a given pixel

------------------------------------------------------------------------------------------------------------------------ */

Vector getColor(Ray& ray, Scene* scene, Vector* light_source, double& I, int max_depth){
    if (max_depth == 0){
        return Vector(0,0,0); // if maximum recusive depth is reached
    }

    // Find the intersection with the closest object and get the intersection point, P, and the vector normal to it, N
    Scene_Intersection scene_intersection = scene->intersection(ray);
    Geometry* object = scene_intersection.object;
    
    Vector P = scene_intersection.P; 
    Vector N = scene_intersection.N;

    double eps = pow(10, -4); P  = P + eps*N; // Offset P (in case of noise)
    
    // Check if the the closest object is a mirror or not
    if (object->is_mirror){
        Vector reflected_dir = ray.u - 2*(dot(ray.u, N)*N);
        Ray reflected_ray = Ray(P, reflected_dir);
        // find the color of the object that the reflected ray ends up hitting
        return getColor(reflected_ray, scene, light_source, I, max_depth-1);
    }
    
    // Case where you have a colored mesh
    /*
    if (object -> is_mesh){


    }
    */

    // if it is just some colored sphere
    else {   
        // DIRECT LIGHTING
        // Find omega and the light ray. Omega is the vector pointing from P in the direction of the light source. d is the distance from P to the light source
        Vector omega = *light_source - P ; omega.normalize();
        Ray omega_ray = Ray(P, omega); 
        double d = norm(*light_source-P); 

        Vector L0; 

        // if omega_ray has no intersection closer than ||light_source-P||, i.e. nothing is blocking it from reaching the light source
        if (norm(scene->intersection(omega_ray).P - P) > d){
            // return the color of the pixel
            L0 =  I / (4*M_PI*pow(d, 2)) * (1/M_PI) * object->albedo * std::max(0., dot(N, omega));
        }
        else {
            // if it is not visible, the pixel is black (shadow)
            L0 = Vector(0,0,0);
        }

        // Combine with indirect lighting
        Ray random_ray = Ray(P, random_cos(N));
        L0 = L0 + mult(object->albedo, getColor(random_ray, scene, light_source, I, max_depth-1)); 
    
        return L0;
    }
};


/* ------------------------------------------------------------------------------------------------------------------------

Main function

------------------------------------------------------------------------------------------------------------------------ */

int main(){
    auto start = high_resolution_clock::now(); // to measure execution time
    srand(time(NULL)); // for generating random numbers

    int W = 500 ; int H = 500; // width and height of image
    unsigned char image[W * H * 3]; //image using row major ordering
    
    // Define the walls (pointers to big spheres)
    Sphere* wall1 = new Sphere(Vector(0,-1000,0), 990, Vector(0, 0, 1),  "none", 0);
    Sphere* wall2 = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0),  "none", 1);
    Sphere* wall3 = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0),  "none", 2);
    Sphere* wall4 = new Sphere(Vector(0, 0, 1000), 940, Vector(0.5, 0, 0.5),  "none", 3);
    Sphere* wall5 = new Sphere(Vector(1000, 0, 0), 940, Vector(0.5, 0.6, 0.9),  "none", 4);
    Sphere* wall6 = new Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0),  "none", 5);

    // Define the light source
    Vector* light_source = new Vector(-10, 20, 40);
    double I = 2*pow(10,5) ; // light intensity

    // gamma for gamma correction
    double gamma = 2.2; 

    // initialize camera location
    Vector camera = Vector(0, 0, 70);

    // Set the horizontal field of view
    double fov = 60.0;

    // Get the distance from the camera to the screen (pixel grid)
    double dist_from_screen = W / (2*tan( (fov * M_PI / 180) / 2));

    // put some spheres in the room
    Sphere* sphere1 = new Sphere(Vector(0,0,0), 10, Vector(1, 1, 1),  "mirror", 6);
    Sphere* sphere2 = new Sphere(Vector(25, 0,0), 10, Vector(1, 1, 1),  "none", 7);
    Sphere* sphere3 = new Sphere(Vector(-15, -5, 30), 5, Vector(1, 1, 0), "none", 8);

    // Load the cat object
    const char* filename = "Cat/cat.obj";
    TriangleMesh* cat = new TriangleMesh();
    cat -> readOBJ(filename);
    cat -> is_mirror = true; 
    cat -> resize(0.5, Vector(0, -10, 0));
 
    // Set the scene
    std::vector<Geometry*> objects_in_room; 
    objects_in_room.push_back(wall1);
    objects_in_room.push_back(wall2);
    objects_in_room.push_back(wall3);
    objects_in_room.push_back(wall4);
    objects_in_room.push_back(wall5);
    objects_in_room.push_back(wall6);
    //objects_in_room.push_back(sphere1);
    //objects_in_room.push_back(sphere2);
    //objects_in_room.push_back(sphere3);
    objects_in_room.push_back(cat);
    Scene* scene = new Scene(objects_in_room);

    int max_depth = 2; 
    int num_paths = 2; 
    double stdv = 0.5; 

    #pragma omp parallel for schedule(dynamic, 1)
    for (int x = 0 ; x < W ; x++){
        for (int y = 0 ; y < H ; y++){
            // get the pixel's color
            Vector pixel_color;

            for (int k = 0 ; k < num_paths ; k++){
                // For each pixel (x,y)
                double eps ; double eta; boxMuller(stdv, eps, eta); // add random fluctuation
                Vector random_direction = pixel_to_coordinate(camera, x+eps, y+eta, dist_from_screen, H, W); random_direction.normalize();

                // find the ray that passes through the pixel, starting from the camera
                Ray ray = Ray(camera, random_direction);
                pixel_color = pixel_color + getColor(ray, scene, light_source, I, max_depth);
            }

            pixel_color = pow(num_paths,-1) * pixel_color; 
            // gamma correction
            pixel_color[0] = pow(pixel_color[0], 1/gamma);
            pixel_color[1] = pow(pixel_color[1], 1/gamma);
            pixel_color[2] = pow(pixel_color[2], 1/gamma);
            // scaling
            pixel_color = 255 * pixel_color; 
            // update the pixel on the image
            image[y*W*3 + x*3 + 0] = std::min(255., pixel_color[0]);
            image[y*W*3 + x*3 + 1] = std::min(255., pixel_color[1]);
            image[y*W*3 + x*3 + 2] = std::min(255., pixel_color[2]);
        }
    }
    // save the image
    stbi_write_jpg("Images/cat.jpg", W, H, 3, image, W * sizeof(int));

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    std::cout << "Rendering time: " << duration.count()*pow(10, -6) << " seconds\n"; 
    return 0;
};
