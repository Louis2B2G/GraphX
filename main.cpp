#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "math.h"
#include <string>
#include <iostream>
#include <vector> 
#include <stdlib.h> 
#include <time.h>
#include <chrono> 
using namespace std::chrono; 
using namespace std;

// Vector Class Definition
class Vector { 
    private:
        double coords[3];

    public:
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

        double norm(){
            return sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]);
        }

        Vector& normalize() {
            double norm = this->norm();
            coords[0] = coords[0] / norm;
            coords[1] = coords[1] / norm;
            coords[2] = coords[2] / norm;
            return *this ;
        };
};

// Vector Operators

// Addition
Vector operator+(const Vector& v, const Vector& w){
    return Vector(v[0] + w[0], v[1] + w[1], v[2] + w[2]);
};

//Substraction
Vector operator-(const Vector& v, const Vector& w){
    return Vector(v[0] - w[0], v[1] - w[1], v[2] - w[2]);
};

// Scalar Product
double dot(const Vector& v, const Vector& w){
    return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
};

// Croos product
Vector cross(const Vector& v, const Vector& w){
    return Vector(v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0]);
};

// Pointwise multiplication
Vector mult(const Vector& v, const Vector& w){
    return Vector(v[0]*w[0], v[1]*w[1], v[2]*w[2]);
};

// Norm
double norm(const Vector& v){
    return sqrt(dot(v, v));
};

// Multiplication by a scalar
Vector operator*(const Vector& v, double lambda) {
    return Vector(v[0]*lambda, v[1]*lambda, v[2]*lambda);
};
Vector operator*(double lambda, const Vector& v) {
    return Vector(v[0]*lambda, v[1]*lambda, v[2]*lambda);
};

// printing
ostream& operator<<(ostream &out, Vector const& v) {
    out << "Vector(" << v[0] << "," << v[1] << "," << v[2] << ")";
    return out;
};

// Intersection structure 
struct Intersection{
    bool intersection_exists; // bool saying if there's an intersection
    Vector P; // Closest point of intersection
    Vector N; // vector normal to point of intersection
};

// Ray Class Definition
class Ray { 
    public:
        Vector O; // Origin
        Vector u; // Direction (unit vector)

        explicit Ray(Vector origin, Vector direction){
            direction.normalize();
            O = origin;
            u = direction;
        };

        // default constructor if no arguments are given
        Ray(){
            O = Vector();
            u = Vector(); u.normalize();
        }
};

//printing
ostream& operator<<(ostream &out, Ray const& ray) {
    out << "Ray(" << ray.O << "," << ray.u << ")";
    return out;
};

// Sphere Class Definition
class Sphere { 
    public:
        Vector C; // center
        double R; // radius
        Vector albedo; // rgb color of the sphere
        bool is_mirror; // whether it is a mirror
        int ID; // identification number
        // double n; refraction index

        explicit Sphere(Vector center, double r, Vector color, bool is_a_mirror, int id){ 
            C = center;
            R = r;
            albedo = color;
            is_mirror = is_a_mirror;
            ID = id; 
        };
        
        // default constructor if no arguments are offered
        Sphere(){
            C = Vector(0,0,0);
            R = 0; 
            albedo = Vector(0,0,0); 
            is_mirror = false; 
            ID = 0;
        }

        // Where does the sphere intersect with a given ray?
        Intersection intersection(Ray ray){
            struct Intersection inter;
            double discr = pow(dot(ray.u, ray.O - C), 2) - (dot(ray.O - C, ray.O -C) - pow(R,2));
            inter.intersection_exists = false;

            // if the discriminant is negative, there are no intersections
            if (discr < 0) return inter;
            
            double t1 = dot(ray.u, C - ray.O) + sqrt(discr);
            double t2 = dot(ray.u, C - ray.O) - sqrt(discr);

            // we only consider the root with t > 0
            if (std::max(t1,t2) < 0.) return inter;
        
            double t = std::min(t1,t2); // keep the smallest positive root (we want the closest one)
            inter.intersection_exists = true;
            Vector P = ray.O + t*ray.u;

            // Find the normal vector to the sphere at P
            Vector N = P - C ; N.normalize(); 

            inter.P = P;
            inter.N = N; 

            return inter;
        };
};

// Printing
ostream& operator<<(ostream &out, Sphere const& s) {
    out << "Sphere(ID=" << s.ID <<")";
    return out;
};

struct Scene_Intersection : Intersection{
    Sphere s; // the object being intersected
};

// Scene Class Definition
class Scene { 
    public:
        std::vector<Sphere> spheres;

        explicit Scene(std::vector<Sphere>& sphere_list){
            spheres = sphere_list;
        };

        Scene_Intersection intersection(Ray ray){
            struct Scene_Intersection intersection;
            intersection.intersection_exists = false; // default
            double smallest_distance = std::numeric_limits<double>::infinity(); //initalize as infinity
            Sphere closest_sphere = Sphere();
            Vector P; //closest point
            Vector N; 

            for (auto s : spheres){
                Intersection inter = s.intersection(ray);
                if (inter.intersection_exists){
                    intersection.intersection_exists = true;
                    double distance = norm(ray.O - inter.P); 
                    if (distance < smallest_distance){
                        smallest_distance = distance;
                        P = inter.P;
                        N = inter.N;
                        closest_sphere = s; 
                    }
                }
            }

            intersection.P = P;
            intersection.N = N;
            intersection.s = closest_sphere; 
            return intersection; 
        }
};

ostream& operator<<(ostream &out, Scene const& scene) {
    out << "Scene{";
    for (std::vector<Sphere>::const_iterator i = scene.spheres.begin(); i != scene.spheres.end(); ++i)
        out << *i << ',' << "\n";
    out << "}\n";
    return out;
};

Vector pixel_to_coordinate(Vector camera, double x, double y, double dist_from_screen, int H, int W){
    return Vector(camera[0] + x + 0.5 - W/2, camera[1] + y + 0.5 - H/2, camera[2] - dist_from_screen);
};

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

Vector getColor(Ray& ray, Scene& scene, Vector& light_source, double& I, int max_depth){
    if (max_depth == 0){
        return Vector(0,0,0); // if maximum recusive depth is reached
    }

    // Find the intersection with the closest sphere, s, and get the point, P, and the vector normal to it, N
    Scene_Intersection scene_intersection = scene.intersection(ray);
    Sphere s = scene_intersection.s;
    Vector P = scene_intersection.P; 
    Vector N = scene_intersection.N;

    double eps = pow(10, -4); P  = P + eps*N; // Offset P (in case of noise)
    
    // Check if the the closest sphere is a mirror or not
    if (s.is_mirror){
        Vector reflected_dir = ray.u - 2*(dot(ray.u, N)*N);
        Ray reflected_ray = Ray(P, reflected_dir);
        // find the color of the object that the reflected ray ends up hitting
        return getColor(reflected_ray, scene, light_source, I, max_depth-1);
    }
    else {   
        // DIRECT LIGHTING
        // Find omega and the light ray. Omega is the vector pointing from P in the direction of the light source. d is the distance from P to the light source
        Vector omega = light_source - P ; omega.normalize();
        Ray omega_ray = Ray(P, omega); 
        double d = norm(light_source-P); 

        Vector L0; 

        // if omega_ray has no intersection closer than ||light_source-P||, i.e. nothing is blocking it from reaching the light source
        if (norm(scene.intersection(omega_ray).P - P) > d){
            // return the color of the pixel
            L0 =  I / (4*M_PI*pow(d, 2)) * (1/M_PI) * s.albedo * std::max(0., dot(N, omega));
        }
        else {
            // if it is not visible, the pixel is black (shadow)
            L0 = Vector(0,0,0);
        }

        // COMBINE WITH INDIRECT LIGHTING
        // Create  a bunch of random rays (many different paths)
        int num_paths = 5; 
        Ray random_ray;
        for (int i = 0 ; i < num_paths ; i++){
            random_ray = Ray(P, random_cos(N));
            L0 += mult(s.albedo, getColor(random_ray, scene, light_source, I, max_depth-1)); 
        }
    
        return L0;
    }
};

int main(){
    auto start = high_resolution_clock::now(); // to measure execution time
    srand(time(NULL)); // for generating random numbers

    int W = 500 ; int H = 500; // width and height of image
    unsigned char image[W * H * 3]; //image using row major ordering
    
    // Define the walls
    Sphere wall1 = Sphere(Vector(0,-1000,0), 990, Vector(0, 0, 1), false, 1);
    Sphere wall2 = Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0), false, 2);
    Sphere wall3 = Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0), false, 3);
    Sphere wall4 = Sphere(Vector(0, 0, 1000), 940, Vector(0.5, 0, 0.5), false, 4);
    Sphere wall5 = Sphere(Vector(1000, 0, 0), 940, Vector(0.5, 0.6, 0.9), false, 5);
    Sphere wall6 = Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0), false, 6);

    // Define the light source
    Vector light_source = Vector(-10, 20, 40);
    double I = 2*pow(10,5) ; // light intensity

    // gamma for gamma correction
    double gamma = 2.2; 

    // initialize camera location
    Vector camera = Vector(0, 0, 70);

    // Set the horizontal field of view
    double fov = 60.0;

    // Get the distance from the camera to the screen (pixel grid)
    double dist_from_screen = W / (2*tan( (fov * M_PI / 180) / 2));

    // put some objects in the room
    Sphere sphere1 = Sphere(Vector(0,0,0), 10, Vector(1, 1, 1), true, 7);
    Sphere sphere2 = Sphere(Vector(25, 0,0), 10, Vector(1, 1, 1), false, 8);
    Sphere sphere3 = Sphere(Vector(-15, -5, 30), 5, Vector(1, 1, 0), false, 9);

    // Set the scene
    std::vector<Sphere> objects_in_room; 
    objects_in_room.push_back(sphere1);
    objects_in_room.push_back(sphere2);
    objects_in_room.push_back(sphere3);
    objects_in_room.push_back(wall1);
    objects_in_room.push_back(wall2);
    objects_in_room.push_back(wall3);
    objects_in_room.push_back(wall4);
    objects_in_room.push_back(wall5);
    objects_in_room.push_back(wall6);
    Scene scene = Scene(objects_in_room);

    int max_depth = 3; 
    #pragma omp parallel for schedule(dynamic, 1)
    for (int x = 0 ; x < W ; x++){
        for (int y = 0 ; y < H ; y++){
            // For each pixel (x,y)
            Vector direction = pixel_to_coordinate(camera, x, y, dist_from_screen, H, W);
            // find the ray that passes through the pixel, starting from the camera
            Ray ray = Ray(camera, direction);
            // get the pixel's color
            Vector pixel_color = getColor(ray, scene, light_source, I, max_depth);
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
    stbi_write_jpg("room.jpg", W, H, 3, image, W * sizeof(int));

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Rendering time: " << duration.count()*pow(10, -6) << " seconds\n"<< endl; 

    return 0;
};

