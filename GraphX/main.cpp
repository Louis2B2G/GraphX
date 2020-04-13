#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "math.h"
#include <string>
#include <iostream>
#include <vector> 
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


// Intersection structure 
struct Intersection{
    bool intersection_exists; // bool saying if there's an intersection
    std::vector<Vector> points; // list of points of intersection
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
};

// Sphere Class Definition
class Sphere { 
    public:
        Vector C; // center
        double R; // radius
        Vector albedo; // rgb color of the sphere
        bool is_mirror; // whether it is a mirror

        explicit Sphere(Vector center, double r, Vector color, bool is_a_mirror){ 
            C = center;
            R = r;
            albedo = color;
            is_mirror = is_a_mirror;
        };

        // Where does the sphere intersect with a given ray?
        Intersection intersection(Ray ray){
            struct Intersection inter;
            double discr = pow(dot(ray.u, ray.O - C), 2) - (dot(ray.O - C, ray.O -C) - pow(R,2));
            // if the discriminant is negative, there are no intersections
            if (discr < 0){
                inter.intersection_exists = false;
                return inter;
            }

            inter.intersection_exists = true;

            double t1 = dot(ray.u, C - ray.O) - sqrt(discr);
            double t2 = dot(ray.u, C - ray.O) - sqrt(discr);

            Vector P1 = ray.O + t1 * ray.u;
            Vector P2 = ray.O + t2 * ray.u;

            std::vector<Vector> points; 
            points.push_back(P1);

            if (discr != 0){
                points.push_back(P2);
            }
            
            inter.points = points;
            return inter;
        };
};

// Scene Class Definition
class Scene { 
    public:
        std::vector<Sphere> spheres;

        explicit Scene(std::vector<Sphere>& sphere_list){
            spheres = sphere_list;
        };
};

// Printing Functions
ostream& operator<<(ostream &out, Vector const& v) {
    out << "Vector(" << v[0] << "," << v[1] << "," << v[2] << ")";
    return out;
};
ostream& operator<<(ostream &out, Sphere const& s) {
    out << "Sphere(" << s.C << "," << s.R << "," << s.albedo << "," << s.is_mirror << ")";
    return out;
};
ostream& operator<<(ostream &out, Ray const& ray) {
    out << "Ray(" << ray.O << "," << ray.u << ")";
    return out;
};
ostream& operator<<(ostream &out, Scene const& scene) {
    out << "Scene(";
    for (std::vector<Sphere>::const_iterator i = scene.spheres.begin(); i != scene.spheres.end(); ++i)
        std::cout << *i << ',';
    return out;
};

Vector pixel_to_coordinate(Vector camera, double x, double y, double dist_from_screen, int H, int W){
    return Vector(camera[0] + x + 0.5 - W/2, camera[1] + y + 0.5 - H/2, camera[2] - dist_from_screen);
};

// determine if P can be reached by the light source
bool is_visible(Ray light_ray, Vector P, Scene scene){ 
    // gotta loop through the other spheres and check their intersection with light_ray and make sure it's null OR that its further away than P
    for (auto s : scene.spheres){
        struct Intersection inter = s.intersection(light_ray);
        if (inter.intersection_exists){
            // look through all of the points of intersection
            for (auto Q : inter.points){
                // check if the point of intersection is closer to the light source than P
                if (norm(light_ray.O - Q) < norm(light_ray.O - P)){
                    // make sure it's not in the wrong direction (on the other side of the ray)
                    if (norm(Q - P) < norm(light_ray.O - P)){
                        return false;
                    }
                }
            }
        }
    }
    return true;
};

Vector getColor(Ray& ray, Scene& scene, Vector& light_source, double& I, int max_depth){

    if (max_depth == 0){
        return Vector(1,1,1); // if maximum recusive depth is reached
    }

    // Find the closest sphere
    double smallest_distance_to_pixel =  std::numeric_limits<double>::infinity(); //initalize as infinity
    Sphere closest_sphere = Sphere(Vector(0,0,0), 0, Vector(0,0,0), false);
    Vector P; //closest point on a sphere

    // loop through all spheres / intersections of those spheres with the ray and find the closest one that's
    // in front of the camera
    for (auto s: scene.spheres){
        struct Intersection inter = s.intersection(ray);
        for (auto Q : inter.points){
            // check if the point of intersection is closer to the source than P
            if (norm(ray.O - Q) < smallest_distance_to_pixel){
                // check that it's not behind the camera
                Vector t = (ray.O - Q) ; t.normalize();
                if (norm(ray.u - t) == 0){
                    closest_sphere = s; 
                    P = Q; 
                    smallest_distance_to_pixel = norm(ray.O - Q);
                }
            }
        } 
    }

    // Find the normal vector to the sphere at P
    Vector N = P - closest_sphere.C ; N.normalize(); 

    // Find omega and the light ray. Omega is the vector pointing from the light source in the direction of P
    Vector omega = light_source - P ; omega.normalize();
    Ray light_ray = Ray(light_source, omega);

    double eps = 0.01; 

    // Check if the the closest sphere is a mirror or not
    if (closest_sphere.is_mirror){
        Vector reflected_dir = ray.u - 2*dot(ray.u, N)*N;
        P += eps*N;
        Ray reflected_ray = Ray(P, reflected_dir);
        // find the color of the object that the reflected ray ends up hitting
        if (max_depth < 5){
            cout << "Closest sphere:" << closest_sphere << "\n";
        }
        return getColor(reflected_ray, scene, light_source, I, max_depth-1);
    }
    else {    
        if (max_depth < 5){
            cout << "Closest sphere:" << closest_sphere << "\n \n \n";
        }
        P += eps*N; // Offset P (in case of noise)
        Vector lambertian_reflection_intensity = I / (4*M_PI*dot(light_source-P, light_source-P)) * (1/M_PI) * closest_sphere.albedo * is_visible(light_ray, P, scene) * std::max(0., dot(N, omega));
        return lambertian_reflection_intensity;
    }
};

int main(){
    int W = 500 ; int H = 500; // width and height of image
    unsigned char image[W * H * 3]; //image using row major ordering

    /*
    Red pixel value at (x,y):   image[y*W*3 + x*3 + 0]
    Green pixel value at (x,y): image[y*W*3 + x*3 + 1]
    Blue pixel value at (x,y):  image[y*W*3 + x*3 + 2]
    */

    // Define the walls
    Sphere wall1 = Sphere(Vector(0,-1000,0), 990, Vector(0, 0, 1), false);
    Sphere wall2 = Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0), false);
    Sphere wall3 = Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0), false);
    Sphere wall4 = Sphere(Vector(0, 0, 1000), 940, Vector(0.5, 0, 0.5), false);

    // Define the light source
    Vector light_source = Vector(-10, 50, 40);
    double I = pow(10,5) ; // light intensity

    // gamma for gamma correction
    double gamma = 2.9; 

    // initialize camera at (0,0,55)
    Vector camera = Vector(0, 0, 55);

    // Set the horizontal field of view to 60 degrees
    double fov = 60.0;

    // Get distance from the camera to the screen (pixel grid)
    double dist_from_screen = W / (2*tan( (fov * M_PI / 180) / 2));

    // put a sphere in the middle of the room at (0, 0, 0)
    Sphere sphere = Sphere(Vector(0,0,0), 10, Vector(1, 1, 1), true);

    // Set the scene
    std::vector<Sphere> objects_in_room; 
    objects_in_room.push_back(sphere);
    objects_in_room.push_back(wall1);
    objects_in_room.push_back(wall2);
    objects_in_room.push_back(wall3);
    objects_in_room.push_back(wall4);
    Scene scene = Scene(objects_in_room);

    for (int x = 0 ; x < W ; x++){
        for (int y = 0 ; y < H ; y++){
            // For each pixel (x,y)
            Vector direction = pixel_to_coordinate(camera, x, y, dist_from_screen, H, W);
            // find the ray that passes through the pixel
            Ray ray = Ray(camera, direction);
        
            // get the pixel's color
            int max_depth = 5; 
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

    return 0;
};

