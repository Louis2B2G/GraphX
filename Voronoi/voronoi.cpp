/* ------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------- 

CSE 306 - Computer Graphics : Louis de Benoist
Project #2

------------------------------------------------------------------------------------------------------------------------ 
------------------------------------------------------------------------------------------------------------------------ */

#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb_image_write.h"
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

const double INF = std::numeric_limits<double>::infinity(); 

/* ------------------------------------------------------------------------------------------------------------------------

Vector Class Definition (redefined to be 2D instead of 3D)

------------------------------------------------------------------------------------------------------------------------ */

class Vector { 
    public:
        ~Vector() {}

        double x;
        double y; 

        explicit Vector(double x_ = 0. , double y_ = 0.){ 
            x = x_;
            y = y_;
        };
        Vector& operator+=(const Vector& v) {
            x += v.x; 
            y += v.y; 
            return *this ;
        };
        Vector& operator-=(const Vector& v) {
            x -= v.x; 
            y -= v.y;
            return *this ;
        };
        Vector& normalize() {
            double norm = sqrt(x*x + y*y);
            x = x / norm;
            y = y / norm;
            return *this ;
        };
        std::string to_string(){
           return "Vector(" + std::to_string(x) + "," + std::to_string(y) + ")";
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
    return Vector(v.x + w.x, v.y + w.y);
};
Vector operator-(const Vector& v, const Vector& w){
    return Vector(v.x - w.x, v.y - w.y);
};
double dot(const Vector& v, const Vector& w){
    return v.x *w.x + v.y*w.y;
};
Vector mult(const Vector& v, const Vector& w){
    return Vector(v.x*w.x, v.y*w.y);
};
double norm(const Vector& v){
    return sqrt(dot(v, v));
};
Vector operator*(const Vector& v, double lambda) {
    return Vector(v.x*lambda, v.y*lambda);
};
Vector operator*(double lambda, const Vector& v) {
    return Vector(v.x*lambda, v.y*lambda);
};

/* ------------------------------------------------------------------------------------------------------------------------

Polygon Class Definition

------------------------------------------------------------------------------------------------------------------------ */

class Polygon{
    public:
        std::vector<Vector> vertices; // set of vertices (Vectors)
        std::vector<std::pair<Vector, Vector>> edges; // set of edges (pairs of Vectors)

        explicit Polygon(std::vector<Vector> v, std::vector<std::pair<Vector, Vector>> e){ 
            vertices = v; 
            edges = e; 
        };

        // default constructor if no arguments are given
        Polygon(){}

};

/* ------------------------------------------------------------------------------------------------------------------------

Intersect function

------------------------------------------------------------------------------------------------------------------------ */

Vector intersect(Vector A, Vector B, std::pair<Vector, Vector> edge){
    Vector u = edge.first, v = edge.second; // extract two points that define the infinite line
    Vector N = Vector(v.y - u.y, u.x - v.x);
    double t = dot(u-A, N) / dot(B-A, N);
    return A + t*(B - A);
}

/* ------------------------------------------------------------------------------------------------------------------------

Inside function

------------------------------------------------------------------------------------------------------------------------ */

bool isInside(Vector P, std::pair<Vector, Vector> edge){
    Vector u = edge.first, v = edge.second; // extract two points that define the infinite line
    Vector N = Vector(v.y - u.y, u.x - v.x);
    
    return dot(P - u, N) <= 0;
}

/* ------------------------------------------------------------------------------------------------------------------------

Sutherland-Hodgman’s polygon clipping algorithm

------------------------------------------------------------------------------------------------------------------------ */

Polygon* clip(Polygon* subjectPolygon, Polygon* clipPolygon){

    Polygon* outPolygon;
    
    // go through each edge of the clipPolygon
    for (auto clipEdge : clipPolygon->edges){
        outPolygon = new Polygon();

        // for each vertex of the subject polygon
        for (int i = 0; i < subjectPolygon->vertices.size(); i++){
            // Test the subject polygon edge with vertices (i-1, i)
            Vector curVertex = subjectPolygon->vertices[i];
            Vector prevVertex = subjectPolygon->vertices[(i > 0)?(i-1):(subjectPolygon->vertices.size()-1)];

            // Compute intersection between the infinite line supported by clipEdge and edge (i-1, i)
            Vector intersection = intersect(prevVertex, curVertex, clipEdge);
            if (isInside(curVertex, clipEdge)){
                if (!isInside(prevVertex, clipEdge)){
                    // The subject polygon edge crosses the clip edge, and we leave the clipping area
                    outPolygon->vertices.push_back(intersection);
                }
                outPolygon->vertices.push_back(curVertex);
            }
            else if (isInside(prevVertex, clipEdge)){
                // The subject polygon edge crosses the clip edge, and we enter the clipping area
                outPolygon->vertices.push_back(intersection);
            }
        }
        subjectPolygon = outPolygon; 
    }
    return outPolygon;
}
/* ------------------------------------------------------------------------------------------------------------------------

Get a bounding box for some data

------------------------------------------------------------------------------------------------------------------------ */

Polygon* getBoundingBox(std::vector<Vector> data){
    std::vector<double> xList, yList;
    for (int i = 0; i < data.size(); i++){
        xList.push_back(data[i].x);
        yList.push_back(data[i].y);
    }
    double smallest_x = *min_element(xList.begin(), xList.end()) - 1; 
    double biggest_x = *max_element(xList.begin(), xList.end()) + 1; 
    double smallest_y = *min_element(yList.begin(), yList.end()) - 1; 
    double biggest_y = *max_element(yList.begin(), yList.end()) + 1; 

    Vector v1 = Vector(smallest_x, smallest_y), v2 = Vector(smallest_x, biggest_y), v3 = Vector(biggest_x, smallest_y), v4 = Vector(biggest_x, biggest_y); 
    std::pair<Vector, Vector> e1 = {v1, v3}, e2 = {v2, v4}, e3 = {v1, v2}, e4 = {v3, v4};

    return new Polygon({v1, v3, v4, v2}, {e1, e2, e3, e4});
}

/* ------------------------------------------------------------------------------------------------------------------------

Voronoi Parallel Linear Enumeration (Power Diagram)
 - takes a bunch of data points (Vectors) and returns a vector of Polygons such that the ith Polygon corresponds 
   to the Voronoi polygon of the ith data point
- Also takes into account weights for each cell (Power Diagram)
------------------------------------------------------------------------------------------------------------------------ */

std::vector<Polygon> getVoronoi(std::vector<Vector> P, std::vector<double> weights){
    std::vector<Polygon> VoronoiDiagram;
    Polygon* boundingBox = getBoundingBox(P); 

    // get the voronoi Polygon for each data point
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < P.size(); i++){
        // initialize the big "quad" that we will then proceed to clip
        Polygon subjectPolygon = *boundingBox; 

        for (int j = 0 ; j < P.size(); j++){
            if (i == j) continue; 
            Vector M = 0.5*(P[i] + P[j]); // vector in the middle of P[i] and P[j]
            M = M + 0.5*(weights[i] - weights[j])/dot(P[i]-P[j], P[i]-P[j]) * (P[j] - P[i]); // power diagram modification

            // apply a modified Sutherland-Hodgman’s polygon clipping algorithm 
            Polygon outPolygon = Polygon();
            for (int k = 0; k < subjectPolygon.vertices.size(); k++){
                Vector curVertex = subjectPolygon.vertices[k];
                Vector prevVertex = subjectPolygon.vertices[(k > 0)?(k-1):(subjectPolygon.vertices.size()-1)];

                Vector intersection = prevVertex + dot(M - prevVertex, P[i]-P[j])/dot(curVertex - prevVertex, P[i]-P[j])*(curVertex - prevVertex);
                
                if (dot(curVertex - M, P[j] - P[i]) < 0){
                    if (dot(prevVertex - M, P[j] - P[i]) >= 0){
                        outPolygon.vertices.push_back(intersection);
                    }
                    outPolygon.vertices.push_back(curVertex);
                }
                else if (dot(prevVertex - M, P[j] - P[i]) < 0){
                    outPolygon.vertices.push_back(intersection);
                }
            }
            subjectPolygon = outPolygon; 
        }
        VoronoiDiagram.push_back(subjectPolygon);
    }

    delete boundingBox;

    return VoronoiDiagram;
}

/* ------------------------------------------------------------------------------------------------------------------------

To save a set of Polygons as an svg file

------------------------------------------------------------------------------------------------------------------------ */

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j].x * 1000), (1000 - polygons[i].vertices[j].y * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}
 
 
// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
    FILE* f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j].x * 1000), (1000-polygons[i].vertices[j].y * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

/* ------------------------------------------------------------------------------------------------------------------------

Rescales the data so that the vertices are between 0 and 1

------------------------------------------------------------------------------------------------------------------------ */


void rescale(std::vector<Polygon> &data){
    double avx = 0, avy = 0, n = 0;
    for (auto polygon : data){
        for (int i = 0 ; i < polygon.vertices.size(); i++){
            avx += polygon.vertices[i].x;
            avy += polygon.vertices[i].y;
        }
        n += 1; 
    }
    for (int i = 0; i < n; i++){
        for (int k = 0 ; k < data[i].vertices.size(); k++){
            data[i].vertices[k].x *= 1/avx;
            data[i].vertices[k].y *= 1/avy; 
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------------

Main function

------------------------------------------------------------------------------------------------------------------------ */

int main(){
    /*
    std::vector<Vector> data1 = {};
    std::vector<double> weights1 = {};
    for (int i = 2 ; i < 5 ; i++){
        for (int j = 2; j < 5 ; j++){
            data1.push_back(Vector(i, j));
            weights1.push_back(1);
        }
    }

    std::vector<Vector> data2 = {};
    std::vector<double> weights2 = {};
    for (int i = 30 ; i < 67 ; i++){
        for (int j = 30; j < 67 ; j++){
            data2.push_back(Vector(i, j));
            weights2.push_back(1);
        }
    }

    Polygon box1 = *getBoundingBox(data1); 
    Polygon box2 = *getBoundingBox(data2); 

    std::vector<Polygon> test = {box1, box2};
    */


    std::vector<Vector> data = {};
    std::vector<double> weights = {};
    for (int i = 0 ; i < 5 ; i++){
        for (int j = 0; j < 5 ; j++){
            data.push_back(Vector(i, j));
            weights.push_back(1);
        }
    }
    
    std::vector<Polygon> diagram = getVoronoi(data, weights);
    std::cout<< diagram[1].vertices[3].to_string();
    rescale(diagram);
    std::cout<< diagram[1].vertices[3].to_string();

    save_svg(diagram, "diagram.svg");
    return 0;
}