/* ------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------- 

CSE 306 - Computer Graphics : Louis de Benoist
Project #2

------------------------------------------------------------------------------------------------------------------------ 
------------------------------------------------------------------------------------------------------------------------ */

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb_image_write.h"
#include <bits/stdc++.h>
#include "lbfgs.c"
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

        Vector(double x_ = 0. , double y_ = 0.){ 
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
double cross(const Vector& v, const Vector& w){
    return v.x*w.y - v.y*w.x; 
}
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

    //Vector v1 = Vector(smallest_x, smallest_y), v2 = Vector(smallest_x, biggest_y), v3 = Vector(biggest_x, smallest_y), v4 = Vector(biggest_x, biggest_y); 
    Vector v1 = {0,0}, v2 = {0,1}, v3 = {1,0}, v4 = {1,1}; 
    std::pair<Vector, Vector> e1 = {v1, v3}, e2 = {v2, v4}, e3 = {v1, v2}, e4 = {v3, v4};

    return new Polygon({v1, v3, v4, v2}, {e1, e2, e3, e4});
}

/* ------------------------------------------------------------------------------------------------------------------------

Voronoi Parallel Linear Enumeration (Power Diagram)
 - takes a bunch of data points (Vectors) and returns a vector of Polygons such that the ith Polygon corresponds 
   to the Voronoi polygon of the ith data point
- Also takes into account weights for each cell (Power Diagram)
------------------------------------------------------------------------------------------------------------------------ */

std::vector<Polygon> getPowerDiagram(std::vector<Vector> P, std::vector<double> weights){
    std::vector<Polygon> VoronoiDiagram;
    Polygon* boundingBox = getBoundingBox(P); 

    // get the voronoi Polygon for each data point
    //#pragma omp parallel for schedule(dynamic)
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
    double minx = INF, maxx = -INF; 
    double miny = INF, maxy = -INF;
    int n = 0;  

    for (auto polygon : data){
        for (int i = 0 ; i < polygon.vertices.size(); i++){
            if (polygon.vertices[i].x < minx) minx = polygon.vertices[i].x;
            if (polygon.vertices[i].y < miny) miny = polygon.vertices[i].y;
            if (polygon.vertices[i].x > maxx) maxx = polygon.vertices[i].x;
            if (polygon.vertices[i].y > maxy) maxy = polygon.vertices[i].y;
        }
        n += 1; 
    }
    for (int i = 0; i < n; i++){
        for (int k = 0 ; k < data[i].vertices.size(); k++){
            data[i].vertices[k].x = (data[i].vertices[k].x - minx) / (maxx - minx) ;
            data[i].vertices[k].y = (data[i].vertices[k].y - miny) / (maxy - miny) ;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------------

Get the area of a polygon

------------------------------------------------------------------------------------------------------------------------ */

double getArea(Polygon P){
    double A = 0; 
    for (int i = 0; i < P.vertices.size(); i++){
        A += cross(P.vertices[i%P.vertices.size()], P.vertices[(i+1)%P.vertices.size()]);
    }
    return A / 2;
}

/* ------------------------------------------------------------------------------------------------------------------------

Semi-Optimal Transport Implementation -- influenced by the example detailed in the documentation website

------------------------------------------------------------------------------------------------------------------------ */

static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls){
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

struct powerArguments {
    const std::vector<Vector> data;
    const std::vector<double> lambdas;
};

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *nabla_g, const int n, const lbfgsfloatval_t step){
    auto [data, lambdas] = *(powerArguments *) instance;
    std::vector<double> w; 
    for (int i = 0; i < data.size(); i++) w.push_back(x[i]);

    std::vector<Polygon> diagram = getPowerDiagram(data, w);
    lbfgsfloatval_t g = 0.0;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        // we compute g analytically
        auto& P = diagram[i].vertices;
        int N = P.size(); 
        double I = 0;
        for (int k = 1; k <= N; k++){
           // Matthias Hasler gave me the idea to reformulate the original equation in terms of dot and cross products
           Vector a = P[k-1], b = P[k%N];
           I += cross(a,b) * (dot(a,a) + dot(a,b) + dot(b,b) - 4*dot(data[i],a+b) + 6*dot(data[i],data[i])) / 12;
        }
        double A = getArea(diagram[i]); 
        #pragma omp atomic
        g += I - A*w[i] + lambdas[i]*w[i];
        // we compute nabla g analytically
        nabla_g[i] = A - lambdas[i];
    }
    return -g;
}

std::vector<double> getOptimalWeights(std::vector<Vector> data, std::vector<double> lambdas){
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(data.size());

    // initialize x
    std::fill(x, x + data.size(), 0); 

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    powerArguments args = {data, lambdas};

    int ret = lbfgs(data.size(), x, &fx, evaluate, progress, &args, &param);
    std::cout << "RET:" << ret << "\n";
    std::vector<double> weights; 
    for (int i = 0; i < data.size(); i++) weights.push_back(x[i]);
    return weights;
}

/* ------------------------------------------------------------------------------------------------------------------------

 Gets the centroid of a polygon (influenced by Matthias Hasler)

------------------------------------------------------------------------------------------------------------------------ */
Vector getCentroid(Polygon P, Vector point){
    if (P.vertices.size() <= 0){
        return Vector(std::clamp(point.x,0.,1.), std::clamp(point.y,0.,1.)); 
    }

    Vector C = Vector(0,0);
    for (int i = 0; i < P.vertices.size(); i++){
        Vector &a = P.vertices[i];
        Vector &b = P.vertices[(i+1)%P.vertices.size()];
        C += (a+b)* cross(a,b);
    }

    return (1 / (6 * getArea(P))) * C;
}

/* ------------------------------------------------------------------------------------------------------------------------

 Gallouet Merigot scheme for fluid simulation

------------------------------------------------------------------------------------------------------------------------ */
void save_frame(const std::vector<Polygon> &cells, std::string filename, int N, int frameid) ;

void fluidSimulation(int N, int num_frames){
    std::vector<Vector> points;
    std::vector<Vector> velocities(2*N);
    std::vector<double> lambdas;

    std::mt19937 rng(time(NULL)); // for generating random numbers
    std::uniform_real_distribution dis(0.,1.);

    // N is the number of water molecules

    for (int i = 0 ; i < 2*N ; i++){
        double r1 = dis(rng), r2 = dis(rng);
        // if it is water
        if (i < N){
            auto v = Vector(r1/4 + 0.375, r2/4 + 0.375);
            lambdas.push_back(0.3/N);
            points.push_back(v);
        }
        // if it is air
        else {
            auto v = Vector(r1, r2); 
            lambdas.push_back(0.7/N);
            points.push_back(v);
        }
    }

    double m = 200;
    Vector Fg = m*Vector(0, -9.8);
    double dt = 0.01; 
    double eps = 0.004;
    
    for (int i = 0 ; i < num_frames; i++){
        auto optimal_weights = getOptimalWeights(points, lambdas);
        auto power_diagram = getPowerDiagram(points, optimal_weights);

        save_frame(power_diagram, "frame", N, i);

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < 2*N; j++){
            Vector F = pow(eps, -2) * (getCentroid(power_diagram[j], points[j]) - points[j]);

            // if it's water, we add gravitational force
            if (j < N) F = F + Fg; 

            velocities[j] += (1/m)*F * dt; 
            points[j] += velocities[j] * dt;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------------

Save frame (Helped by Matthias Hasler)

------------------------------------------------------------------------------------------------------------------------ */

void save_frame(const std::vector<Polygon> &cells, std::string filename, int N, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {
        int n = cells[i].vertices.size();
        if (n==0) continue;
        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j].x);
            bminy = std::min(bminy, cells[i].vertices[j].y);
            bmaxx = std::max(bmaxx, cells[i].vertices[j].x);
            bmaxy = std::max(bmaxy, cells[i].vertices[j].y);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::min(0., W * bmaxx));
        bmaxy = std::max(H-1., std::min(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                Polygon cell = cells[i];
                Vector v = (1./W) * Vector(x,y);
                bool in = 1;
                double dst = 1e9;
                for (int j = 0; j < n; j++) {
                    Vector a = cell.vertices[j], b = cell.vertices[(j+1) % n];
                    b = b-a;
                    a = v-a;

                    double alt = cross(b,a) / norm(b);
                    if (alt < -1e-3) { in = 0; break; }
                    dst = std::min(dst, alt);
                }
                if (!in) continue;

                int idx = ((H - y - 1)*W + x) * 3;
                if (dst <= 2e-3) { // border
                    image[idx + 0] = 0;
                    image[idx + 1] = 0;
                    image[idx + 2] = 0;
                } else if (i < N) { // is water
                    image[idx + 0] = 0;
                    image[idx + 1] = 0;
                    image[idx + 2] = 255;
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

/* ------------------------------------------------------------------------------------------------------------------------

Main function

------------------------------------------------------------------------------------------------------------------------ */

int main(){
    auto start = high_resolution_clock::now(); // to measure execution time
    fluidSimulation(200, 40);
    auto stop = high_resolution_clock::now(); auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Rendering time: " << duration.count()*pow(10, -6) << " seconds\n"; 

   /*
   std::mt19937 rng(time(NULL)); // for generating random numbers
    std::uniform_real_distribution dis(0.,1.);

    std::vector<Vector> data;
    std::vector<double> lambdas;
    for (int i = 0 ; i < 1000 ; i++){
        double r1 = dis(rng);
        double r2 = dis(rng);
        auto v = Vector(r1, r2); 
        data.push_back(v);
        lambdas.push_back(exp(-5 * norm(v - Vector(0.5, 0.5))));
    }
    {
        double s = std::accumulate(lambdas.begin(), lambdas.end(), 0.);
        std::for_each(lambdas.begin(), lambdas.end(), [=](double& x){ x /= s; });
    }

    auto start = high_resolution_clock::now();
    std::vector<double> optimal_weights = getOptimalWeights(data, lambdas);
    std::vector<Polygon> power_diagram = getPowerDiagram(data, optimal_weights);
    auto stop = high_resolution_clock::now(); auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Rendering time: " << duration.count()*pow(10, -6) << " seconds\n"; 
    //rescale(power_diagram);
    save_svg(power_diagram, "power_diagram_with_semiOptimalTransport.svg");
    */


    /*
    // Power diagram with weights 1 (i.e. a Voronoi diagram )
    std::vector<Vector> data = {Vector(0,0), Vector(2, 1), Vector(0,2)};
    std::vector<double> weights = {1,1,1};
    auto start = high_resolution_clock::now(); // to measure execution time
    std::vector<Polygon> voronoi_diagram = getPowerDiagram(data, weights);
    auto stop = high_resolution_clock::now(); auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Rendering time: " << duration.count()*pow(10, -6) << " seconds\n"; 
    rescale(voronoi_diagram);
    save_svg(voronoi_diagram, "power_diagram.svg");

    // Semi-Optimal Transport
    std::vector<Vector> data2= {Vector(0,0), Vector(2, 1), Vector(0,2)};
    std::vector<double> lambdas(3);
    for (int i = 0 ; i < data2.size(); i++){
        auto v = data2[i];
        auto C = Vector(1,1);
        lambdas[i] = exp(-norm(v - C)/0.02);
    }


    lambdas = {0, 0, 9};

    start = high_resolution_clock::now();
    std::vector<double> optimal_weights = getOptimalWeights(data2, lambdas);
    std::vector<Polygon> power_diagram = getPowerDiagram(data2, optimal_weights);
    stop = high_resolution_clock::now(); duration = duration_cast<microseconds>(stop - start);
    std::cout << "Rendering time: " << duration.count()*pow(10, -6) << " seconds\n"; 
    rescale(power_diagram);
    save_svg(power_diagram, "power_diagram_with_semiOptimalTransport.svg");
    */
    
    return 0;
}