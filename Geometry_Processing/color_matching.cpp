/* ------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------- 

CSE 306 - Computer Graphics : Louis de Benoist
Color matching an image based on a template

------------------------------------------------------------------------------------------------------------------------ 
------------------------------------------------------------------------------------------------------------------------ */

#define STB_IMAGE_IMPLEMENTATION
#include "main.cpp"
#include "stb_image.h"

/* ------------------------------------------------------------------------------------------------------------------------

Returns random vector

------------------------------------------------------------------------------------------------------------------------ */

Vector random_direction(){
    double r1 = ((double) rand() / (RAND_MAX)), r2 = ((double) rand() / (RAND_MAX));
    double x = cos(2*M_PI*r1)*sqrt(r2*(1-r2));
    double y = sin(2*M_PI*r1)*sqrt(r2*(1-r2));
    double z = 1-2*r2;
    return Vector(x,y,z);
}

/* ------------------------------------------------------------------------------------------------------------------------

Sliced optimal transport color transfer algorithm

------------------------------------------------------------------------------------------------------------------------ */

void sliced_optimal_transport(unsigned char I[], unsigned char M[], int H, int W, int num_iters){
    for (int iter = 0 ; iter < num_iters ; iter++){
        // get random vector on which we will project the pixels
        Vector v = random_direction(); 

        std::vector<std::pair<double, std::pair<int, int>>> projI, projM;
        
        // for each pixel, we store the dot product and indices (i,j) of the pixel
        for (int i = 0 ; i < W ; i++){
            for (int j = 0 ; j < H ; j++){
                Vector pixel_I = Vector(I[j*W*3 + i*3 + 0], I[j*W*3 + i*3 + 1], I[j*W*3 + i*3 + 2]);
                Vector pixel_M= Vector(M[j*W*3 + i*3 + 0], M[j*W*3 + i*3 + 1], M[j*W*3 + i*3 + 2]);
                projI.push_back(std::make_pair(dot(pixel_I, v),std::make_pair(i,j)));
                projM.push_back(std::make_pair(dot(pixel_M, v),std::make_pair(i,j)));
            }
        }

        // sort according to the dot product
        std::sort(projI.begin(), projI.end()); 
        std::sort(projM.begin(), projM.end()); 

        // modify the original point cloud
        for (int i = 0 ; i < W ; i++){
            for (int j = 0 ; j < H ; j++){
                if (j == 175){
                    return;
                }
                // for r value
                int i_coord = std::get<0>(std::get<1>(projI[j*W*3 + i*3 + 0])); 
                int j_coord = std::get<1>(std::get<1>(projI[j*W*3 + i*3 + 0])); 

                I[j_coord*W*3 + i_coord*3 + 0] += (std::get<0>(projM[j*W*3 + i*3 + 0]) - std::get<0>(projI[j*W*3 + i*3 + 0]))*v.x;

                // for g value
                i_coord = std::get<0>(std::get<1>(projI[j*W*3 + i*3 + 1])); 
                j_coord = std::get<1>(std::get<1>(projI[j*W*3 + i*3 + 1])); 

                I[j_coord*W*3 + i_coord*3 + 1] += (std::get<0>(projM[j*W*3 + i*3 + 1]) - std::get<0>(projI[j*W*3 + i*3 + 1]))*v.y;

                // for b value
                i_coord = std::get<0>(std::get<1>(projI[j*W*3 + i*3 + 2])); 
                j_coord = std::get<1>(std::get<1>(projI[j*W*3 + i*3 + 2])); 

                I[j_coord*W*3 + i_coord*3 + 2] += (std::get<0>(projM[j*W*3 + i*3 + 2]) - std::get<0>(projI[j*W*3 + i*3 + 2]))*v.z;
            }
        }
    }
    stbi_write_jpg("Outputs/color_matched.jpg", W, H, 3, I, W * sizeof(int));
}

int main(){
    srand(time(NULL)); // for generating random numbers
    int W, H, bpp;
    int num_iters = 100; 

    unsigned char* I = stbi_load("Outputs/shrekdonkey.jpg", &W, &H, &bpp, 3);
    unsigned char* M = stbi_load("Other/template.jpg", &W, &H, &bpp, 3);

    sliced_optimal_transport(I, M, H, W, num_iters);
    return 0;
}