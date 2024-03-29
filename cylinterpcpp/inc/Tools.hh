#ifndef CYLINTERP_TOOLS_HH
#define CYLINTERP_TOOLS_HH

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include "pcg_random.hh"
#include <omp.h>

template<typename T>
static void PrintArr2d(T tArray[], int iNx, int iNy){
    std::cout << "[";
    for (int i = 0; i < iNy; i++){
        std::cout << "[";
        for (int j = 0; j<iNx; j++){
            std::cout << tArray[iNx*i + j];
            if (j!=iNx-1)
                std::cout << ",";
        }
        std::cout << "]";
        if (i!=iNy-1)
            std::cout<<",\n";
    }
    std::cout<<"]\n";
}

static size_t FileSize(std::string& filename){
    std::ifstream input_data(filename.c_str(), std::ios::binary);
    input_data.seekg(0, std::ios_base::end);
    size_t size = input_data.tellg();
    return size;
}

template<typename T>
static void LoadArrayFile(std::string& filename, T arr[], size_t size){
    std::ifstream input_data(filename.c_str(), std::ios::binary);
    input_data.seekg(0, std::ios_base::beg);

    input_data.read(reinterpret_cast<char*>(arr), size);
    input_data.close();
}

template<typename T>
static void LoadArrayFile(std::string& filename, std::vector<T> &arr){
    std::ifstream input_data(filename.c_str());
    input_data.seekg(0, std::ios_base::end);
    std::size_t size = input_data.tellg();
    input_data.seekg(0, std::ios_base::beg);

    arr.resize(size / sizeof(T), 0);
    input_data.read((char*)(&arr[0]), size);
    input_data.close();
}

static std::vector<double> Linspace(double a0, double a1, int n){
    std::vector<double> a_out;
    for (int i=0; i<n; i++){
        a_out.emplace_back(a0 + i*(a1-a0)/(n-1));
    }

    return a_out;
}

static void ToCartesian(double polar_points[],
    double cartesian_points[], int i, int i_out){
    
    //polar_points: 1-d array of polar points which is actually hidden as a 2-d array of 3-d points.
    //  This "3" is hard coded, you're not doing anything in anything other than 3 dimensions anyway
    //cartesian_points: the 1-d array of cartesian points corresponding to the polar points
    //i: The point being looked at

    //0th index of cartesian is x, 1st index of polar is r
    cartesian_points[3*i_out] = polar_points[3*i+1] * std::cos(polar_points[3*i+2]);
    cartesian_points[3*i_out+1] = polar_points[3*i+1] * std::sin(polar_points[3*i+2]);
    cartesian_points[3*i_out+2] = polar_points[3*i];
}

static void ToPolar(double polar_points[],
    double cartesian_points[], int i, int i_out){
    
    //polar_points: 1-d array of polar points which is actually hidden as a 2-d array of 3-d points.
    //  This "3" is hard coded, you're not doing anything in anything other than 3 dimensions anyway
    //cartesian_points: the 1-d array of cartesian points corresponding to the polar points
    //i: The point being looked at

    double angle = std::atan2(cartesian_points[3*i+1], cartesian_points[3*i]);
    if (angle < 0) angle += 2*M_PI;

    polar_points[3*i_out] = cartesian_points[3*i+2];
    polar_points[3*i_out+1] = std::sqrt(cartesian_points[3*i]*cartesian_points[3*i]+
        cartesian_points[3*i+1]*cartesian_points[3*i+1]);
    polar_points[3*i_out+2] = angle;
}

static void CrossProduct(double a[], double b[], double a_cross_b[],
    int i_in, int i_out){
    // a: an array with at least 3 elements
    // b: ''
    // i_in: the 3-vector in question to be extracted from a and b
    // a_cross_b: axb
    // i_out: the 3-vector in question that lies on the output of a_cross_b

    a_cross_b[3*i_out] = a[3*i_in+1]*b[3*i_in+2]-a[3*i_in+2]*b[3*i_in+1];
    a_cross_b[3*i_out+1] = a[3*i_in+2]*b[3*i_in]-a[3*i_in]*b[3*i_in+2];
    a_cross_b[3*i_out+2] = a[3*i_in]*b[3*i_in+1]-a[3*i_in+1]*b[3*i_in];
}

static double ScalarTripleProduct(double a[], double b[], double c[], int i){
    double axb[3];

    CrossProduct(a, b, axb, i, 0);

    double abc = axb[0]*c[3*i]+axb[1]*c[3*i+1]+axb[2]*c[3*i+2];

    return abc;
}

static void GenerateRandPointsPipe(int n_points,
    double r_min, double r_max,
    double z_min, double z_max,
    double th_min, double th_max,
    double points[], pcg32& rng){
    
    std::uniform_real_distribution<double> uni_dist_z(z_min, z_max);
    std::uniform_real_distribution<double> uni_dist_r(r_min*r_min, r_max*r_max);
    std::uniform_real_distribution<double> uni_dist_theta(th_min, th_max);

    #pragma omp parallel for num_threads(4)
    for (int i=0; i < n_points; i++){
        points[3*i] = uni_dist_z(rng);
        points[3*i+1] = std::sqrt(uni_dist_r(rng));
        points[3*i+2] = uni_dist_theta(rng);
    }
}

#endif