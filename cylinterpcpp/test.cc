#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include "pcg_random.hh"
#include "Tools.hh"
#include "Geometry.hh"

int main()
{
    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Make a random number engine
    pcg32 rng(seed_source);

    int n_points = 10;
    double points[3*n_points];
    double points2[3*n_points];
    double points3[3*n_points];

    double cartesian_points[3*n_points];
    double cartesian_points2[3*n_points];
    double cartesian_points3[3*n_points];

    GenerateRandPointsPipe(n_points, 0., 2.5, -1, 11, points, rng);
    GenerateRandPointsPipe(n_points, 0., 2.5, -1, 11, points2, rng);
    GenerateRandPointsPipe(n_points, 0., 2.5, -1, 11, points3, rng);

    for (int i=0; i<n_points; i++){
        ToCartesian(points, cartesian_points, i);
    }
    for (int i=0; i<n_points; i++){
        ToCartesian(points2, cartesian_points2, i);
    }
    for (int i=0; i<n_points; i++){
        ToCartesian(points3, cartesian_points3, i);
    }
    
    std::cout << "Randomly generated polar coordinates: \n";
    PrintArr2d(points, 3, n_points);
    std::cout << "\n";
    PrintArr2d(points2, 3, n_points);
    std::cout << "\n";
    PrintArr2d(points3, 3, n_points);
    std::cout << "\n";

    std::cout << "\nRandomly generated cartesian coordinates: \n";
    PrintArr2d(cartesian_points, 3, n_points);
    std::cout << "\n";
    PrintArr2d(cartesian_points2, 3, n_points);
    std::cout << "\n";
    PrintArr2d(cartesian_points3, 3, n_points);
    std::cout << "\n";


    double axb[3*n_points];
    for (int i=0; i<n_points; i++){
        CrossProduct(cartesian_points, cartesian_points2, axb, i, i);
    }
    std::cout<<"a x b: \n";
    PrintArr2d(axb, 3, n_points);

    std::cout<< "\nScalar triple products: \n[";
    for (int i=0; i<n_points; i++){
        double stp = ScalarTripleProduct(cartesian_points, cartesian_points2, cartesian_points3, i);
        std::cout << stp<< ",";
    }
    std::cout<<"]\n";

}