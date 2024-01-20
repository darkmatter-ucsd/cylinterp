#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include <fstream>
#include <omp.h>

#include "pcg_random.hh"
#include "Tools.hh"
#include "Geometry.hh"
#include "Interpolator.hh"
#include "Physics.hh"
#include "Timer.hh"

int main()
{
    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Make a random number engine
    pcg32 rng(seed_source);

    std::cout <<"Num threads: "<< omp_get_max_threads()<<"\n";

    int n_points = 100;
    // double fixed_rand_points[3*n_points];

    // GenerateRandPointsPipe(n_points, 0., 2.5, -1, 11, fixed_rand_points, rng);


    double fixed_rand_points[3*n_points] = {0.91738,1.96182,1.48023,
        4.73132,1.34186,1.42465,
        9.67598,1.71061,3.84124,
        4.60419,2.2572,2.93833,
        3.42724,1.89636,3.41235,
        4.64924,1.75566,2.1089,
        8.66973,1.8854,2.83619,
        3.70464,0.830971,1.4201,
        5.48373,1.55557,3.87441,
        0.147654,0.303774,3.44808,
        1.37607,1.54288,3.43999,
        8.16804,1.70794,0.614368,
        6.86776,1.57014,3.55508,
        4.58794,1.31836,0.00244411,
        7.02584,1.85233,1.68586,
        6.50964,2.31881,3.0541,
        -0.568133,0.516047,1.70015,
        10.4431,2.0342,1.03111,
        2.6066,0.773119,3.13455,
        1.62602,2.02589,0.702002,
        1.34281,0.567153,4.00771,
        6.11748,0.506676,2.29367,
        8.37495,2.33224,1.68078,
        0.569711,2.23358,0.911663,
        10.0738,2.42018,1.83573,
        6.79702,1.61558,0.813459,
        7.28601,1.22748,2.93556,
        0.541302,2.2572,2.35181,
        6.45421,1.96408,5.84607,
        3.81909,1.86442,3.0078,
        9.02254,1.02257,2.08277,
        7.68493,2.36589,3.22923,
        3.98321,1.91034,1.65569,
        7.95458,1.66197,5.27991,
        8.37182,2.26756,3.27332,
        -0.183441,1.78167,5.20063,
        4.29273,0.455087,5.42704,
        7.35216,1.61583,4.03161,
        1.86336,2.04108,0.611553,
        6.21618,1.72,4.36145,
        9.41344,1.21072,2.37982,
        0.44772,1.16608,3.78846,
        0.307001,1.92489,4.46445,
        10.3757,1.43508,1.438,
        3.69252,1.66611,1.28029,
        0.247658,0.379611,5.19705,
        6.37379,2.38721,3.02738,
        3.34126,1.30414,1.84264,
        5.4515,1.55237,1.30855,
        5.43054,2.4259,5.90639,
        1.84525,1.05208,3.24216,
        4.02466,1.87538,4.16545,
        0.241712,1.80176,4.65268,
        8.51467,2.23657,5.04718,
        2.136,1.48767,0.576337,
        3.77137,2.0877,3.32904,
        1.22555,2.26071,4.35569,
        3.43676,1.93393,1.89517,
        7.62383,1.1309,0.612481,
        1.10662,2.47674,4.19127,
        9.39617,1.25062,6.27891,
        9.07329,1.25975,0.175064,
        10.6264,1.21311,2.09354,
        -0.54946,2.26438,4.20871,
        0.385467,0.232229,1.29286,
        1.41057,1.18742,4.68084,
        -0.744761,2.21812,0.0061976,
        10.9738,0.964124,2.31921,
        1.30796,0.994178,1.6718,
        9.18013,1.55393,1.32054,
        4.20236,1.73994,2.78742,
        1.06572,0.192482,1.42858,
        3.94253,1.29344,4.1197,
        6.10133,1.90587,2.90578,
        1.47779,1.25761,1.12234,
        0.992317,2.37758,0.425441,
        -0.439995,1.5531,5.98091,
        9.85029,2.23534,0.400166,
        0.803315,0.920057,0.480759,
        -0.000144681,2.01357,0.656627,
        3.97599,1.05916,1.22822,
        8.93875,2.15883,3.98676,
        8.90582,0.787646,5.77091,
        0.299833,0.927657,1.36629,
        6.50496,2.02011,4.24571,
        4.47529,1.76245,1.13156,
        6.80437,2.2365,1.02469,
        5.96441,2.02272,2.31146,
        7.79264,1.3089,0.330132,
        3.88085,2.00258,5.07904,
        7.83826,0.761,2.26563,
        9.71582,1.97275,2.77507,
        6.8991,0.0938258,6.25366,
        8.01934,1.53503,1.86274,
        9.52879,2.14658,4.22747,
        1.47811,1.52122,2.88096,
        8.01448,1.48852,3.51634,
        0.889599,1.49401,3.06632,
        -0.0288105,2.29321,1.52589,
        3.3734,0.984517,4.55499};

    double cart_points[3*n_points];
    for (int i=0; i<n_points; i++) ToCartesian(fixed_rand_points, cart_points, i, i);

    UniformCylindricalGrid *unicylgrid = new UniformCylindricalGrid(5E-4, 4., -1., 11., 60, 100, 8);
    
    // int tetra_indices[4*n_points];
    // for (int i=0; i<n_points; i++){
    //     unicylgrid->TetraIndices(fixed_rand_points, tetra_indices, i, i);
    //     // std::cout << "Worked for index "<<i<<"\n";
    // }

    // std::cout << "Tetrahedral indices for the points: \n";
    // PrintArr2d(tetra_indices, 4, n_points);
    
    // //Get the field data
    std::string fieldfile = "4500V_1mmSag_Run29_Field.bin";
    size_t fieldfile_size = FileSize(fieldfile);
    double* field = new double[fieldfile_size/sizeof(double)];
    LoadArrayFile(fieldfile, field, fieldfile_size);
    std::vector<double> vd_E;
    std::vector<double> vd_vd;
    std::string vd_vd_file = "vd_vd.bin";
    std::string vd_E_file = "vd_E.bin";
    LoadArrayFile(vd_vd_file, vd_vd);
    LoadArrayFile(vd_E_file, vd_E);

    // std::cout << "\n";
    // for (int i=0; i<10; i++){
    //     std::cout << vd_E[i] << " ";
    // }
    // std::cout << "\n";
    // for (int i=0; i<10; i++){
    //     std::cout << vd_vd[i] << " ";
    // }

    Interpolator *interpol = new Interpolator(5E-4, 4., -1., 11., 60, 100, 8, field, 3);

    // double interp_field[3*n_points];
    // for (int i=0; i<n_points; i++){
    //     interpol->Interpolate(fixed_rand_points, cart_points, interp_field, i, i);
    // }

    // std::cout << "\nInterpolated fields: \n";
    // PrintArr2d(interp_field, 3, n_points);

    RTPC *rtpc = new RTPC(interpol, vd_E, vd_vd, 0.1*cm, 20, 2.5*cm);
    // std::cout << "\nInterpolated v_drift: \n[";
    // for (int i=0; i<n_points; i++){
    //     double E_norm = std::sqrt(interp_field[3*i]*interp_field[3*i]+
    //         interp_field[3*i+1]*interp_field[3*i+1]+
    //         interp_field[3*i+2]*interp_field[3*i+2]);
    //     std::cout << rtpc->VdInterp(E_norm);
    //     if (i!=n_points-2) std::cout <<", ";
    // }
    // std::cout << "]\n";

    // std::cout << "\nNearest cathode distances: \n[";
    // for (int i=0; i<n_points; i++){
    //     std::cout << rtpc->NearestCathodeDistance(fixed_rand_points, cart_points, i);
    //     if (i!=n_points-2) std::cout <<", ";
    // }
    // std::cout << "]\n";

    // double tracks = {0.};
    // double tracks[3] = {0.};
    int tracksize = 3*n_points*rtpc->m_irecursion_limit;
    double *tracks = new double[tracksize];
    std::fill(tracks, tracks + tracksize, 0.);
    double end_time[n_points];
    double end_pos[3*n_points];
    int status[n_points];

    std::cout << "Drifting took: \n";
    {
        Timer timer;
        rtpc->Drift(fixed_rand_points, n_points,
            -1., 11., rtpc->m_field->m_dunique_r_coords[0], 2.5,
            true, true, rng,
            end_time, end_pos, status, tracks);
    }
    //Output the tracks
    std::ofstream TracksOut("cpp_cylinterp_tracks.bin", std::ios::out | std::ios::binary);
    TracksOut.write((char*)&tracks[0], sizeof(double) * 3*n_points*rtpc->m_irecursion_limit);
    TracksOut.close();

    std::ofstream EndPosOut("cpp_cylinterp_endpos.bin", std::ios::out | std::ios::binary);
    EndPosOut.write((char*)&end_pos[0], sizeof(double) * 3*n_points);
    EndPosOut.close();

    // delete tracks;
    delete interpol; //Deleting the interpolator automatically deletes the field
    delete rtpc;
    delete unicylgrid;
}