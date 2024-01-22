/*
Simulate 10k electrons for a bunch of random positions within SanDiX with diffusion.
There is an 8-fold symmetry, so we only need to simulate 1/8th of the detector
The goal is to get:
 - S2 electron arrival time distributions for waveform simulation
 - Charge insensitive volume information
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <omp.h>

#include "pcg_random.hh"
#include "Tools.hh"
#include "Geometry.hh"
#include "Interpolator.hh"
#include "Physics.hh"
#include "Timer.hh"

int main(int argc, char* argv[])
{
    int c=0;
    int n_threads = 1;
    std::string fieldfile = "4500V_1mmSag_Run29_Field.bin";
    std::string vd_E_file = "vd_E.bin";
    std::string vd_vd_file = "vd_vd.bin";
    std::string outfile = "./";
    int n_points = 10000;
    int n_elec_per_s2 = 1000;
    double r_min_gen = 0.03439406779661017*cm;
    double r_max_gen = 2.5*cm;
    double z_min_gen = 5.;
    double z_max_gen = 11.;
    double theta_min_gen = 0.;
    double theta_max_gen = 2*M_PI/8;
    int chunk_size = 1000; //Number of electrons to produce in a batch

    while ((c = getopt(argc, argv, "t:E:v:V:o:n:e:r:R:z:Z:u:U:M:")) != -1)
	{
		switch (c)
		{
		case 't':
			n_threads = atoi(optarg);;
			break;

		case 'E':
			fieldfile.assign(optarg);
			break;

		case 'v':
			vd_E_file.assign(optarg);
			break;

		case 'V':
			vd_vd_file.assign(optarg);
			break;

		case 'o':
			outfile.assign(optarg);
			break;
        
        case 'n':
            n_points = atoi(optarg);
            break;

        case 'e':
            n_elec_per_s2 = atoi(optarg);
            break;
        
        case 'r':
            r_min_gen = atof(optarg);
            break;
        
        case 'R':
            r_max_gen = atof(optarg);
            break;

        case 'z':
            z_min_gen = atof(optarg);
            break;

        case 'Z':
            z_max_gen = atof(optarg);
            break;

        case 'u':
            theta_min_gen = atof(optarg);
            break;

        case 'U':
            theta_max_gen = atof(optarg);
            break;

        case 'M':
            chunk_size = atoi(optarg);
            break;

		default:
			break;
		}
	}


    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Make a random number engine
    pcg32 rng(seed_source);
    
    // Get the field data
    // NOTE: The field file is kind of big and will exceed github's maximum,
    // Contact Jianyang (jiq019@ucsd.edu) if you want it!
    size_t fieldfile_size = FileSize(fieldfile);
    double* field = new double[fieldfile_size/sizeof(double)];
    LoadArrayFile(fieldfile, field, fieldfile_size);
    std::vector<double> vd_E;
    std::vector<double> vd_vd;
    LoadArrayFile(vd_vd_file, vd_vd);
    LoadArrayFile(vd_E_file, vd_E);

    //Make the interpolator and RTPC classes
    Interpolator *E_field = new Interpolator(5E-4, 4., -1., 11., 60, 100, 8, field, 3);
    RTPC *rtpc = new RTPC(E_field, vd_E, vd_vd, 0.1*cm, 20, 2.5*cm);
    rtpc->m_inum_threads = n_threads;

    int n_chunks = n_points/chunk_size;
    for (int c=0; c<n_chunks; c++){
        double *points_unique = new double[3*chunk_size];
        GenerateRandPointsPipe(chunk_size, r_min_gen, r_max_gen, z_min_gen, z_max_gen, theta_min_gen, theta_max_gen, points_unique, rng);

        double *points = new double[3*chunk_size*n_elec_per_s2];
        //Copy the points unique over
        int s2_i;
        int r_i;
        int d;
        #pragma omp parallel for shared(points, points_unique) private(s2_i, r_i, d) num_threads(n_threads) collapse(3)
        for (s2_i=0; s2_i<chunk_size; s2_i++){
            for (r_i=0; r_i<n_elec_per_s2; r_i++){
                for (d=0; d<3; d++){
                    points[3*n_elec_per_s2*s2_i+3*r_i+d] = points_unique[3*s2_i + d];
                }
            }
        }
    

        double tracks[3] = {0.};
        // int tracksize = 3*chunk_size*n_elec_per_s2*rtpc->m_irecursion_limit;
        // double *tracks = new double[tracksize];
        // std::fill(tracks, tracks + tracksize, 0.);
        double *end_time = new double[chunk_size*n_elec_per_s2];
        double *end_pos = new double[3*chunk_size*n_elec_per_s2];
        int *status = new int[chunk_size*n_elec_per_s2];
        std::cout << "Drifting chunk "<< c <<" took: \n";
        {
            Timer timer;
            rtpc->Drift(points, chunk_size*n_elec_per_s2,
                -1., 11., rtpc->m_field->m_dunique_r_coords[0], 2.5,
                false, true, rng,
                end_time, end_pos, status, tracks);
        }

        // std::string tracks_file = outfile + "tracks_" + std::to_string(c) + ".bin";
        // std::ofstream TracksOut(tracks_file, std::ios::out | std::ios::binary);
        // TracksOut.write((char*)&tracks[0], sizeof(double) * 3*chunk_size*n_elec_per_s2*rtpc->m_irecursion_limit);
        // TracksOut.close();
        // delete[] tracks;

        std::string init_pos_file = outfile + "init_pos_" + std::to_string(c) + ".bin";
        std::ofstream InitPos(init_pos_file, std::ios::out | std::ios::binary);
        InitPos.write((char*)&points_unique[0], sizeof(double) * 3*chunk_size);
        InitPos.close();
        delete[] points_unique;

        std::string end_time_file = outfile + "end_time_" + std::to_string(c) + ".bin";
        std::ofstream EndTimeOut(end_time_file, std::ios::out | std::ios::binary);
        EndTimeOut.write((char*)&end_time[0], sizeof(double) * chunk_size*n_elec_per_s2);
        EndTimeOut.close();
        delete[] end_time;

        std::string end_pos_file = outfile + "end_pos_" + std::to_string(c) + ".bin";
        std::ofstream EndPosOut(end_pos_file, std::ios::out | std::ios::binary);
        EndPosOut.write((char*)&end_pos[0], sizeof(double) * 3*chunk_size*n_elec_per_s2);
        EndPosOut.close();
        delete[] end_pos;

        std::string status_file = outfile + "status_" + std::to_string(c) + ".bin";
        std::ofstream StatusOut(status_file, std::ios::out | std::ios::binary);
        StatusOut.write((char*)&status[0], sizeof(int) * chunk_size*n_elec_per_s2);
        StatusOut.close();
        delete[] status;

        //Delete points for the next chunk
        delete[] points;
    }

    // delete tracks;
    delete E_field; //Deleting the interpolator automatically deletes the field
    delete rtpc;
}