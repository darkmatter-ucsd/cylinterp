#ifndef CYLINTERP_PHYSICS_HH
#define CYLINTERP_PHYSICS_HH

#include "Units.hh"
#include "Interpolator.hh"
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <omp.h>

class RTPC {
    public:
        RTPC(Interpolator* field,
            std::vector<double> &vd_E,
            std::vector<double> &vd_vd,
            double sag,
            int n_cath_wires,
            double r_max_det);
        ~RTPC();

        double NearestCathodeDistance(double polar_points[], double cart_points[], int i);
        double VdInterp(double E_norm);
        double LongDiffusion(double E_norm);
        void Drift(double r[], int n_pts,
            double dr_zmin, double dr_zmax,
            double dr_rmin, double dr_rmax,
            bool tracking, bool diffusion, pcg32& rng,
            double end_time[], double end_pos[], int status[],
            double tracks[]);

        /*
        Key differences from the python:
        - Do the sampling of the points r in main
        - We're never not drifting electrons, so q = -1 always
        - Cathode threshold, dt, and recursion limits are defined as members
        */
        double m_ddt = 10.*ns;
        int m_irecursion_limit = 3000;
        double m_dcath_thresh = 300.*um;
    
    // protected:
        Interpolator *m_field;

        //Drift velocity params
        std::vector<double> m_dvd_E; //E_field
        std::vector<double> m_dvd_vd; //V_drift

        //Diffusion params
        double m_dlxe_trans_diff = 55*cm*cm/s;
        double m_dlongdiff_A = 15.31;
        double m_dlongdiff_B = 53.61;
        double m_dlongdiff_C = 26.21;
        double m_dlongdiff_E0 = 99.72156937;

        //Cathode params
        double m_dsag;
        int m_in_cath_wires;
        double m_dr_max_det;
        double m_dcath_dth;

        std::vector<double> m_dcath_angles;

        bool m_bregular_vd_grid = false;
        double m_ddE;
};

#endif