#ifndef CYLINTERP_GEOMETRY_HH
#define CYLINTERP_GEOMETRY_HH

#include "Tools.hh"
#include <vector>
#include <fstream>

class UniformCylindricalGrid {
    public:
        UniformCylindricalGrid(double r0, double r1,
            double z0, double z1,
            int nr, int nz,
            int n_first_ring);
        ~UniformCylindricalGrid();

        void TetraIndices(double points[], int tetra_indices[], int i, int i_out);
        void ExportClass();

    // protected:
        //Setup methods
        void MakePolarCoords();
        void SetGrid();

        double m_dr0;
        double m_dr1;
        double m_dz0;
        double m_dz1;
        int m_inr;
        int m_inz;
        int m_in_first_ring;

        //Radial points
        std::vector<double> m_dr_grid;
        std::vector<double> m_dr_coords;
        std::vector<double> m_dtheta_coords;
        std::vector<int> m_in_theta_in_r_slice;
        std::vector<double> m_dd_theta_in_r_slice;
        std::vector<int> m_icum_n_theta_below_r;
        std::vector<double> m_dunique_r_coords;
        double m_ddr;

        //z points
        std::vector<double> m_dz_coords;
        double m_ddz;

        //Combined grids
        double *m_dpolar_cs_grid = nullptr;
        double *m_dcart_cs_grid = nullptr;

        //Hard coded extra thing
        int m_itetra_iter[8][3] = {{0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1}};
};

#endif