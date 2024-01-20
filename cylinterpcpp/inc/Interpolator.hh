#ifndef CYLINTERP_INTERPOLATOR_HH
#define CYLINTERP_INTERPOLATOR_HH

#include "Geometry.hh"
#include "Tools.hh"

class Interpolator : public UniformCylindricalGrid{
    public:
        Interpolator(double r0, double r1,
            double z0, double z1,
            int nr, int nz,
            int n_first_ring,
            double *map_values, int map_dim)
             : UniformCylindricalGrid(r0, r1, z0, z1, nr, nz, n_first_ring){
            
            m_dpmap_values = map_values;
            m_imap_dim = map_dim;
        };
        
        ~Interpolator();

        void Interpolate(double points[],
            double cart_points[],
            double interp_values[], int i, int i_out);
    
    // protected:
        double *m_dpmap_values;
        int m_imap_dim;
};

#endif