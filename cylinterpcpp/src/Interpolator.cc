#include "Interpolator.hh"

Interpolator::~Interpolator(){
    delete[] m_dpmap_values;
}

void Interpolator::Interpolate(double points[], double cart_points[], double interp_values[], int i, int i_out){
    // Find the tetrahedral indices
    int tetra_indices[4];
    TetraIndices(points, tetra_indices, i, 0);
    
    // Cartesian coordinate of the point of interest
    // double cart_points[3];
    ToCartesian(points, cart_points, i, i);

    // Find the cartesian coordinates of those tetrahedral indices
    double cart_tetra_points[12];
    for (int j=0; j<4; j++){
        for (int dim_j=0; dim_j<3; dim_j++){
            cart_tetra_points[3*j+dim_j] = m_dcart_cs_grid[3*tetra_indices[j]+dim_j];
        }
    }

    // Subtract off the first corner
    for (int k=1; k<4; k++){
        for (int l=0; l<3; l++){
            cart_tetra_points[3*k+l] -= cart_tetra_points[l];
        }
    }

    // Subtract off the first corner of the tetrahedron from the point in question
    double cart_points_shifted[3];
    for (int m=0; m<3; m++)
        cart_points_shifted[m] = cart_points[3*i+m] - cart_tetra_points[m];

    // Barycentric coefficients
    double tetra_stp = ScalarTripleProduct(&cart_tetra_points[3], &cart_tetra_points[6], &cart_tetra_points[9], 0);
    double bary_c[4];
    bary_c[1] = ScalarTripleProduct(&cart_points_shifted[0], &cart_tetra_points[6], &cart_tetra_points[9], 0)/tetra_stp;
    bary_c[2] = -ScalarTripleProduct(&cart_points_shifted[0], &cart_tetra_points[3], &cart_tetra_points[9], 0)/tetra_stp;
    bary_c[3] = ScalarTripleProduct(&cart_points_shifted[0], &cart_tetra_points[3], &cart_tetra_points[6], 0)/tetra_stp;
    bary_c[0] = 1-bary_c[1]-bary_c[2]-bary_c[3];

    for (int d=0; d<m_imap_dim; d++){
        interp_values[m_imap_dim*i_out + d]=0;
        for (int bc=0; bc<4; bc++){
            interp_values[m_imap_dim*i_out + d]+=bary_c[bc]*m_dpmap_values[m_imap_dim*tetra_indices[bc]+d];
        }
    }
}