#include "Geometry.hh"

UniformCylindricalGrid::UniformCylindricalGrid(double r0, double r1,
    double z0, double z1,
    int nr, int nz,
    int n_first_ring) {
    
    m_dr0 = r0;
    m_dr1 = r1;
    m_dz0 = z0;
    m_dz1 = z1;
    m_inr = nr;
    m_inz = nz;
    m_in_first_ring = n_first_ring;

    MakePolarCoords();
    m_ddr = m_dr_grid[1] - m_dr_grid[0];
    for (int i=0; i<m_dr_grid.size()-1; i++){
        m_dunique_r_coords.emplace_back((m_dr_grid[i] + m_dr_grid[i+1])/2);
    }

    m_dz_coords = Linspace(m_dz0, m_dz1, m_inz);
    m_ddz = m_dz_coords[1]-m_dz_coords[0];

    SetGrid();
}

UniformCylindricalGrid::~UniformCylindricalGrid(){
    delete m_dpolar_cs_grid;
    delete m_dcart_cs_grid;
}

void UniformCylindricalGrid::MakePolarCoords(){
    m_dr_grid = Linspace(m_dr0, m_dr1, m_inr);

    double da = M_PI*(std::pow(m_dr_grid[1], 2) - std::pow(m_dr_grid[0], 2))/m_in_first_ring;
    m_in_theta_in_r_slice.emplace_back(0);
    for (int i=0; i<m_inr-1; i++){
        //Loop over each ring (i.e. step in dr)
        double r_center = (m_dr_grid[i+1]+m_dr_grid[i])/2;
        double area_ring = M_PI * (std::pow(m_dr_grid[i+1],2)-std::pow(m_dr_grid[i],2));

        //The number of angular segments is approximately the area of the ring divided by the initial da
        int n_theta = std::round(area_ring/da);

        //For indexing, we also need to know how many theta points are within each r-slice
        m_in_theta_in_r_slice.emplace_back(n_theta);
        //...as well as d_theta
        m_dd_theta_in_r_slice.emplace_back(2*M_PI/n_theta);

        //...and the theta coordinates for each r-slice as well
        std::vector<double> theta_ranges = Linspace(0., 2*M_PI*(n_theta-1)/n_theta, n_theta);
        for (int j=0; j<n_theta; j++){
            m_dr_coords.emplace_back(r_center);
            m_dtheta_coords.emplace_back(theta_ranges[j]);
        }
    }

    int cumulative = 0;
    for (int i=0; i<m_inr-1; i++){
        cumulative+=m_in_theta_in_r_slice[i];
        m_icum_n_theta_below_r.emplace_back(cumulative);
    }
}

void UniformCylindricalGrid::SetGrid(){
    int total_grid_size = m_dr_coords.size()*m_dz_coords.size();

    m_dpolar_cs_grid = new double[3*total_grid_size];
    m_dcart_cs_grid = new double[3*total_grid_size];

    for (int i=0; i<m_inz; i++){
        for (int j=0; j<m_dr_coords.size(); j++){
            int event_index = 3*(m_dr_coords.size()*i+j);
            m_dcart_cs_grid[event_index] = m_dr_coords[j]*std::cos(m_dtheta_coords[j]);
            m_dcart_cs_grid[event_index+1] = m_dr_coords[j]*std::sin(m_dtheta_coords[j]);
            m_dcart_cs_grid[event_index+2] = m_dz_coords[i];

            m_dpolar_cs_grid[event_index] = m_dz_coords[i];
            m_dpolar_cs_grid[event_index+1] = m_dr_coords[j];
            m_dpolar_cs_grid[event_index+2] = m_dtheta_coords[j];
        }
    }
}

void UniformCylindricalGrid::TetraIndices(double points[], int tetra_indices[], int i, int i_out){
    /*
    :param points: a 1-d array with length 3 x n_points arranged like
        [z0, r0, phi0, ..., zN-1, rN-1, phiN-1]
    :param tetra_indices: a 1-d array with length 4 x n of indices of the
        grid points which lie in a tetrahedron around the point of interest
    :param i: the index of the point of interest 
    */

    double z_ind_float = (points[3*i] - m_dz0)/m_ddz;
    double r_ind_float = (points[3*i+1] - m_dr_coords[0])/m_ddr;

    int z_ind_nearest = std::round(z_ind_float);
    int z_ind_bin = int(z_ind_float);
    int r_ind_nearest = std::round(r_ind_float);
    int r_ind_bin = int(r_ind_float);

    int r_nearest_offset = m_icum_n_theta_below_r[r_ind_nearest];

    double theta_ind_float = points[3*i+2]/m_dd_theta_in_r_slice[r_ind_nearest];
    int theta_ind_nearest = std::round(theta_ind_float);
    int theta_ind_bin_0 = int(theta_ind_float);
    int theta_ind_bin_1 = theta_ind_bin_0 + 1;

    // In case the next index will bring the theta past 2pi
    bool theta_inb = (points[3*i+2] + m_dd_theta_in_r_slice[r_ind_nearest]) < (2*M_PI);

    // If the nearest index is greater than the bin index (i.e. it's closer to the next bin AND it's oob, make the nearest theta index 0
    if (!theta_inb) theta_ind_bin_1 = 0;
    if (!(!((theta_ind_nearest > theta_ind_bin_0)&&(!theta_inb)))) theta_ind_nearest = 0;

    // Two points closest to the closest ring
    tetra_indices[4*i_out] = z_ind_nearest * m_dr_coords.size() + r_nearest_offset + theta_ind_bin_0;
    tetra_indices[4*i_out+1] = z_ind_nearest * m_dr_coords.size() + r_nearest_offset + theta_ind_bin_1;

    // Now we will choose to step either backwards or forwards in r
    int r_step = 1 - 2 * (r_ind_nearest - r_ind_bin);
    int r_ind_stepped = r_ind_nearest + r_step;

    int r_step_offset = m_icum_n_theta_below_r[r_ind_stepped];

    int theta_ind_nearest_rstep = std::round(points[3*i+2]/m_dd_theta_in_r_slice[r_ind_stepped]);
    if (!((points[3*i+2] + m_dd_theta_in_r_slice[r_ind_stepped])< (2*M_PI)))
        theta_ind_nearest_rstep = 0;

    tetra_indices[4*i_out+2] = z_ind_nearest * m_dr_coords.size() + r_step_offset + theta_ind_nearest_rstep;

    // Lastly we need to choose a step in z
    int z_step = 1 - 2*(z_ind_nearest - z_ind_bin);
    tetra_indices[4*i_out+3] = (z_ind_nearest + z_step) * m_dr_coords.size() + r_nearest_offset + theta_ind_nearest;
}

void UniformCylindricalGrid::ExportClass(){
    std::ofstream CartGridOut("cartesian_points.bin", std::ios::out | std::ios::binary);
    std::ofstream PolarGridOut("polar_points.bin", std::ios::out | std::ios::binary);
    std::ofstream ZGridOut("z_points.bin", std::ios::out | std::ios::binary);
    std::ofstream RGridOut("r_points.bin", std::ios::out | std::ios::binary);
    std::ofstream ThetaGridOut("theta_points.bin", std::ios::out | std::ios::binary);

    CartGridOut.write((char*)&m_dcart_cs_grid[0], sizeof(double) * 3*m_dr_coords.size() * m_dz_coords.size());
    PolarGridOut.write((char*)&m_dpolar_cs_grid[0], sizeof(double) * 3*m_dr_coords.size() * m_dz_coords.size());
    ZGridOut.write((char*)&m_dz_coords[0], sizeof(double) * m_dz_coords.size());
    RGridOut.write((char*)&m_dr_coords[0], sizeof(double) * m_dr_coords.size());
    ThetaGridOut.write((char*)&m_dtheta_coords[0], sizeof(double) * m_dtheta_coords.size());

    CartGridOut.close();
    PolarGridOut.close();
    ZGridOut.close();
    RGridOut.close();
    ThetaGridOut.close();
}
