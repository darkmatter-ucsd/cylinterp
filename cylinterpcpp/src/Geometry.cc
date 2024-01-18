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
        double r_center = (m_dr_grid[i+1]-m_dr_grid[i])/2;
        double area_ring = M_PI * (std::pow(m_dr_grid[i+1],2)-std::pow(m_dr_grid[i],2))

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
            m_dtheta_coords.emplace_back(theta_ranges[j])
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