#include "Physics.hh"

RTPC::RTPC(Interpolator* field,
            std::vector<double> &vd_E,
            std::vector<double> &vd_vd,
            double sag,
            int n_cath_wires,
            double r_max_det){

    m_field = field;
    m_dvd_E = vd_E;
    m_dvd_vd = vd_vd;
    m_dsag = sag;
    m_in_cath_wires = n_cath_wires;
    m_dr_max_det = r_max_det;

    m_dcath_dth = 2*M_PI/m_in_cath_wires;

    for (int i=0; i<m_in_cath_wires; i++){
        m_dcath_angles.emplace_back(2*M_PI*i/m_in_cath_wires);
    }

    if ((m_dvd_E[1]-m_dvd_E[0])==(m_dvd_E[m_dvd_E.size()-1]-m_dvd_E[m_dvd_E.size()-2]))
        m_bregular_vd_grid = true;
        m_ddE = m_dvd_E[1]-m_dvd_E[0];
}

RTPC::~RTPC(){};

double RTPC::NearestCathodeDistance(double polar_points[], double cart_points[], int i){
    double zmid = (m_field->m_dz_coords[0]+m_field->m_dz_coords[m_field->m_inz-1])/2;
    double r = m_dsag * (polar_points[3*i]-m_field->m_dz1) * (polar_points[3*i]-m_field->m_dz0);
    r = r/((zmid-m_field->m_dz1) * (zmid-m_field->m_dz0));
    r = m_dr_max_det - r;

    int th_ind = std::round(polar_points[3*i+2]/m_dcath_dth);
    if ((polar_points[3*i+2]>2*M_PI-m_dcath_dth/2)||(std::isnan(polar_points[3*i+2])))
        th_ind = 0;
    double cath_pos[2] = {r*std::cos(m_dcath_angles[th_ind]),
        r*std::sin(m_dcath_angles[th_ind])};

    return std::sqrt((cart_points[3*i]-cath_pos[0])*(cart_points[3*i]-cath_pos[0]) +
        (cart_points[3*i+1]-cath_pos[1])*(cart_points[3*i+1]-cath_pos[1]));
}

double RTPC::VdInterp(double E_norm) {
    if (E_norm < m_dvd_E[0]) return 0;
    else if (E_norm > m_dvd_E[m_dvd_E.size()-1]) return m_dvd_vd[m_dvd_vd.size()-1];
    else{
        int vd_ind;
        if (m_bregular_vd_grid){
            vd_ind = int((E_norm - m_dvd_E[0])/m_ddE);
        }
        else{
            auto it = std::upper_bound(m_dvd_E.begin(), m_dvd_E.end(), E_norm);
            vd_ind = std::distance(m_dvd_E.begin(), it) -1;
        }

        double vd = (m_dvd_vd[vd_ind+1]-m_dvd_vd[vd_ind])/(m_dvd_E[vd_ind+1]-m_dvd_E[vd_ind]);
        vd = vd*(E_norm - m_dvd_E[vd_ind]) + m_dvd_vd[vd_ind];
        return vd;
    }
}

double RTPC::LongDiffusion(double E_norm){
    // E_norm MUST be in V/cm!!!
    double long_diff = m_dlongdiff_C+m_dlongdiff_A*std::exp(-(E_norm - m_dlongdiff_E0)/m_dlongdiff_B);
    if (long_diff>m_dlongdiff_C+m_dlongdiff_A) return (m_dlongdiff_C+m_dlongdiff_A)*cm*cm/s;
    else if (long_diff < m_dlongdiff_C) return m_dlongdiff_C*cm*cm/s;
    else return long_diff*cm*cm/s;

}

void RTPC::Drift(double r[], int n_pts,
    double dr_zmin, double dr_zmax,
    double dr_rmin, double dr_rmax,
    bool tracking, bool diffusion, pcg32& rng,
    double end_time[], double end_pos[], int status[],
    double tracks[]){

    //Load cartesian point buffer
    double *cart_r = new double[3*n_pts];
    //Load electric field buffer
    double *E_interp = new double[3*n_pts];

    for (int p=0; p<n_pts; p++){
        for (int c=0; c<m_irecursion_limit; c++){
            if (tracking){
                for (int t_d=0; t_d<3; t_d++) {
                    tracks[p*3*m_irecursion_limit + 3*c +t_d] = r[3*p+t_d];
                }
            }

            //This is a 2-in-1 step, calculate the electric field...
            //and the cartesian point
            m_field->Interpolate(r, cart_r, E_interp, p, p);
            double E_interp_norm = std::sqrt(E_interp[3*p]*E_interp[3*p]+
                E_interp[3*p+1]*E_interp[3*p+1]+
                E_interp[3*p+2]*E_interp[3*p+2]);
            
            if (std::isnan(E_interp_norm)){
                status[p] = 2;
                end_time[p] = (c+1)*m_ddt;
                for (int ep_i=0; ep_i<3; ep_i++)
                    end_pos[3*p+ep_i] = r[3*p+ep_i];
                break;
            }

            double v_p = VdInterp(E_interp_norm);
            double dl_normvec[3] = {E_interp[3*p]/E_interp_norm,
                E_interp[3*p+1]/E_interp_norm,
                E_interp[3*p+2]/E_interp_norm};
            
            if (diffusion){
                double v_trans_1[3] = {0, 1, -dl_normvec[1]/dl_normvec[2]};
                double v_trans_1_norm = std::sqrt(v_trans_1[0]*v_trans_1[0]+
                    v_trans_1[1]*v_trans_1[1]+
                    v_trans_1[2]*v_trans_1[2]);

                for (int vt_i=0; vt_i<3; vt_i++)
                    v_trans_1[vt_i] = v_trans_1[vt_i]/v_trans_1_norm;

                double v_trans_2[3];
                CrossProduct(dl_normvec, v_trans_1, v_trans_2, 0, 0);

                std::normal_distribution<double> long_norm_dist(v_p*m_ddt, std::sqrt(2*LongDiffusion(E_interp_norm)*m_ddt));
                std::normal_distribution<double> trans_norm_dist(0., std::sqrt(2*m_dlxe_trans_diff*m_ddt));

                double long_step = long_norm_dist(rng);
                double trans_step_1 = trans_norm_dist(rng);
                double trans_step_2 = trans_norm_dist(rng);

                for (int d=0; d<3; d++){
                    cart_r[3*p+d] = cart_r[3*p+d]
                         - long_step*dl_normvec[d]
                         - trans_step_1*v_trans_1[d]
                         - trans_step_2*v_trans_2[d];
                }
            }
            else {
                for (int d=0; d<3; d++){
                    cart_r[3*p+d] = cart_r[3*p+d] - v_p * m_ddt * dl_normvec[d];
                }
            }
            //Update the polar point with the current cartesian point
            ToPolar(r, cart_r, p, p);

            /*
            FLAGS TO STOP
            */

            bool hit_anode = (r[3*p+1]<dr_rmin);
            if (hit_anode){
                status[p] = 0;
                end_time[p] = (c+1)*m_ddt;
                for (int ep_i=0; ep_i<3; ep_i++)
                    end_pos[3*p+ep_i] = r[3*p+ep_i];
                break;
            }

            bool hit_cathode = NearestCathodeDistance(r, cart_r, p) < m_dcath_thresh;
            if (hit_cathode){
                status[p] = 1;
                end_time[p] = (c+1)*m_ddt;
                for (int ep_i=0; ep_i<3; ep_i++)
                    end_pos[3*p+ep_i] = r[3*p+ep_i];
                break;
            }

            bool oob = ((r[3*p+1]>dr_rmax)||(r[3*p]>dr_zmax)||(r[3*p]<dr_zmin));
            if (oob) {
                status[p] = 3;
                end_time[p] = (c+1)*m_ddt;
                for (int ep_i=0; ep_i<3; ep_i++)
                    end_pos[3*p+ep_i] = r[3*p+ep_i];
                break;
            }
        }
    }

    delete[] cart_r;
    delete[] E_interp;
}