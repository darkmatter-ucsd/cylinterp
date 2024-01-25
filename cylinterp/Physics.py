import cylinterp
import numpy as np
import pandas as pd
import itertools
import os
from scipy.interpolate import interp1d

lxe_trans_diff = 55*(cylinterp.centimeters**2)/cylinterp.seconds #cm^2/s (EXO-200)
def lxe_long_diff(E, A=1.531e+01, B=5.361e+01, C=2.621e+01, E0=99.72156937):
    long_diff = C+A*np.exp(-(E-E0)/B)
    max_long_diff = C+A
    return np.clip(long_diff, C, max_long_diff)*(cylinterp.centimeters**2)/cylinterp.seconds

class Interpolator(cylinterp.Geometry.UniformCylindricalGrid):
    def __init__(self, r0, r1, z0, z1, nr, nz, n_first_ring):

        """
        Provides a template for the interpolator function

        :param r0: Lower radius
        :param r1: Upper radius
        :param z0: Lower z
        :param z1: Upper z
        :param nr: Number of radii values
        :param nz: Number of z slices
        :param n_first_ring: Number of sections of the first ring
        :param interp_values: The value to be interpolated
        """

        super().__init__(r0, r1, z0, z1, nr, nz, n_first_ring)
        self.interp_values = None
    
    def Interpolate(self, points, cart_points = None):
        tetra_indices = self.TetraIndices(points)
        cart_tetra = self.cart_cs_grid[tetra_indices]
        if type(cart_points)!=np.ndarray:
            cart_points = cylinterp.Tools.to_cartesian(points)

        A = cart_tetra[:, 0, :]
        B = cart_tetra[:, 1, :]
        C = cart_tetra[:, 2, :]
        D = cart_tetra[:, 3, :]
        BA = B - A
        CA = C - A
        DA = D - A
        PA = cart_points - A

        stp = cylinterp.Tools.ScalarTripleProduct(BA, CA, DA)
        c1 = cylinterp.Tools.ScalarTripleProduct(PA, CA, DA) / stp
        c2 = -cylinterp.Tools.ScalarTripleProduct(PA, BA, DA) / stp
        c3 = cylinterp.Tools.ScalarTripleProduct(PA, BA, CA) / stp
        c0 = 1 - c1 - c2 - c3

        c0 = c0.reshape(len(c0), 1)
        c1 = c1.reshape(len(c1), 1)
        c2 = c2.reshape(len(c2), 1)
        c3 = c3.reshape(len(c3), 1)

        return (c0 * self.interp_values[tetra_indices[:, 0]] +
                c1 * self.interp_values[tetra_indices[:, 1]] +
                c2 * self.interp_values[tetra_indices[:, 2]] +
                c3 * self.interp_values[tetra_indices[:, 3]])


class Field(Interpolator):
    def __init__(self, r0, r1, z0, z1, nr, nz, n_first_ring, file):

        """
        Takes in the result of COMSOL according to a grid made by the
        UniformCylindricalGrid class and provides an interpolation function
        using barycentric coordinates. The electric field MUST be exported in kV/cm,
        and is converted to V/cm inside the class

        :param r0: Lower radius
        :param r1: Upper radius
        :param z0: Lower z
        :param z1: Upper z
        :param nr: Number of radii values
        :param nz: Number of z slices
        :param n_first_ring: Number of sections of the first ring
        :param file: The COMSOL file for the field map
        """

        super().__init__(r0, r1, z0, z1, nr, nz, n_first_ring)
        self.Emap = pd.read_csv(file, sep=' ', header=None)
        self.Emap = self.Emap.rename(columns=dict(zip(self.Emap.columns, ['x', 'y', 'z', 'Enorm', 'Ex', 'Ey', 'Ez'])))
        self.Emap['Enorm'] = 1e3 * self.Emap['Enorm']
        self.Emap['Ex'] = 1e3 * self.Emap['Ex']
        self.Emap['Ey'] = 1e3 * self.Emap['Ey']
        self.Emap['Ez'] = 1e3 * self.Emap['Ez']
        self.Evec = self.Emap[['Ex', 'Ey', 'Ez']].values
        self.interp_values = self.Evec


class OpticalSimulation(Interpolator):
    def __init__(self, r0, r1, z0, z1, nr, nz, n_first_ring, n_initial_photons, file):

        """
        Takes in the result of an optical simulation

        :param r0: Lower radius
        :param r1: Upper radius
        :param z0: Lower z
        :param z1: Upper z
        :param nr: Number of radii values
        :param nz: Number of z slices
        :param n_first_ring: Number of sections of the first ring
        :param file: The COMSOL file for the field map
        """

        super().__init__(r0, r1, z0, z1, nr, nz, n_first_ring)
        self.n_initial_photons = n_initial_photons
        self.hitpatterns = np.load(file)
        self.interp_values = self.hitpatterns


# nc = nestpy.NESTcalc(nestpy.VDetector())
# @np.vectorize
# def nest_vd(field, temp=177.15 , density = 2.94):
#     dv = nc.SetDriftVelocity(temp, # K temp
#                             density,
#                             field)
#     return dv*(1e-1) #cm/us

# egridvd = np.linspace(0,1e6,int(1e6))
vd_E_grid_path = os.path.join(
    os.path.dirname(__file__),
    'data_files/vd_E_grid.npy')
vd_vd_grid_path = os.path.join(
    os.path.dirname(__file__),
    'data_files/vd_vd_grid.npy')

egridvd = np.load(vd_E_grid_path)
vdgrid = np.load(vd_vd_grid_path)
vd_interp = interp1d(egridvd, vdgrid)
# vd_interp = interp1d(egridvd, nest_vd(egridvd))

class RTPC(Field):

    def __init__(self, r0, r1, z0, z1, nr, nz, n_first_ring,
                 file, n_cath_wires, r_max_det,
                 sag = 0,
                 vd = vd_interp):

        """
        A type of field described by a geometry of cathode wires on the outer edge oriented along z

        :param r0: Lower radius
        :param r1: Upper radius
        :param z0: Lower z
        :param z1: Upper z
        :param nr: Number of radii values
        :param nz: Number of z slices
        :param n_first_ring: Number of sections of the first ring
        :param file: The COMSOL file for the field map
        :param n_cath_wires: Number of cathode wires
        :param r_max_det: Maximum radius of the detector
        :param sag: Sag of the cathodes in the middle of the detector
        :param vd: Drift velocity in cm/us
        """
        super().__init__(r0, r1, z0, z1, nr, nz, n_first_ring, file)
        self.sag = sag
        self.n_cath_wires = n_cath_wires
        self.cath_angles = 2 * np.pi * np.arange(n_cath_wires) / n_cath_wires
        self.cath_angles = self.cath_angles.reshape((n_cath_wires, 1))
        self.r_max_det = r_max_det
        self.vd = vd

    def cathode_position(self, z):
        """
        Assuming the sagging cathodes are parabolic, find the distance from the center of the detector
        to the cathode.

        :param z: A 1-d array of z-values
        :return: A 2-d array of [x,y] values for where the cathodes are located at each particular z
        """
        zmid = (self.z_coords[0] + self.z_coords[-1]) / 2
        r = self.r_max_det - (self.sag * (z - self.z_coords[-1]) * (z - self.z_coords[0]) /
                              ((zmid - self.z_coords[-1]) * (zmid - self.z_coords[0])))
        return np.array([r * np.cos(self.cath_angles), r * np.sin(self.cath_angles)]).T

    def nearest_cathode_pos(self, z, theta):
        """
        Assuming the sagging cathodes are parabolic, find the distance from the center of the detector
        to the cathode.

        :param z: A 1-d array of z-values
        :return: A 2-d array of [x,y] values for where the cathodes are located at each particular z
        """
        zmid = (self.z_coords[0] + self.z_coords[-1]) / 2
        r = self.r_max_det - (self.sag * (z - self.z_coords[-1]) * (z - self.z_coords[0]) /
                            ((zmid - self.z_coords[-1]) * (zmid - self.z_coords[0])))
        
        dth = 2*np.pi/self.n_cath_wires
        th_ind = np.round(theta/dth).astype('int')
        th_ind[(theta>2*np.pi - dth/2)|(np.isnan(theta))] = 0
        th_wires = self.cath_angles.flatten()[th_ind]
        return np.array([r*np.cos(th_wires), r*np.sin(th_wires)]).T

    def Drift(self,
              r=None,
              n_pts=100,
              dt=10e-3,
              recursion_limit=3000,
              q=-1,
              cath_thresh=200e-4,
              driftregion=None,
              sampleregion=None,
              tracking=True,
              diffusion=False, **kwargs):

        """
        Gives the average charge cloud path and drift time

        :param r: A 2-D array of polar coordinates [[z1, r1, th1],...,[zn, rn, thn]]
        :param n_pts: The number of starting positions to simulate
        :param dt: The time step for the drift in microseconds
        :param recursion_limit: The maximum number of loop iterations
        :param q: Charge of the particle
        :param cath_thresh: The maximum distance to the cathode without getting a 'hit cathode' flag
        :param driftregion: The region for the drift
        :param sampleregion: The region for which to sample points
        :param tracking: Whether or not to return the tracks of the electron cloud
        :param kwargs:
        :return: end_time of each electron cloud, status, x y and z tracks (optional), and initial starting positions
        """
        if driftregion == None:
            driftregion = {'zmin': self.z_coords[0],
                           'zmax': self.z_coords[-1],
                           'rmin': self.unique_r_coords[0],
                           'rmax': self.r_max_det}
        if sampleregion == None:
            sampleregion = driftregion

        if type(r) != np.ndarray:
            r = self.generate_rand_points(n_pts, zmin=sampleregion['zmin'],
                                          zmax=sampleregion['zmax'],
                                          rmin=sampleregion['rmin'],
                                          rmax=sampleregion['rmax'])
        points = r
        points_cartesian = cylinterp.Tools.to_cartesian(points)

        ids = np.arange(n_pts)
        remaining_ids = ids
        ended = np.repeat(False, n_pts)
        end_time = np.zeros(n_pts)
        end_pos = np.zeros((n_pts, 3))
        status = np.zeros(n_pts)

        #Statuses: 0 = hit anode, 1 = hit cathode, 2 = hit a region of NaN field, 3 = out of drift region

        if tracking:
            #I know that this isn't the best way to deal with tracking, but basically use one loop if we track
            #and use another loop if we don't track

            z_tracks = np.zeros((recursion_limit + 1, n_pts))
            r_tracks = np.zeros((recursion_limit + 1, n_pts))
            theta_tracks = np.zeros((recursion_limit + 1, n_pts))

            z_tracks[0] = r.T[0]
            r_tracks[0] = r.T[1]
            theta_tracks[0] = r.T[2]

        for i in range(recursion_limit):
            if len(points) == 0:
                break
            
            #points is in cylindrical at this point
            E_interp = self.Interpolate(points, cart_points = points_cartesian)

            #The magnitude of the electric field at each point
            E_interp_norm = np.linalg.norm(E_interp, axis=1)

            #The drift velocity for each time, dl is in the direction of the field
            v = self.vd(E_interp_norm)
            dl = (v * dt * E_interp.T / E_interp_norm).T
            
            if diffusion:
                #Set up the longitudinal directions
                dl_norm = np.linalg.norm(dl, axis = 1)
                
                #Set up the transverse directions
                v_trans_1 = np.array([np.zeros(len(dl)),
                                        np.ones(len(dl)),
                                        -dl.T[1]/dl.T[2]]).T
                
                v_trans_1 = v_trans_1/np.linalg.norm(v_trans_1, axis=1).reshape(len(dl), 1)
                v_trans_2 = np.cross(dl/dl_norm.reshape(len(dl), 1), v_trans_1)
                
                long_step = np.random.normal(dl_norm,
                                                np.sqrt(2*lxe_long_diff(E_interp_norm)*dt)).reshape(len(dl), 1)
                trans_step_1 = np.random.normal(0, np.sqrt(2*lxe_trans_diff*dt), len(dl)).reshape(len(dl), 1)
                trans_step_2 = np.random.normal(0, np.sqrt(2*lxe_trans_diff*dt), len(dl)).reshape(len(dl), 1)
                
                dl = long_step*dl/dl_norm.reshape(len(dl), 1)+trans_step_1*v_trans_1+trans_step_2*v_trans_2
            

            # points_cartesian = cylinterp.Tools.to_cartesian(points)
            #Step forward by dl
            points_cartesian = points_cartesian + q * dl
            points = cylinterp.Tools.to_polar(points_cartesian)

            # Flag for successful drift to the anode
            hit_anode = points.T[1] <= driftregion['rmin']

            # Flag for hitting the cathode
            # cath_pos = self.cathode_position(points.T[0])
            # reshapen_cart = np.tile(points_cartesian[:, :2], self.n_cath_wires).reshape(cath_pos.shape)
            # dist_to_cath = np.linalg.norm(reshapen_cart - cath_pos, axis=2)
            cath_pos = self.nearest_cathode_pos(points.T[0], points.T[2])
            dist_to_cath = np.linalg.norm(cath_pos-points_cartesian[:,:2], axis = 1)
            # hit_cath = np.any(dist_to_cath < cath_thresh, axis=1)
            hit_cath = dist_to_cath<cath_thresh

            # Flag for hitting a region of NaN
            nan_field = np.isnan(E_interp_norm)

            # Flag for out of bounds
            oob = ((points.T[1] > driftregion['rmax']) |
                    (points.T[0] > driftregion['zmax']) |
                    (points.T[0] < driftregion['zmin']))

            status[remaining_ids[hit_anode]] = 0
            status[remaining_ids[hit_cath]] = 1
            status[remaining_ids[nan_field]] = 2
            status[remaining_ids[oob]] = 3

            ended = hit_anode | hit_cath | nan_field | oob

            #Only keep the points which haven't ended their tracks
            ended_ids = remaining_ids[ended]
            end_pos[ended_ids] = points[ended]
            points = points[~ended]
            points_cartesian = points_cartesian[~ended]
            end_time[ended_ids] += (i + 1) * dt
            #The remaining indices that haven't ended
            remaining_ids = remaining_ids[~ended]

            if tracking:
                z_tracks[i + 1][remaining_ids] = points.T[0]
                r_tracks[i + 1][remaining_ids] = points.T[1]
                theta_tracks[i + 1][remaining_ids] = points.T[2]

        if tracking:
            return end_time, status, z_tracks.T, r_tracks.T, theta_tracks.T, r, end_pos
        else:
            return end_time, status, r, end_pos

