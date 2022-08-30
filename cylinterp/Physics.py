import cylinterp
import numpy as np
import pandas as pd
import itertools
import nestpy
from scipy.interpolate import interp1d

class Field(cylinterp.Geometry.UniformCylindricalGrid):
    def __init__(self, r0, r1, z0, z1, nr, nz, n_first_ring,
                 file, sag, n_cath_wires, r_max_det):

        """
        Takes in the result of COMSOL according to a grid made by the
        UniformCylindricalGrid class and provides an interpolation function
        using barycentric coordinates

        :param r0: Lower radius
        :param r1: Upper radius
        :param z0: Lower z
        :param z1: Upper z
        :param nr: Number of radii values
        :param nz: Number of z slices
        :param n_first_ring: Number of sections of the first ring
        :param file: The COMSOL file for the field map
        :param sag:
        :param n_cath_wires:
        :param r_max_det:
        """

        super().__init__(r0, r1, z0, z1, nr, nz, n_first_ring)
        self.Emap = pd.read_csv(file, sep=' ', header=None)
        self.Emap = self.Emap.rename(columns=dict(zip(self.Emap.columns, ['x', 'y', 'z', 'Enorm', 'Ex', 'Ey', 'Ez'])))
        self.Emap['Enorm'] = 1e3 * self.Emap['Enorm']
        self.Emap['Ex'] = 1e3 * self.Emap['Ex']
        self.Emap['Ey'] = 1e3 * self.Emap['Ey']
        self.Emap['Ez'] = 1e3 * self.Emap['Ez']
        self.Evec = self.Emap[['Ex', 'Ey', 'Ez']].values

    def Interpolate(self, points):
        tetra_indices = self.TetraIndices(points)
        cart_tetra = Tools.to_cartesian(self.polar_cs_grid[tetra_indices])
        cart_points = Tools.to_cartesian(points)

        A = cart_tetra[:, 0, :]
        B = cart_tetra[:, 1, :]
        C = cart_tetra[:, 2, :]
        D = cart_tetra[:, 3, :]
        BA = B - A
        CA = C - A
        DA = D - A
        PA = cart_points - A

        c1 = Tools.ScalarTripleProduct(PA, CA, DA) / Tools.ScalarTripleProduct(BA, CA, DA)
        c2 = Tools.ScalarTripleProduct(PA, BA, DA) / Tools.ScalarTripleProduct(CA, BA, DA)
        c3 = Tools.ScalarTripleProduct(PA, BA, CA) / Tools.ScalarTripleProduct(DA, BA, CA)
        c0 = 1 - c1 - c2 - c3

        c0 = c0.reshape(len(c0), 1)
        c1 = c1.reshape(len(c1), 1)
        c2 = c2.reshape(len(c2), 1)
        c3 = c3.reshape(len(c3), 1)

        return (c0 * self.Evec[tetra_indices[:, 0]] +
                c1 * self.Evec[tetra_indices[:, 1]] +
                c2 * self.Evec[tetra_indices[:, 2]] +
                c3 * self.Evec[tetra_indices[:, 3]])

nc = nestpy.NESTcalc(nestpy.VDetector())
@np.vectorize
def nest_vd(field, temp=177.15 , density = 2.94):
    dv = nc.SetDriftVelocity(temp, # K temp
                            density,
                            field)
    return dv*(1e-1) #cm/us

egridvd = np.linspace(0,1e6,int(1e6))
vd_interp = interp1d(egridvd, nest_vd(egridvd))

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

    def Drift(self,
              r=None,
              n_pts=100,
              dt=10e-3,
              recursion_limit=3000,
              q=-1,
              cath_thresh=200e-4,
              driftregion=None,
              sampleregion=None,
              tracking=True, **kwargs):

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

        if r == None:
            r = self.generate_rand_points(n_pts, zmin=sampleregion['zmin'],
                                          zmax=sampleregion['zmax'],
                                          rmin=sampleregion['rmin'],
                                          rmax=sampleregion['rmax'])
        points = r
        points_cartesian = Tools.to_cartesian(points)

        ids = np.arange(n_pts)
        remaining_ids = ids
        ended = np.repeat(False, n_pts)
        end_time = np.zeros(n_pts)
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

                E_interp = self.Interpolate(points)

                #The magnitude of the electric field at each point
                E_interp_norm = np.linalg.norm(E_interp, axis=1)

                #The drift velocity for each time, dl is in the direction of the field
                v = self.vd(E_interp_norm)
                dl = (v * dt * E_interp.T / E_interp_norm).T

                points_cartesian = Tools.to_cartesian(points)
                #Step forward by dl
                points_cartesian = points_cartesian + q * dl
                points = Tools.to_polar(points_cartesian)

                # Flag for successful drift to the anode
                hit_anode = points.T[1] <= driftregion['rmin']

                # Flag for hitting the cathode
                cath_pos = self.cathode_position(points.T[0])
                reshapen_cart = np.tile(points_cartesian[:, :2], self.n_cath_wires).reshape(cath_pos.shape)
                dist_to_cath = np.linalg.norm(reshapen_cart - cath_pos, axis=2)
                hit_cath = np.any(dist_to_cath < cath_thresh, axis=1)

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
                points = points[~ended]
                ended_ids = remaining_ids[ended]
                end_time[ended_ids] += (i + 1) * dt
                #The remaining indices that haven't ended
                remaining_ids = remaining_ids[~ended]

                z_tracks[i + 1][remaining_ids] = points.T[0]
                r_tracks[i + 1][remaining_ids] = points.T[1]
                theta_tracks[i + 1][remaining_ids] = points.T[2]

            return end_time, status, z_tracks.T, r_tracks.T, theta_tracks.T, r

        else:
            for i in range(recursion_limit):
                if len(points) == 0:
                    break
                E_interp = self.Interpolate(points)
                E_interp_norm = np.linalg.norm(E_interp, axis=1)

                v = self.vd(E_interp_norm)
                v[v < 0] = 0
                dl = (v * dt * E_interp.T / E_interp_norm).T

                points_cartesian = Tools.to_cartesian(points)
                points_cartesian = points_cartesian + q * dl
                points = Tools.to_polar(points_cartesian)


                # Flag for successful drift to the anode
                hit_anode = points.T[1] <= driftregion['rmin']

                # Flag for hitting the cathode
                cath_pos = self.cathode_position(points.T[0])
                reshapen_cart = np.tile(points_cartesian[:, :2], self.n_cath_wires).reshape(cath_pos.shape)
                dist_to_cath = np.linalg.norm(reshapen_cart - cath_pos, axis=2)
                hit_cath = np.any(dist_to_cath < cath_thresh, axis=1)

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

                points = points[~ended]
                ended_ids = remaining_ids[ended]
                end_time[ended_ids] += (i + 1) * dt
                remaining_ids = remaining_ids[~ended]
            return end_time, status, r