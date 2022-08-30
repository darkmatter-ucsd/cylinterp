import cylinterp
import numpy as np
import pandas as pd
import itertools

class UniformCylindricalGrid():
    """
    Defines a cylindrical grid with roughly equal volumes. This is a controlled
    grid which lets us index the field points quicker
    """

    def __init__(self, r0, r1, z0, z1, nr, nz, n_first_ring):
        
        """
        Gives a uniform cylindrical grid as well as indexing functions for fast access to grid points
        
        :param r0: Lower radius
        :param r1: Upper radius
        :param z0: Lower z
        :param z1: Upper z
        :param nr: Number of radii values
        :param nz: Number of z slices
        :param n_first_ring: Number of sections of the first ring
        """
        
        self.r0 = r0
        self.r1 = r1
        self.z0 = z0
        self.z1 = z1
        self.nr = nr
        self.nz = nz
        self.n_first_ring = n_first_ring

        # Set up the radial points, grid means edges, coords means the actual evaluation points
        self.rgrid = np.linspace(r0, r1, nr)
        self.r_coords, self.theta_coords, self.n_theta_in_r_slice, self.d_theta_in_r_slice, self.cum_n_theta_below_r = self.MakePolarCoords()
        self.unique_r_coords = (self.rgrid[1:] + self.rgrid[:-1]) / 2
        self.dr = self.rgrid[1] - self.rgrid[0]

        # Set up the z coordinates
        self.z_coords = np.linspace(z0, z1, nz)
        self.dz = self.z_coords[1] - self.z_coords[0]

        # Set the COMSOL Grids
        self.cs_grid = self.ComsolGrid()
        self.polar_cs_grid = cylinterp.Tools.to_polar(self.cs_grid.values)
        self.cart_cs_grid = self.cs_grid.values

        # Some helpful tools to save to memory
        self.tetra_iter = list(itertools.product(range(2), repeat=3))

    def MakePolarCoords(self):

        """
        Returns the info that will be useful for indexing a cylindrical grid.
        This first divides the radii into slices in r, then each slice in r is
        divided into roughly equal area bins, which defines the theta bins.


        :return: r_coords: the points at which a radius is evaluated
        theta_coords: the points of the evaluated thetas
        n_theta_in_r_slice: the number of theta coordinates in a particular slice of r
        d_theta_in_r_slice: the dtheta for each slice in r
        tot_n_theta_below_r: the cumulative number of the theta below the slice in r
        """

        #Find the small area of each angular segment within the first ring (from 0 to dr)
        da = np.pi * (self.rgrid[1] ** 2 - self.rgrid[0] ** 2) / self.n_first_ring
        r_coords = []
        theta_coords = []
        n_theta_in_r_slice = [0]
        d_theta_in_r_slice = []

        for i in range(len(self.rgrid) - 1):
            #Loop over each ring (i.e. step in dr)
            r_center = (self.rgrid[i + 1] + self.rgrid[i]) / 2
            area_ring = np.pi * (self.rgrid[i + 1] ** 2 - self.rgrid[i] ** 2)
            
            #the number of angular segments is approximately the area of the ring divided by the initial da
            n_theta = int(np.round(area_ring / da))
            
            #for indexing, we also need to know how many theta points are within each r-slice
            n_theta_in_r_slice.append(n_theta)
            #as well as d_theta
            d_theta_in_r_slice.append(2 * np.pi / n_theta)
            
            #and the theta coordinates for each r-slice as well
            theta_ranges = np.linspace(0, 2 * np.pi * (n_theta - 1) / n_theta, n_theta)
            r_coords.append(np.repeat(r_center, n_theta))
            theta_coords.append(theta_ranges)

        r_coords = np.concatenate(r_coords)
        theta_coords = np.concatenate(theta_coords)
        n_theta_in_r_slice = np.array(n_theta_in_r_slice)
        d_theta_in_r_slice = np.array(d_theta_in_r_slice)
        
        #It's also useful to use 1-d indices by considering the total number of theta values below a certain radius
        tot_n_theta_below_r = np.cumsum(n_theta_in_r_slice)

        return r_coords, theta_coords, n_theta_in_r_slice, d_theta_in_r_slice, tot_n_theta_below_r[:-1]

    def ComsolGrid(self):

        """
        Gives a cartesian grid for COMSOL to evaluate at

        :return: Cartesian grid for COMSOL to evaluate at
        """

        x_coords = self.r_coords * np.cos(self.theta_coords)
        y_coords = self.r_coords * np.sin(self.theta_coords)

        X_coords = np.tile(x_coords, len(self.z_coords))
        Y_coords = np.tile(y_coords, len(self.z_coords))
        Z_coords = np.repeat(self.z_coords, len(x_coords))

        cyl_coords = pd.DataFrame({'x': X_coords,
                                   'y': Y_coords,
                                   'z': Z_coords})

        return cyl_coords

    def ExportGrid(self, outfilename):

        """
        Creates a coordinate file to export for COMSOL

        :param outfilename: Name of the output file
        :return: None
        """

        comsol_grid = self.ComsolGrid()
        comsol_grid.to_csv(outfilename, sep=' ', header=False, index=False)

    def NearestIndex(self, points):
        """
        Finds the grid indices closest to your points of interest

        :param points: a (n_points, n_dim) array of (z,r,theta) coordinates that you want to find the closest point to
        :return: the (flat) index in the coordinate grid
        """

        z_ind = (np.round((points.T[0] - self.z0) / self.dz)).astype('int')
        r_ind = (np.round((points.T[1] - self.r_coords[0]) / self.dr)).astype('int')
        r_offset = self.cum_n_theta_below_r[np.clip(r_ind, 0, len(self.r_coords)).astype('int')]
        theta_ind = (np.round(
            points.T[2] / self.d_theta_in_r_slice[np.clip(r_ind, 0, len(self.r_coords)).astype('int')])).astype('int')

        flat_ind = z_ind * len(self.r_coords) + r_offset + theta_ind

        return flat_ind

    def LowerBinIndex(self, points):
        """
        Finds the bin indices closest to your points of interest.
        Bin indices are the closest values of r, theta, z that are on the grid but less than the
        r, theta, and z of interest (i.e. the input points)

        :param points: a (n_points, n_dim) array of (z,r,theta) coordinates
        :return: the (flat) indices in the field map
        """
        z_ind = ((points.T[0] - self.z0) / self.dz).astype('int')
        
        #Find which ring you're on
        r_ind = ((points.T[1] - self.r_coords[0]) / self.dr).astype('int')
        r_ind = np.clip(r_ind, 0, len(self.r_coords) - 1).astype('int')
        
        #Find the number of theta values below the ring's radius
        r_offset = self.cum_n_theta_below_r[r_ind]
        
        #Find the index of theta within the ring
        theta_ind = ((points.T[2] / self.d_theta_in_r_slice[r_ind])).astype('int')
        
        #z_ind needs to be multiplied by the number of coordinates in each z-sliced disk
        flat_ind = z_ind * len(self.r_coords) + r_offset + theta_ind
        
        #returns the 1-d index (pretty efficient!)
        return flat_ind

    def BoxIndices(self, points):
        """
        Finds the eight points which encompass the point of interest in a box

        :param points: a (n_points, n_dim) array of (z,r,theta) coordinates
        :return: an array of the 8 nearest grid indices to the given points
        """

        z_ind = ((points.T[0] - self.z0) / self.dz).astype('int')

        r_ind = ((points.T[1] - self.r_coords[0]) / self.dr).astype('int')
        # r_ind = np.clip(r_ind,0,len(self.r_coords)).astype('int')

        indices = np.zeros((8, len(points)))
        for i, t in enumerate(self.tetra_iter):
            # t[0]: Whether to stay on the current ring or jump to the next one
            # t[1]: '' theta
            # t[2]: '' z
            r_offset = self.cum_n_theta_below_r[r_ind + t[0]]
            theta_ind = ((points.T[2] / self.d_theta_in_r_slice[r_ind + t[0]]).astype('int'))

            # If you add dtheta to your current theta, and the new theta goes out of bounds, then you reset to the initial point
            theta_reset_mask = ((self.d_theta_in_r_slice[r_ind + t[0]] + points.T[2] < 2 * np.pi) & (t[1] == 1)) | (
                        t[1] == 0)
            indices[i] = (z_ind + t[2]) * len(self.r_coords) + r_offset + (theta_ind + t[1]) * theta_reset_mask

        return indices.astype('int').T

    def TetraIndices(self, points):
        """
        Gives the indices of the tetrahedron that encompasses your point in question.
        Sometimes the actual point is outside of the tetrahedron, but this usually
        means that the point is very close to one of the faces anyway. Which means for the
        interpolation, it doesn't matter. However I do need to fix this at some point...
        
        :param points: a (n_points, n_dim) array of (z,r,theta) coordinates
        :return: an array of the 4 grid indices that form a tetrahedron around the given points
        """
        tetra_indices = np.zeros((4, len(points)))
        z_ind_float = (points.T[0] - self.z0) / self.dz
        r_ind_float = (points.T[1] - self.r_coords[0]) / self.dr

        z_ind_nearest = np.round(z_ind_float).astype('int')
        z_ind_bin = z_ind_float.astype('int')
        r_ind_nearest = np.round(r_ind_float).astype('int')
        r_ind_bin = r_ind_float.astype('int')

        r_nearest_offset = self.cum_n_theta_below_r[r_ind_nearest]

        theta_ind_float = points.T[2] / self.d_theta_in_r_slice[r_ind_nearest]
        theta_ind_nearest = np.round(theta_ind_float).astype('int')
        theta_ind_bin_0 = theta_ind_float.astype('int')
        theta_ind_bin_1 = theta_ind_bin_0 + 1
        # Need a mask for when the next index will bring the theta past 2pi
        theta_inb_mask = ((points.T[2] + self.d_theta_in_r_slice[r_ind_nearest]) < 2 * np.pi)
        theta_ind_bin_1 *= theta_inb_mask

        # If the nearest index is greater than the bin index (i.e. it's closer to the next bin AND it's oob, make the nearest theta index 0
        theta_ind_nearest *= (~((theta_ind_nearest > theta_ind_bin_0) & (~theta_inb_mask)))

        # Two points closest to the closest ring
        tetra_indices[0] = z_ind_nearest * len(self.r_coords) + r_nearest_offset + theta_ind_bin_0
        tetra_indices[1] = z_ind_nearest * len(self.r_coords) + r_nearest_offset + theta_ind_bin_1

        # Now we will choose to step either backwards or forwards in r
        r_step = 1 - 2 * (r_ind_nearest - r_ind_bin)
        r_ind_stepped = r_ind_nearest + r_step

        r_step_offset = self.cum_n_theta_below_r[r_ind_stepped]

        theta_ind_nearest_rstep = np.round(points.T[2] / self.d_theta_in_r_slice[r_ind_stepped]).astype('int')
        theta_inb_rstep_mask = ((points.T[2] + self.d_theta_in_r_slice[r_ind_stepped]) < 2 * np.pi).astype('int')
        theta_ind_nearest_rstep *= theta_inb_rstep_mask

        tetra_indices[2] = z_ind_nearest * len(self.r_coords) + r_step_offset + theta_ind_nearest_rstep

        # Lastly we need to choose a step in z
        z_step = 1 - 2 * (z_ind_nearest - z_ind_bin)
        tetra_indices[3] = (z_ind_nearest + z_step) * len(self.r_coords) + r_nearest_offset + theta_ind_nearest

        return tetra_indices.astype('int').T

    def generate_rand_points(self,
                             n_pts,
                             rmin=None,
                             rmax=None,
                             zmin=None,
                             zmax=None):
        """
        Generates random points uniformly along a cylindrical annulus.
        Final answer is in cylindrical coordinates of (z,r,theta)
        
        :param n_pts: Number of points to generate
        :param rmin: Minimum radius to draw from
        :param rmax: Maximum radius to draw from
        :param zmin: Minimum z to draw from
        :param zmax: Maximum z to draw from
        :return: 
        """

        if (rmin == None):
            rmin = self.unique_r_coords[0]
        if (rmax == None):
            rmax = self.unique_r_coords[-1]
        if (zmin == None):
            zmin = self.z_coords[0]
        if (zmax == None):
            zmax = self.z_coords[-1]

        rand_z = np.random.uniform(zmin, zmax, n_pts)
        rand_r = np.sqrt(np.random.uniform(rmin ** 2, rmax ** 2, n_pts))
        rand_theta = np.random.uniform(0, 2 * np.pi, n_pts)
        rand_points = np.array([rand_z, rand_r, rand_theta]).T
        return rand_points