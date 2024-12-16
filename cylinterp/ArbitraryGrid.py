import cylinterp
import numpy as np
import pandas as pd
import itertools
import os
from scipy.interpolate import interp1d
from scipy.spatial import KDTree


class MeshInterpolator():
    def __init__(mesh_points, interp_values, n_dim, **kwargs):

        """
        This is really just a repeat of the Interpolator class from Physics.py.
        The only difference is that there's no expectation of a cylindrical grid
        which kind of goes against the point of cylinterp but I'm too lazy to
        put this into a new package.

        :param mesh_points: The points
        :param interp_values: The value to be interpolated
        """

        self.mesh_points = mesh_points
        self.interp_values = interp_values
        self.n_dim = n_dim


        self.mesh_kdtree = KDTree(self.mesh_points)
    
    
    def InterpolateArbitraryGrid(self, points,
                                 cart_points=None,
                                 k_n=4):
        """
        Interpolate points on an arbitrary grid using weighted inverse distance.
        Technically the argument of "points" isn't used, but we require the same
        format as above. Is it stupid? Yes. But I can't be bothered right now to
        figure out how to change it.
        """
        if type(cart_points)!=np.ndarray:
            raise ValueError("Argument 'cart_points' must be specified and must be a numpy.ndarray.")

        n_pts = len(cart_points)
        nearest_dist, nearest_ind = self.mesh_kdtree.query(cart_points, k=k_n)
        vec_nearest_dist = nearest_dist.repeat(self.n_dim).reshape(n_pts,k_n,self.n_dim)
        weighted_sum = np.sum(self.interp_values[nearest_ind]*(1/vec_nearest_dist), axis=1)
        sum_nearest_neighbors = np.sum(1/nearest_dist, axis=1).repeat(self.n_dim).reshape(n_pts,self.n_dim)
        interp_vec = weighted_sum/sum_nearest_neighbors

        return interp_vec