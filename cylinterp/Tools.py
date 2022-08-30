import numpy as np

class Tools():
    """
    A class that provides some useful math tools
    """
    @staticmethod
    def to_cartesian(polar_points):
        """
        A function that takes in a (n_pts, 3) array where each element represents a polar
        coordinate in the form of [z, r, theta] and returns those points in cartesian coordinates

        :param polar_points: A (n_pts, 3) array where each entry is a polar coord [z, r, theta]
        :return: A (n_pts, 3) array where each entry is a cartesian coord [x,y,z]
        """

        points_cartesian = np.array([polar_points.T[1] * np.cos(polar_points.T[2]),
                                     polar_points.T[1] * np.sin(polar_points.T[2]),
                                     polar_points.T[0]]).T
        return points_cartesian

    @staticmethod
    def to_polar(cartesian_points):
        """
        A function that takes in a (n_pts, 3) array where each element represents a cartesian
        coordinate in the form of [x,y,z] and returns those points in polar coordinates

        :param cartesian_points: A (n_pts, 3) array where each entry is a polar coord [x,y,z]
        :return: A (n_pts, 3) array where each entry is a cartesian coord [z, r, theta]
        """
        angles = np.arctan2(cartesian_points.T[1],
                            cartesian_points.T[0])
        angles[angles < 0] += 2 * np.pi #Angles need to be from 0 to 2pi

        points_polar = np.array([cartesian_points.T[2],
                                 np.sqrt(cartesian_points.T[0] ** 2 +
                                         cartesian_points.T[1] ** 2),
                                 angles]).T

        return points_polar

    @staticmethod
    def ScalarTripleProduct(A, B, C):

        """
        Returns (AxB) dot C where A, B, C are an array of 3-D vectors.
        This is 1/6th the volume of the tetrahedron spanned by ABC
        :param A: (n, 3) array of 3-D vectors
        :param B: (n, 3) array of 3-D vectors
        :param C: (n, 3) array of 3-D vectors
        :return: Array of (AxB) dot C for each entry in A, B, C
        """

        return np.sum(np.cross(A, B) * C, axis=1)