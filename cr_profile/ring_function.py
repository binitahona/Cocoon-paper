__author__ = 'henrike'

import exceptions
import math

import astropy.units as u
import numpy as np
import warnings

from astromodels.core.units import get_units
from astromodels.utils.angular_distance import angular_distance
from astromodels.functions.function import Function2D, FunctionMeta, ModelAssertionViolation
import warnings



class NaimaNotAvailable(ImportWarning):
    pass


class InvalidUsageForFunction(exceptions.Exception):
    pass


class Ring_on_sphere(Function2D):
    r"""
        description :

            A bidimensional ring/tophat function on a sphere (in spherical coordinates)


        parameters :

            lon0 :

                desc : Longitude of the center of the source
                initial value : 0.0
                min : 0.0
                max : 360.0

            lat0 :

                desc : Latitude of the center of the source
                initial value : 0.0
                min : -90.0
                max : 90.0

            radius_inner :

                desc : Inner radius of the ring
                initial value : 0.5
                min : 0
                max : 20
                
            radius_outer :

                desc : Outer radius of the ring
                initial value : 1.0
                min : 0
                max : 20


        """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.radius_inner.unit = x_unit
        self.radius_outer.unit = x_unit

    def evaluate(self, x, y, lon0, lat0, radius_inner, radius_outer):

        lon, lat = x,y

        angsep = angular_distance(lon0, lat0, lon, lat)

        return np.power(180 / np.pi, 2) * 1. / (np.pi * (radius_outer ** 2 - radius_inner ** 2 )) * (angsep <= radius_outer) * (angsep > radius_inner)

    def get_boundaries(self):

        # Truncate the disk at 2 times the max of radius allowed

        max_radius = self.radius_outer.max_value

        min_lat = max(-90., self.lat0.value - 2 * max_radius)
        max_lat = min(90., self.lat0.value + 2 * max_radius)

        max_abs_lat = max(np.absolute(min_lat), np.absolute(max_lat))

        if max_abs_lat > 89. or 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

            min_lon = 0.
            max_lon = 360.

        else:

            min_lon = self.lon0.value - 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.)
            max_lon = self.lon0.value + 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.)

            if min_lon < 0.:

                min_lon = min_lon + 360.

            elif max_lon > 360.:

                max_lon = max_lon - 360.

        return (min_lon, max_lon), (min_lat, max_lat)

