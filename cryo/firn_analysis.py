import numpy as np
import scipy as sp

def double_linear_exponential(x, a, b, c, d, f):
    return a*(1-np.exp(-b*x)) + c*(1-np.exp(-d*x)) + f*x


def double_linear_exponential_fit(args): # returns fit of data to exponential function inv1
    # args is an iterable of direct arrival pickfile objects
    dist_array = np.array([])
    tmin_array = np.array([])
    for i in args:
        dist_array = np.append(dist_array, i.dist)
        tmin_array = np.append(tmin_array, i.tmin)
    return sp.optimize.curve_fit(
        f=double_linear_exponential,
        xdata=dist_array,
        ydata=tmin_array,
        #p0=np.array([0.85/100,0.035,1/1000,1.4, 1/3900]),
        p0=([0.011,0.04, 0.35, 0.005, 1/3630]),
        bounds=([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf, np.inf]),
        maxfev=10000000
    )


def double_linear_exponential_slope(x, params):  # This will accept a tuple of parameters
    # returns slowness u=1/dtdx based on the derivative of the exponential distance-time function inv1
    a, b, c, d, f = params
    dtdx = a*b*np.exp(-b*x) + c*d*np.exp(-d*x) + f
    return 1/dtdx


def firn_depth_vs_velocity(dist, params):
    '''
    This function calculates the maximum depth reached by rays traveling a given source-receiver offset
    It then gives the velocity at that depth
    :param dist:
    :param params:
    :return:
    '''
    # dist is the source-receiver offset for which you want to calculate
    # the maximum depth of the raypath, and the velocity at that maximum depth
    # If dist is an array, this will output an array of depths and velocities
    # params is the output from double_linear_exponential_fit
    vel_grad = double_linear_exponential_slope(dist, params)
    # vel_apparent = primary.dist_no_outliers/primary.tmin_no_outliers
    vel_apparent = dist / double_linear_exponential(dist, *params)
    z = np.array([])
    for i in range(len(dist)):
        z_int = sp.integrate.quad(lambda x: 1/(np.arccosh(vel_grad[i]/vel_apparent[i])), 0, dist[i])
        z_temp = 1/np.pi * z_int[0]
        z = np.append(z, z_temp)
    return z, vel_grad


def inc_angle(primary, params):
    vmax = double_linear_exponential_slope(primary.dist, params)
    vmin = double_linear_exponential_slope(0, params)
    return np.arcsin(vmin/vmax)