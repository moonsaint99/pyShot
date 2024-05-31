import numpy as np
import scipy as sp

def ref_trans_array(mat1, mat2, theta1): #lowermat = [rho, vp, vs]
    theta2 = np.arcsin(mat2[1]/mat1[1] * np.sin(theta1))
    phi1 = np.arcsin(mat1[2]/mat1[1] * np.sin(theta1))
    phi2 = np.arcsin(mat2[2]/mat1[1] * np.sin(theta1))
    matrix = np.array([[-np.sin(theta1), -np.cos(phi1), np.sin(theta2), np.cos(phi2)],
                       [np.cos(theta1), -np.sin(phi1), np.cos(theta2), -np.sin(phi2)],
                       [np.sin(2*theta1), mat1[1]/mat1[2]*np.cos(2*phi1),
                        mat2[0]*mat2[2]**2*mat1[1]/(mat1[0]*mat1[2]**2*mat2[1])*np.sin(2*theta2),
                        mat2[0]*mat2[2]*mat1[1]/(mat1[0]*mat1[2]**2)*np.cos(2*phi2)],
                       [-np.cos(2*phi1), mat1[2]/mat1[1]*np.sin(2*phi1),
                        mat2[0]*mat2[1]/(mat1[0]*mat1[1])*np.cos(2*phi2),
                        -mat2[0]*mat2[2]/(mat1[0]*mat1[1])*np.sin(2*phi2)]])
    matrix = np.linalg.inv(matrix)
    array = np.array([np.sin(theta1), np.cos(theta1), np.sin(2*theta1), np.cos(2*phi1)])
    return np.matmul(matrix, array)[0] # remove the [0] to get the full reflectivity vector, not just P wave


def refl_optimize_kernel(lowermat, uppermat, theta, minmax):
    ice_rho = 1187
    ice_vp = 3840
    ice_vs = 1930
    return minmax*ref_trans_array(uppermat, lowermat, theta)


def reflectivity_optimizer(mat1min, mat1max, uppermat, theta, optimum):
    mat1_rhomin, mat1_vpmin, mat1_vsmin = mat1min
    mat1_rhomax, mat1_vpmax, mat1_vsmax = mat1max
    if optimum == 'min':
        minmax = 1
    elif optimum == 'max':
        minmax = -1
    else:
        raise ValueError('Optimum must be either "min" or "max"')
    optimum = sp.optimize.minimize(refl_optimize_kernel, np.array([mat1_rhomin, mat1_vpmin, mat1_vsmin]), args=(uppermat, theta, minmax),
                                   bounds=[(mat1_rhomin, mat1_rhomax), (mat1_vpmin, mat1_vpmax), (mat1_vsmin, mat1_vsmax)])

    return optimum.x


def reflectivity_window(mat1min, mat1max, uppermatmin, uppermatmax, theta_array):
    max_refl = []
    min_refl = []
    for theta in theta_array:
        max_params = reflectivity_optimizer(mat1min, mat1max, uppermatmin, theta, 'max')
        max_refl.append(-1*refl_optimize_kernel(max_params, uppermatmin, theta, -1))
        min_params = reflectivity_optimizer(mat1min, mat1max, uppermatmax, theta, 'min')
        min_refl.append(refl_optimize_kernel(min_params, uppermatmax, theta, 1))
    return np.array(max_refl), np.array(min_refl)


def reflectivity_array(lowermat, uppermat, theta_array):
    """
    Compute the reflectivity for an array of angles assuming
    :param lowermat: [rho, vp, vs] for the lower layer
    :param uppermat: [rho, vp, vs] for the upper layer
    :param theta_array:
    :return: Reflectivity Array
    """
    reflectivity_array = []
    for theta in theta_array:
        reflectivity_array.append(refl_optimize_kernel(lowermat, uppermat, theta, 1))
    return np.array(reflectivity_array)


def shuey_reflectivity(mat1, mat2, theta):
    """
    Compute the Shuey reflectivity for a single angle
    :param mat1: [rho, vp, vs] for the upper layer
    :param mat2: [rho, vp, vs] for the lower layer
    :param theta: angle of incidence in radians
    :return: Shuey reflectivity
    """
    # Compute the P-wave reflectivity
    delta_vp = mat2[1] - mat1[1]
    delta_vs = mat2[2] - mat1[2]
    delta_rho = mat2[0] - mat1[0]
    r0 = 0.5 * (delta_vp / mat1[1] + delta_rho / mat1[0])
    g = 0.5 * delta_vp / mat1[1] - 2 * mat1[2] ** 2 / mat1[1] ** 2 * (delta_rho / mat1[0] + 2 * delta_vs / mat1[2])
    rtheta = r0 + g * np.sin(theta) ** 2
    return rtheta

def shuey_kernel(lowermat, uppermat, theta, minmax):
    ice_rho = 917
    ice_vp = 3840
    ice_vs = 1930
    return minmax*shuey_reflectivity([ice_rho, ice_vp, ice_vs], lowermat, theta)


def shuey_array(lowermat, uppermat, theta_array):
    """
    Compute the Shuey reflectivity for an array of angles assuming
    a single lower layer below ice
    :param lowermat: [rho, vp, vs] for the lower layer
    :param theta_array: array of angles of incidence in radians
    :return: Shuey reflectivity
    """
    # Compute the P-wave reflectivity
    reflectivity_array = []
    for theta in theta_array:
        reflectivity_array.append(shuey_kernel(lowermat, uppermat, theta, 1))
    return np.array(reflectivity_array)


def normal_reflectivity(mat1, mat2, theta_array):
    """
    Compute the normal incidence reflectivity
    :param mat1: [rho, vp, vs] for the upper layer
    :param mat2: [rho, vp, vs] for the lower layer
    :return: Normal incidence reflectivity
    """
    r = (mat2[0] * mat2[1] - mat1[0] * mat1[1]) / (mat2[0] * mat2[1] + mat1[0] * mat1[1])
    reflectivity_array = []
    for theta in theta_array:
        reflectivity_array.append(r)
    return np.array(reflectivity_array)