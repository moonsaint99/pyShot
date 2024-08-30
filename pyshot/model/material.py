import numpy as np
import scipy as sp

class Material:
    def __init__(self, rho, vp, vs):
        self.rho = rho
        self.vp = vp
        self.vs = vs

        # Define elastic moduli
        self.shear_modulus = self.vs**2 * self.rho
        self.bulk_modulus = self.rho*self.vp**2 - 4/3*self.shear_modulus

    def __getitem__(self, item):
        return [self.rho, self.vp, self.vs][item]


def mat_from_moduli(bulk, shear, rho):
    """
    Compute the material properties from bulk and shear moduli
    :param bulk: bulk modulus
    :param shear: shear modulus
    :param rho: density
    :return: Material
    """
    vp = np.sqrt((bulk + 4/3*shear)/rho)
    vs = np.sqrt(shear/rho)
    return Material(rho, vp, vs)

def voigt_modulus(mat1mod, mat2mod, mat1frac):
    """
    Compute the Voigt average of two materials
    :param mat1mod:
    :param mat2mod:
    :param mat1frac:
    :return: voigt average
    """
    return mat1frac*mat1mod + (1-mat1frac)*mat2mod

def reuss_modulus(mat1mod, mat2mod, mat1frac):
    """
    Compute the Reuss average of two materials
    :param mat1mod:
    :param mat2mod:
    :param mat1frac: fraction of mat1
    :return: reuss average
    """
    return 1/(mat1frac/mat1mod + (1-mat1frac)/mat2mod)

def vrh_mod(mod1, mod2, frac):
    """
    Compute the Voigt-Reuss-Hill average of two materials
    :param mod1:
    :param mod2:
    :param frac: fraction of mod1
    :return: vrh average
    """
    voigt = voigt_modulus(mod1, mod2, frac)
    reuss = reuss_modulus(mod1, mod2, frac)
    return (voigt + reuss)/2

def vrh_mat(mat1, mat2, frac):
    """
    Compute the Voigt-Reuss-Hill average of two materials
    :param mat1:
    :param mat2:
    :param frac: fraction of mat1
    :return: vrh average
    """
    bulk = vrh_mod(mat1.bulk_modulus, mat2.bulk_modulus, frac)
    shear = vrh_mod(mat1.shear_modulus, mat2.shear_modulus, frac)
    rho = mat1.rho*frac + mat2.rho*(1-frac)
    return mat_from_moduli(bulk, shear, rho)