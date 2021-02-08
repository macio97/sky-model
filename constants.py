# Libraries
from math import pi
import numpy as np
from properties import look


# Constants
# Planck's constant (J*s)
PLANCK_CONSTANT = 6.62607004e-34
# Boltzmann constant (J*K^-1)
BOLTZMANN_CONSTANT = 1.38064852e-23
# molecular number density of air (m^-3)
MOLECULAR_NUMBER = 0.02504e27
# Dobson unit (molecules/m^2)
DOBSON_UNIT = 2.687e20
# speed of light (m/s)
LIGHT_SPEED = 299792458
# sun's Temperature (K)
SUN_TEMPERATURE = 5778
# IOR of air
AIR_IOR = 1.0002926
# Rayleigh scale height (m)
RAYLEIGH_SCALE = 8000
# Mie scale height (m)
MIE_SCALE = 1200
# aerosols anisotropy
MIE_G = 0.76
# squared mie_G
SQR_G = MIE_G * MIE_G
# average distance Earth-Sun (m)
DISTANCE_EARTH_SUN = 149.6e9
# radius of Sun (m)
SUN_RADIUS = 695500e3
# radius of Earth (m)
EARTH_RADIUS = 6360e3
# radius of atmosphere (m)
ATMOSPHERE_RADIUS = 6420e3
# number of sampled wavelengths
NUM_WAVELENGTHS = 21
# lowest sampled wavelength (nm)
MIN_WAVELENGTH = 380
# highest sampled wavelength (nm)
MAX_WAVELENGTH = 780
# step between each sampled wavelength (nm)
WAVELENGTHS_STEP = (MAX_WAVELENGTH - MIN_WAVELENGTH) / (NUM_WAVELENGTHS - 1)
# sampled wavelengths (m)
LAMBDA = np.arange(MIN_WAVELENGTH, MAX_WAVELENGTH +
                   1, WAVELENGTHS_STEP) * 10**-9
# blackbody spectral irradiance of sun (W*m^-2*nm^-1)
SUN_BLACKBODY = (2 * pi * PLANCK_CONSTANT * LIGHT_SPEED * LIGHT_SPEED) / (LAMBDA**5 *
                                                                          (np.exp((PLANCK_CONSTANT * LIGHT_SPEED) / (BOLTZMANN_CONSTANT * SUN_TEMPERATURE * LAMBDA)) - 1)) * 10**-9
# irradiance on top of atmosphere (W*m^-2*nm^-1)
SUN_IRRADIANCE = SUN_BLACKBODY * \
    ((SUN_RADIUS * SUN_RADIUS) / (DISTANCE_EARTH_SUN * DISTANCE_EARTH_SUN))
# Rayleigh scattering coefficient (m^-1)
RAYLEIGH_COEFFICIENT = ((8 * pi**3) * (AIR_IOR * AIR_IOR - 1)
                        ** 2) / (3 * MOLECULAR_NUMBER * LAMBDA**4)
# Mie scattering coefficient (m^-1)
MIE_COEFFICIENT = 2e-5
# maximum number density of ozone molecules (m^-3)
OZONE_MAX = 300 * DOBSON_UNIT / 15e3
# Ozone cross section (cm^2/molecules)
OZONE_CROSS = np.loadtxt('data/ozone_cross_section.csv', usecols=(1))
# Ozone absorption coefficient (m^-1)
OZONE_COEFFICIENT = OZONE_CROSS * 10**-4 * OZONE_MAX

# CIE XYZ color matching functions
COLOR_MATCHING_FUNCTIONS = np.loadtxt('data/cie_xyz.csv', usecols=(1, 2, 3))
# illuminants
ILLUMINANT_D65 = np.array([[3.2404542, -1.5371385, -0.4985314],
                           [-0.9692660, 1.8760108, 0.0415560],
                           [0.0556434, -0.2040259, 1.0572252]])
ILLUMINANT_E = np.array([[2.3706743, -0.9000405, -0.4706338],
                         [-0.5138850, 1.4253036, 0.0885814],
                         [0.0052982, -0.0146949, 1.0093968]])


def read_filmic_look(path):
    nums = []
    with open(path) as filmic_file:
        for line in filmic_file:
            nums.append(float(line))
    return nums


# filmic contrast
FILMIC_LOOK = read_filmic_look("looks/" + look + ".txt")
