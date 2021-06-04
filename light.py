# Libraries
import functions as fun
import numpy as np
from constants import EARTH_RADIUS, SUN_IRRADIANCE, MIE_COEFFICIENT, NUM_WAVELENGTHS, OZONE_COEFFICIENT, RAYLEIGH_COEFFICIENT
from math import cos, exp, pi, radians, sin, sqrt, dist
from properties import air_density, altitude, dust_density, ozone_density, steps, sun_lat
from random import uniform as random


# Definitions
# convert altitude from km to m and clamp to avoid intersection issues
cam_altitude = 1000 * max(min(altitude, 59.999), 0.001)
camera_position = np.array([0, 0, EARTH_RADIUS + cam_altitude])
# convert sun latitude and longitude to vector
sun_direction = fun.geographical_to_direction(radians(sun_lat), 0)
# scattering and absorption coefficients
coefficients = np.array(
    [RAYLEIGH_COEFFICIENT, 1.11 * MIE_COEFFICIENT, OZONE_COEFFICIENT], dtype=object)
density_multipliers = np.array([air_density, dust_density, ozone_density])


def ray_optical_depth(ray_origin, ray_direction, ray_length):
    # step along the ray in segments and accumulate the optical depth along each segment
    segment_length = ray_length / steps
    segment = segment_length * ray_direction
    optical_depth = np.zeros(3)
    # the density of each segment is evaluated at its middle
    middle_point = ray_origin + 0.5 * segment

    for _ in range(steps):
        # height above sea level
        height = sqrt(np.dot(middle_point, middle_point)) - EARTH_RADIUS
        # accumulate optical depth of this segment
        density = np.array([fun.density_rayleigh(
            height), fun.density_mie(height), fun.density_ozone(height)])
        optical_depth += density
        # advance along ray
        middle_point += segment

    return optical_depth * segment_length


def multiple_scattering(ray_direction):
    '''
    still using 1 single sample for now,
    when ready we can use multiple samples per pixel as Monte Carlo method.
    the algorithm is (correct me if i'm wrong):
    - start ray from camera, end at random position, accumulate light from sun -> end point -> camera
    - then shoot a ray at random direction from previous end point, and accumulate sun -> new end point -> old end point -> camera
    - continue until the number of bounces end
    '''
    n_bounces = 4
    ray_origin = camera_position
    light = 0
    throughput = 1

    # light bounces per sample
    for i in range(n_bounces):
        # if Earth surface isn't hit then hit atmosphere (one or the other needs to be hit)
        distance = fun.surface_intersection(ray_origin, ray_direction)
        if distance < 0:
            distance = fun.atmosphere_intersection(ray_origin, ray_direction)
        # pick a random point between origin and end of ray
        ray_length = random(0, 1) * distance
        # send ray to sun, it has to end in atmosphere otherwise return black (?)
        ray_end = ray_origin + ray_direction * ray_length
        ray_length_sun = fun.surface_intersection(ray_end, sun_direction)
        if ray_length_sun < 0:
            ray_length_sun = fun.atmosphere_intersection(
                ray_end, sun_direction)
        else:
            break
        # calculate optical depths
        optical_depth = ray_optical_depth(
            ray_origin, ray_direction, ray_length)
        optical_depth_sun = ray_optical_depth(
            ray_end, sun_direction, ray_length_sun)
        # attenuation of light
        transmittance = np.exp(-np.sum(coefficients * optical_depth))
        transmittance_sun = np.exp(-np.sum(coefficients * optical_depth_sun))
        # phase function (sr^-1)
        mu = np.dot(ray_direction, sun_direction)
        phase_R = fun.phase_rayleigh(mu)
        phase_M = fun.phase_mie(mu)
        # densities on scattering point
        height = sqrt(np.dot(ray_end, ray_end)) - EARTH_RADIUS
        density_R = fun.density_rayleigh(height)
        density_M = fun.density_mie(height)
        # compute scattering
        scattering_R = coefficients[0] * density_R * phase_R
        scattering_M = coefficients[1] * density_M * phase_M
        scattering = scattering_R * phase_R + scattering_M * phase_M
        # accumulate light
        light += SUN_IRRADIANCE * throughput * scattering * transmittance_sun
        # change ray direction and origin for a new path
        sphere_lat = random(0, pi)
        sphere_lon = random(0, 2 * pi)
        new_ray_direction = fun.normalize_vector(sphere_lat, sphere_lon)
        ray_origin = ray_end

        if i != n_bounces-1:
            # if Earth surface isn't hit then hit atmosphere (one or the other needs to be hit)
            distance = fun.surface_intersection(ray_origin, new_ray_direction)
            if distance < 0:
                distance = fun.atmosphere_intersection(
                    ray_origin, new_ray_direction)
            # pick a random point between origin and end of ray
            ray_length = random(0, 1) * distance
            ray_end = ray_origin + ray_direction * ray_length
            # calculate optical depths
            optical_depth = ray_optical_depth(
                ray_origin, new_ray_direction, ray_length)
            # attenuation of light
            transmittance = np.exp(-np.sum(coefficients * optical_depth))
            # phase function (sr^-1)
            mu = np.dot(ray_direction, new_ray_direction)
            phase_R = fun.phase_rayleigh(mu)
            phase_M = fun.phase_mie(mu)
            # densities on scattering point
            height = sqrt(np.dot(ray_end, ray_end)) - EARTH_RADIUS
            density_R = fun.density_rayleigh(height)
            density_M = fun.density_mie(height)
            # compute scattering
            scattering_R = coefficients[0] * density_R * phase_R
            scattering_M = coefficients[1] * density_M * phase_M
            scattering = scattering_R * phase_R + scattering_M * phase_M

            # update throughput
            throughput *= transmittance * scattering

    return light


def monte_carlo(ray_direction):
    n_samples = 1
    srgb = np.array([0.0, 0.0, 0.0])

    # samples
    for _ in range(n_samples):
        # get spectrum from camera direction
        spectrum = multiple_scattering(ray_direction)
        # convert spectrum to xyz
        xyz = fun.spectrum_to_xyz(spectrum)
        # convert xyz to srgb
        srgb += fun.xyz_to_srgb(xyz)

    # average
    return srgb / n_samples
