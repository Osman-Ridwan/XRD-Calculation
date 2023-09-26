#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import matplotlib.pyplot as plt


# # lattice parameters a1, a2, and a3

# In[13]:


direct_lattice_a1 = np.array([5.599559, 0.0, 0.0])
direct_lattice_a2 = np.array([0.0, 5.599559, 0.0])
direct_lattice_a3 = np.array([0.0, 0.0, 5.599559])
wavelength = 1.54184


# # maximum Miller indices

# In[11]:


h_max = 4
k_max = 4
l_max = 4


# # fractional atomic coordinates for NaCl

# In[14]:


fractional_coords = np.array([
    [0.0, 0.0, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.5],
    [0.0, 0.0, 0.5],
    [0.0, 0.5, 0.0]
])


# In[15]:


atomic_coefficients_na = (4.763, 3.174, 1.267, 1.113, 3.285, 8.842, 0.314, 129.424, 0.676)


# In[16]:


two_theta_values = []
intensity_values = []


# # compute b1, b2, b3

# In[18]:


def reciprocal_lattice_vectors(a1, a2, a3):
    b1 = np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
    b2 = np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
    b3 = np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))
    return b1, b2, b3


# # compute Dhkl

# In[20]:


def compute_dhkl(h, k, l, b1, b2, b3):
    return np.linalg.norm(h * b1 + k * b2 + l * b3)


# # compute atomic scattering factors

# In[21]:


def atomic_scattering_factor(sin_over_lambda, coefficients):
    a1, a2, a3, a4, b1, b2, b3, b4, c = coefficients
    return (a1 * np.exp(-b1 * (sin_over_lambda**2)) + 
            a2 * np.exp(-b2 * (sin_over_lambda**2)) + 
            a3 * np.exp(-b3 * (sin_over_lambda**2)) + 
            a4 * np.exp(-b4 * (sin_over_lambda**2)) + c)


# # Loop over h,k,l

# In[23]:


for h in range(h_max + 1):
    for k in range(k_max + 1):
        for l in range(l_max + 1):
            if h == 0 and k == 0 and l == 0:
                continue  

            # Compute b1, b2, b3
            reciprocal_vectors_b1, reciprocal_vectors_b2, reciprocal_vectors_b3 = reciprocal_lattice_vectors(direct_lattice_a1, direct_lattice_a2, direct_lattice_a3)

            #Compute Dhkl
            dhkl = compute_dhkl(h, k, l, reciprocal_vectors_b1, reciprocal_vectors_b2, reciprocal_vectors_b3)

            #Check the condition
            if 1 / dhkl < 2 / wavelength:

            #alculate 2 Theta
                sin_over_lambda = wavelength / (2 * dhkl)
                two_theta = 2 * np.arcsin(sin_over_lambda)

            #Compute the structure factor
                structure_factor_Fhkl = 0
                for j in range(len(fractional_coords)):
                    fractional_xj, fractional_yj, fractional_zj = fractional_coords[j]
                    phase = 2 * np.pi * (h * fractional_xj + k * fractional_yj + l * fractional_zj)
                    fj = atomic_scattering_factor(sin_over_lambda, atomic_coefficients_na)
                    structure_factor_Fhkl += fj * np.exp(2j * phase)

            #Calculate intensity and store values
                intensity = np.abs(structure_factor_Fhkl)**2
                two_theta_values.append(np.degrees(two_theta))
                intensity_values.append(intensity)


# # Plot the XRD pattern

# In[27]:



plt.figure(figsize=(10, 6))
plt.bar(two_theta_values, intensity_values, width=0.1, color='blue', label='XRD Pattern (FCC-NaCl)')
plt.xlabel('2θ (degrees)')
plt.ylabel('Intensity')
plt.title('X-ray Diffraction Pattern for FCC-NaCl')
plt.legend()
plt.grid()
plt.show()

