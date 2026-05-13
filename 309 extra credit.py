#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import plotly.offline as ploff
import plotly.graph_objects as go
from scipy import constants

import plotly.io as pio
ploff.init_notebook_mode()
# pio.renderers.default = "notebook_connected"

a0 = constants.physical_constants['Bohr radius'][0]
Angstrom = constants.angstrom


# In[3]:


def r_phi_theta(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x) + np.pi
    theta = np.arccos(z/r)
    return r, phi, theta


# In[4]:


X, Y, Z = np.mgrid[-20*a0:20*a0:48j, -20*a0:20*a0:48j, -20*a0:20*a0:48j]
r_grid, phi_grid, theta_grid = r_phi_theta(X, Y, Z)


# In[5]:


def probfunc(r, wavefunc):
    p = r**2 * np.conj(wavefunc) * wavefunc
    return p.real / np.max(p.real)
# The .real syntax picks off the real part of a complex number. Even though the 
# imaginary part of the probability is zero, variable p is still represented in
# the computer as a complex number. Our plotting routines will only work on
# arrays of real numbers.


# In[6]:


n1l0m0 = lambda r: np.exp(-r/a0) / (np.sqrt(np.pi) * a0**1.5)

n2l0m0 = lambda r: (2 - r/a0) * np.exp(-r/(2*a0)) / (4 * np.sqrt(2*np.pi) * a0**1.5)

#Only depends on r because there is no angular behavior, l = m = 0.


# In[7]:


prob_n1l0m0 = probfunc(r_grid, n1l0m0(r_grid))
prob_n2l0m0 = probfunc(r_grid, n2l0m0(r_grid))


# In[9]:


def visualize_prob(probability_function):

    fig = go.Figure(data=go.Isosurface(
        x=X.flatten() / Angstrom, # x values of the grid in Angstroms
        y=Y.flatten() / Angstrom, # y values of the grid in Angstroms
        z=Z.flatten() / Angstrom, # z values of the grid in Angstroms
        value=probability_function.flatten(), # independent variable
        isomin=0.05, # Minimum normalized probability density to render in an isosurface
        isomax=0.95, # Maximum normalized probability density to render in an isosurface
        opacity=0.4, # Set a low opacity so each surface is partially transparent
        colorscale='Plotly3_r', # Nice-looking color table
        surface_count=8, # number of isosurfaces to plot (2 by default: only min and max)
        colorbar_nticks=8, # colorbar ticks correspond to isosurface values
        caps=dict(x_show=False, y_show=False, z_show=False))
        )

    # Change axis lables and make the plot larger than the default
    fig.update_layout(scene = dict(
                        xaxis_title='x (Angstroms)',
                        yaxis_title='y (Angstroms)',
                        zaxis_title='z (Angstroms)'),
                        width=700,
                        margin=dict(r=10, b=10, l=10, t=10))

    fig.show()


# In[14]:


visualize_prob(prob_n1l0m0) #higher probability of electron being close to the proton, thus more centrally condensed.
visualize_prob(prob_n2l0m0)


# In[12]:


def visualize_prob_slice(probability_function):

    fig = go.Figure(data=go.Isosurface(
        x=X.flatten() / Angstrom, # x values of the grid in Angstroms
        y=Y.flatten() / Angstrom, # y values of the grid in Angstroms
        z=Z.flatten() / Angstrom, # z values of the grid in Angstroms
        value=probability_function.flatten(), # independent variable
        isomin=0.05, # Minimum normalized probability density to render in an isosurface
        isomax=0.95, # Maximum normalized probability density to render in an isosurface
        opacity=0.4, # Set a low opacity so each surface is partially transparent
        colorscale='Plotly3_r', # Nice-looking color table
        surface_count=8, # number of isosurfaces to plot (2 by default: only min and max)
        colorbar_nticks=8, # colorbar ticks correspond to isosurface values
        slices_x=dict(show=True, locations=[0]),
        caps=dict(x_show=False, y_show=False, z_show=False))
        )

    # Change axis lables and make the plot larger than the default
    fig.update_layout(scene = dict(
                        xaxis_title='x (Angstroms)',
                        yaxis_title='y (Angstroms)',
                        zaxis_title='z (Angstroms)'),
                        width=700,
                        margin=dict(r=10, b=10, l=10, t=10))

    fig.show()


# In[13]:


visualize_prob_slice(prob_n1l0m0) #matches expectation.
visualize_prob_slice(prob_n2l0m0)


# In[ ]:


#Based on the visualizations, we can safely say that the two visualizations are fairly different. 
#when n increases, the probability of the electron being at further distances away from the proton increases
#i.e. r triples bewtween n = 1 and n = 2.


# In[15]:


n2l1m0 = lambda r, theta: (r/a0) * np.exp(-r/(2*a0)) * np.cos(theta) \
                          / (4 * np.sqrt(2*np.pi) * a0**1.5)

visualize_prob_slice(probfunc(r_grid, n2l1m0(r_grid, theta_grid)))


# In[16]:


n2l1m1 = lambda r, phi, theta: (r/a0) * np.exp(-r/(2*a0)) * np.sin(theta) \
                               * np.exp(1j * phi) / (8 * np.sqrt(3*np.pi) * a0**1.5)
n3l2m1 = lambda r, phi, theta: (r**2/a0**2) * np.exp(-r/(3*a0)) * np.sin(theta) * \
                               np.cos(theta) * (np.cos(phi) + 1j * np.sin(phi)) / \
                               (81 * np.sqrt(np.pi) * a0**1.5)
visualize_prob_slice(probfunc(r_grid, n3l2m1(r_grid, phi_grid, theta_grid)))


# In[18]:


n3l4m2 = lambda r, phi, theta: ((r/a0)**4 * np.exp(-r/(5*a0)) * (np.sin(theta)**2) * (7*np.cos(theta)**2 - 1)* np.exp(2j * phi)
    / (18750 * np.sqrt(7*np.pi) * a0**1.5))

visualize_prob_slice(probfunc(r_grid, n5l4m2(r_grid, phi_grid, theta_grid)))


# In[ ]:




