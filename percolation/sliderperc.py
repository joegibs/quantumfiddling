# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:10:56 2023

@author: jogib
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import ndimage


L = 50
# pv = [0.2,0.3,0.4,0.5,0.6,0.7]
# z = rand(L,L)
# for i in range(len(pv)):
#     z = rand(L,L)
#     p = pv[i]
#     m = z<p
#     lw, num = ndimage.label(m)
#     area = ndimage.sum(m, lw, index=arange(lw.max() + 1))
#     areaImg = area[lw]
#     subplot(2,3,i+1)
#     tit = 'p='+str(p)
#     imshow(areaImg, origin='lower',cmap='gist_heat_r')
#     title(tit)
#     axis()
    
# The parametrized function to be plotted
def f(p):
    z = np.random.rand(L,L*2)
    m = z<p
    lw, num = ndimage.label(m)
    b = arange(lw.max() + 1)
    shuffle(b[1:])
    shuffledLw = b[lw]
    # area = ndimage.sum(m, lw, index=arange(lw.max() + 1))
    # areaImg = area[lw]
    return shuffledLw


# Define initial parameters
init_p = 0.2
init_frequency = 3

# Make a horizontal slider to control the frequency.
axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
freq_slider = Slider(
    ax=axfreq,
    label='Probability',
    valmin=0.0,
    valmax=1,
    valinit=init_frequency,
)

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
im1 = ax.imshow(f(freq_slider.val),cmap='gist_heat_r',origin='lower')

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control the frequency.
axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
freq_slider = Slider(
    ax=axfreq,
    label='Probability',
    valmin=0.0,
    valmax=1,
    valinit=init_frequency,
)

# Make a vertically oriented slider to control the amplitude
# axamp = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
# amp_slider = Slider(
#     ax=axamp,
#     label="Amplitude",
#     valmin=0,
#     valmax=10,
#     valinit=init_amplitude,
#     orientation="vertical"
# )


# The function to be called anytime a slider's value changes
def update(val):
    im1 = ax.imshow(f(freq_slider.val),cmap='gist_heat_r',origin='lower')
    fig.canvas.draw_idle()


# register the update function with each slider
freq_slider.on_changed(update)

# # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
# resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
# button = Button(resetax, 'Reset', hovercolor='0.975')


# def reset(event):
#     freq_slider.reset()
# button.on_clicked(reset)

plt.show()