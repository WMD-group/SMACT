#!/usr/bin/env python
"""
Draws Hinton diagrams using matplotlib ( http://matplotlib.sf.net/ ).
Hinton diagrams are a handy way of visualizing weight matrices, using
colour to denote sign and area to denote magnitude.

By David Warde-Farley -- user AT cs dot toronto dot edu (user = dwf)
  with thanks to Geoffrey Hinton for providing the MATLAB code off of
  which this is modeled.

Redistributable under the terms of the 3-clause BSD license
(see http://www.opensource.org/licenses/bsd-license.php for details)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import smact.core as core
import pylab

def convert_to_hex(rgba_color) :
    '''Conversion of a colour in rgba format to hex'''
    red = int(rgba_color[0]*255)
    green = int(rgba_color[1]*255)
    blue = int(rgba_color[2]*255)
    return '#{r:02x}{g:02x}{b:02x}'.format(r=red,g=green,b=blue)

def _blob(x, y, area, colour):
    """
    Draws a square-shaped blob with the given area (< 1) at
    the given coordinates.
    """
    hs = np.sqrt(area) / 2
    xcorners = np.array([x - hs, x + hs, x + hs, x - hs])
    ycorners = np.array([y - hs, y - hs, y + hs, y + hs])
    plt.fill(xcorners, ycorners, colour, edgecolor='black')

def hinton_2D(W,V,R,S, maxweight=None):
    """
    W defines the colour of the squares, V defines the size
    Hinton diagram with 4 sets of data
    """
    search_space = core.ordered_elements(1,len(W))
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    jet = cm = plt.get_cmap('jet')
    cool = cm = plt.get_cmap('PuOr')
    minw = np.min(W)
    maxw = np.max(W)
    minr = np.min(R)
    maxr = np.max(R)
    cNorm  = colors.Normalize(vmin=minw, vmax=maxw)
    cNorm_r  = colors.Normalize(vmin=minr, vmax=maxr)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    scalarMap_r = cmx.ScalarMappable(norm=cNorm_r, cmap=cool)

    #cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=jet,norm=cNorm,orientation='horizontal')
    #plt.savefig('Cb1.pdf')
    #plt.show()
    cb2 = mpl.colorbar.ColorbarBase(ax1, cmap=cool,norm=cNorm_r,orientation='horizontal')
    plt.savefig('Cb2.pdf')
    plt.show()


    reenable = True
    if plt.isinteractive():
        plt.ioff()
    plt.clf()
    height, width = V.shape
    if not maxweight:
        maxweight = 2**np.ceil(np.log(np.max(np.abs(V)))/np.log(2))
        maxweight_s = 2**np.ceil(np.log(np.max(np.abs(S)))/np.log(2))
    else:
	maxweight_s = maxweight

    x1 = np.arange(0,width+1)
    y  = -x1 + height
    plt.fill_between(x1,height, y,
            facecolor='#7A7A87')
    plt.fill_between(x1,0, y,
            facecolor='#A19798')

    plt.axis('off')
    #plt.axis('equal')
    print width
    for x in range(0,width):
        for y in xrange(x):
            _x = x+1
            _y = y+1
            w = W[y,x]
            v = V[y,x]
	    if v > 0:
		weight = 0.75
	    else:
		weight = 0
            colorVal = scalarMap.to_rgba(w)
            colorVal_hex = convert_to_hex(colorVal)
            _blob(_x - 0.5,
                  height - _y + 0.5,
		  weight,
                  #min(1, maxweight),
                  colorVal_hex)

    for y in xrange(height):
	for x in xrange(y):
            _x = x+1
            _y = y+1
            r = R[y, x]
            s = S[y, x]
	    if s > 0:
		weight = 0.75
	    else:
		weight = 0
            colorVal = scalarMap_r.to_rgba(r)
            colorVal_hex = convert_to_hex(colorVal)
            _blob(_x - 0.5,
                  height - _y + 0.5,
		  weight,
                  #min(1, maxweight_s),   Turn this on if you want size to scale to a certain variable
                  colorVal_hex)

    if reenable:
        plt.ion()

def hinton(W,V, maxweight=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    Temporarily disables matplotlib interactive mode if it is on,
    otherwise this takes forever.

    W defines the colour of the squares, V defines the size
    """
    print maxweight
    jet = cm = plt.get_cmap('jet')
    minx = np.min(W)
    #miny = np.min(W[1])
    #min_val=min(minx,miny)
    min_val=minx
    maxx = np.max(W)
    #maxy = np.max(W[1])
    #max_val=max(maxx,maxy)
    max_val=maxx

    fig = plt.figure(figsize=(12,12))
    ax = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    cNorm  = colors.Normalize(vmin=0, vmax=4)
    cNorm_r  = colors.Normalize(vmin=0, vmax=4)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    reenable = True
    if plt.isinteractive():
        plt.ioff()

    cb2 = mpl.colorbar.ColorbarBase(ax, cmap=jet,norm=cNorm_r,orientation='horizontal')
    plt.savefig('colourbar.pdf')
    plt.clf()
    height, width = V.shape
    if not maxweight:
        maxweight = np.max(V)
        maxweight = 2**np.ceil(np.log(np.max(np.abs(V)))/np.log(2))
    print 'Maxweight', maxweight
    plt.fill(np.array([0, width, width, 0]),
             np.array([0, 0, height, height]),
             'gray')

    #plt.axis('off')
    plt.axis('equal')
    for x in xrange(width):
        for y in xrange(height):
            _x = x+1
            _y = y+1
            w = W[y, x]
	    v = V[y, x]
	    colorVal = scalarMap.to_rgba(w)
	    colorVal_hex = convert_to_hex(colorVal)
            _blob(_x - 0.5,
                  height - _y + 0.5,
                  1 - min(1, v/maxweight),
                  colorVal_hex)
    if reenable:
        plt.ion()

if __name__ == "__main__":
    hinton(np.random.randn(20, 20),np.random.randn(20,20))
    plt.title('Example Hinton diagram - 20x20 random normal')
    plt.show()
