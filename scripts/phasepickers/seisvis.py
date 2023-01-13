#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 12:42:04 2020

@author: spri902
"""
import warnings
import numpy
import matplotlib.pyplot as plt

def seismic_wiggle(section, dt, picktimes=None, ranges=None, scale=1., color=None,
                   globnormalize=False,trnormalize=False):
    """
    Plot a seismic section (numpy 2D array matrix) as wiggles.
    Parameters:
    * picktimes:
        phase arrival picktimes in seconds
    * section :  2D array
        matrix of traces (first dimension time, second dimension traces)
    * dt : float
        sample rate in seconds
    * ranges : (x1, x2)
        min and max horizontal values (default trace number)
    * scale : float
        scale factor multiplied by the section values before plotting
    * color : tuple of strings
        Color for marking picktimes.
    * normalize :
        True to normalizes all trace in the section using global max/min
        data will be in the range (-0.5, 0.5) zero centered
    .. warning::
        Slow for more than 200 traces, in this case decimate your
        data or use ``seismic_image``.
    """
    npts, ntraces = section.shape  # time/traces
    if ntraces < 1:
        raise IndexError("Nothing to plot")
    if npts < 1:
        raise IndexError("Nothing to plot")
    t = numpy.linspace(0, dt*npts, npts)
    amp = 1.  # normalization factor
    gmin = 0.  # global minimum
    toffset = 0.  # offset in time to make 0 centered
    if globnormalize:
        gmax = section.max()
        gmin = section.min()
        amp = (gmax - gmin)
        toffset = 0.5
    if trnormalize:
        gmax = numpy.max(numpy.abs(section),0)
    #plt.ylim(max(t), 0)
    plt.ylim(t[400],0)
    if ranges is None:
        ranges = (0, ntraces)
    x0, x1 = ranges
    # horizontal increment
    dx = (x1 - x0)/ntraces
    #sc=numpy.abs([section.max(), section.min()]).max()
    plt.xlim(x0-1,x1)
    if globnormalize:
        if picktimes is None:
            for i, trace in enumerate(section.transpose()):
                tr = (((trace - gmin)/amp) - toffset)*scale*dx
                x = x0 + i*dx  # x positon for this trace
                plt.plot(x + tr, t, 'k')
        else:
            for i, trace in enumerate(section.transpose()):
                tr = (((trace - gmin)/amp) - toffset)*scale*dx
                x = x0 + i*dx  # x positon for this trace
                plt.plot(x + tr, t, 'k')
                #pyplot.fill_betweenx(t, x + tr, x, tr > 0, color=color)
                plt.plot(x,picktimes[i][0],'.',color=color)
    
    elif trnormalize:
        if picktimes is None:
            for i, trace in enumerate(section.transpose()):
                tr = ((trace/gmax[i]) - toffset)*scale*dx
                x = x0 + i*dx  # x positon for this trace
                plt.plot(x + tr, t, 'k')
                #pyplot.fill_betweenx(t, x + tr, x, tr > 0, color=color)
        else:
            for i, trace in enumerate(section.transpose()):
                tr = ((trace/gmax[i]) - toffset)*scale*dx
                x = x0 + i*dx  # x positon for this trace
                plt.plot(x + tr, t, 'k')
                #pyplot.fill_betweenx(t, x + tr, x, tr > 0, color=color)
                plt.plot(x,picktimes[i][0],'.',color=color)
    else:
        for i, trace in enumerate(section.transpose()):
            tr = trace*scale*dx
            x = x0 + i*dx  # x positon for this trace
            plt.plot(x + tr, t,'k')
        
def seismic_image(section, dt, ranges=None, cmap=plt.cm.gray,
                  aspect=None, vmin=None, vmax=None):
    """
    Plot a seismic section (numpy 2D array matrix) as an image.
    Parameters:
    * section :  2D array
        matrix of traces (first dimension time, second dimension traces)
    * dt : float
        sample rate in seconds
    * ranges : (x1, x2)
        min and max horizontal coordinate values (default trace number)
    * cmap : colormap
        color map to be used. (see pyplot.cm module)
    * aspect : float
        matplotlib imshow aspect parameter, ratio between axes
    * vmin, vmax : float
        min and max values for imshow
    """
    npts, maxtraces = section.shape  # time/traces
    if maxtraces < 1:
        raise IndexError("Nothing to plot")
    if npts < 1:
        raise IndexError("Nothing to plot")
    t = numpy.linspace(0, dt*npts, npts)
    data = section
    if ranges is None:
        ranges = (0, maxtraces)
    x0, x1 = ranges
    extent = (x0, x1, t[-1], t[0])
    if aspect is None:  # guarantee a rectangular picture
        aspect = numpy.round((x1 - x0)/numpy.max(t))
        aspect -= aspect*0.2
    if vmin is None and vmax is None:
        scale = numpy.abs([section.max(), section.min()]).max()
        vmin = -scale
        vmax = scale
    plt.imshow(data, aspect=aspect, cmap=cmap, origin='upper',
                  extent=extent, vmin=vmin, vmax=vmax)
    
def shot_gather(section, dt, numSamp=1000,ranges=None, scale=1., color='k',
                   globnormalize=False,trnormalize=False):
    """
    Plot a seismic section (numpy 2D array matrix) as wiggles.
    Parameters:

    * section :  2D array
        matrix of traces (first dimension time, second dimension traces)
    * dt : float
        sample rate in seconds
    * ranges : (x1, x2)
        min and max horizontal values (default trace number)
    * scale : float
        scale factor multiplied by the section values before plotting
    * color : tuple of strings
        Color for filling the wiggle, positive  and negative lobes.
    * normalize :
        True to normalizes all trace in the section using global max/min
        data will be in the range (-0.5, 0.5) zero centered
    .. warning::
        Slow for more than 200 traces, in this case decimate your
        data or use ``seismic_image``.
    """
    npts, ntraces = section.shape  # time/traces
    if ntraces < 1:
        raise IndexError("Nothing to plot")
    if npts < 1:
        raise IndexError("Nothing to plot")
    t = numpy.linspace(0, dt*npts, npts)
    amp = 1.  # normalization factor
    gmin = 0.  # global minimum
    toffset = 0.  # offset in time to make 0 centered
    if globnormalize:
        gmax = section.max()
        gmin = section.min()
        amp = (gmax - gmin)
        toffset = 0.5
    if trnormalize:
        gmax = numpy.max(numpy.abs(section),0)
    plt.ylim(numSamp*dt, 0)
    if ranges is None:
        ranges = (0, ntraces)
    x0, x1 = ranges
    # horizontal increment
    dx = (x1 - x0)/ntraces
    plt.xlim(x0-1, x1+1)
    if globnormalize:
        for i, trace in enumerate(section.transpose()):
            tr = (((trace - gmin)/amp) - toffset)*scale*dx
            x = x0 + i*dx  # x positon for this trace
            plt.plot(x + tr, t, color=color)
    if trnormalize:
        for i, trace in enumerate(section.transpose()):
            tr = ((trace/gmax[i]) - toffset)*scale*dx
            x = x0 + i*dx  # x positon for this trace
            plt.plot(x + tr, t, color=color)
def plot_section(section, dt, ranges=None):
    npts, ntraces = section.shape  # time/traces
    t = numpy.linspace(0, dt*npts, npts)
    plt.ylim(max(t), 0)
    if ranges is None:
        ranges = (0, ntraces)
    x0, x1 = ranges
    # horizontal increment
    dx = (x1 - x0)/ntraces
    #sc=numpy.abs([section.max(), section.min()]).max()
    plt.xlim(x0,x1)
    for i, trace in enumerate(section.transpose()):
        tr = trace*dx
        x = x0 + i*dx  # x positon for this trace
        plt.plot(x + tr, t,'k')