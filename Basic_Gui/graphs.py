#-----------------------------------------------------------------------
"""
Display program for PIC Gui
This is the graphics (View) component of the PIC GUI.
It displays multiple plots requested by the physics (model) component.
There can be up to 4 plots per window and mutiple windows.  The
dictionary comms.plot_loc contains the plot names and identifiers.
The dictionary comms.sim_data contains additional informationb that
may be needed by the plots.

The physics component communicates with the Control Panel by generating
a custom event <<PlotChange>> to request a plot.  The Control Panel then
calls a plot function in the graphics component.  The graphics component 
communicates with the physics component by writing to the queue plotq to
synchronize the plots.

Main graphic procedures are:
dscaler1   displays 1d scalar field in real space.
displayfv1 displays 1d velocity distribution functions
grasp1     displays 1d (iyp-ixp) phase space
displayw1  displays time history of electric field, kinetic, and total
           energies

written by Viktor K. Decyk, UCLA, with contributions from Aileen Wang
update: September 21, 2023
"""
from __future__ import print_function
import sys
import tkinter as tk
from tkinter import ttk
import math
import cmath
import numpy
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
if (sys.version_info.major==2):
   from matplotlib.backends.backend_tkagg import \
      NavigationToolbar2TkAgg as NavigationToolbar2Tk
else:
   from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import comms

global win, num_wins, frame, fig, canvas, toolbar
num_wins = 0
win = []; frame = []; fig = []; canvas = []; toolbar = []

def main():
   global win, num_wins, frame, fig, canvas, toolbar

#Find location of control panel
   x = comms.root.winfo_x() + comms.root.winfo_width()
   y = comms.root.winfo_y()
#Find number of plots required from plot_loc dictionary
   num_plots = len(comms.plot_loc)
   num_wins = int((num_plots - 1)/4) + 1
#Default vertical and horizontal Figure size, in inches
   vfsz = 6; hfsz = 4
#Create windows
   for j in range(0,num_wins):
      nfs = min(4,max(0,num_plots-4*j))
      win.append(tk.Tk())
      win[j].title('Plot Window %d'%(j))
#Move window to the right of control panel
      win[j].geometry('+%d+%d'%(x+10+30*j,y+30*j))
      win[j].columnconfigure(0,weight=1)
      win[j].rowconfigure(0,weight=1)
#Find number of plots in window
      nfs = min(4,max(0,num_plots-4*j))
#Create frames in window
      for i in range(0,nfs):
         ii = i + 4*j
         frame.append(ttk.Frame(win[j]))
         nc = int(i/2)
         nr = i - 2*nc
         win[j].rowconfigure(nr, weight=1)
         win[j].columnconfigure(nc, weight=1)
         fig.append(Figure(figsize=(vfsz,hfsz),dpi=100))
         canvas.append(FigureCanvasTkAgg(fig[ii],master=frame[ii]))
         toolbar.append(NavigationToolbar2Tk(canvas[ii],frame[ii]))
         frame[ii].grid(row=nr,column=nc,columnspan=1,sticky='news')
         canvas[ii].get_tk_widget().pack(expand=True,fill=tk.BOTH)
         toolbar[ii].pack(expand=True,fill=tk.BOTH)
      win[j].update_idletasks()

#Send a message to start physics time loop
   comms.set_run_status()


#----------------------------------------------------------------------
def find_scales(fmax,fmin,ist):
   """
   This procedure calculates the plot scales, based on finding the
   minimum value power of 2 which will contain the plots, determined by
   the maximum and minimum (fmax,fmin) of data to be plotted.
   isc = power of 2 scale of y coordinate for plot
   In addition, the plot scales may be adjusted by the value of ist
   ist = flag for choosing positive and/or negative values
   if ist = 0, then ymax = 2**isc and ymin = -2**isc.
   if ist = 1, then ymax = 2**isc and ymin = 0.
   if ist = -1, then ymax = 0 and ymin = -2**isc.
   if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
   where fmin/fmax are the function minimum/maximum, 
   and ir = power of 2 scale for (fmax - fmin)
   function returns ymax, ymin
   """
   dv = 2.0
   algdvi = 1.0/math.log(dv)
# find maximum range
   rmax = fmax - fmin
   if (rmax==0.0):
         rmax = 1.0e-35
   rmin = fmin
# find power of 2 scale for maximum/minimum values
   ymax = abs(fmax)
   isc = int(math.log(ymax)*algdvi)
   if (ymax >= 1.0):
      isc += 1
   if (ymax <= dv**(isc-1)):
      isc -= 1
   ymin = abs(fmin)  
   if (ymin > 0.0):
      it = int(math.log(ymin)*algdvi)
      if (ymin >= 1.0):
         it += 1
      if (ymin <= dv**(it-1)):
         it -= 1
# find maximum/minimum values for plot
   if (fmax > 0.0):
      if (fmin > 0.0):
         fmin = dv**(it - 1)
      elif (fmin < 0.0):
         fmin = -dv**it
      fmax = dv**isc
   else:
      fmax = -dv**(isc - 1)
      fmin = -dv**it
# adjust for ist display option
   if (ist==0):
      if (ymin > ymax):
         fmax = dv**it
      else:
         fmin = -dv**isc
   elif (ist==2):
      ir = int(math.log(rmax)*algdvi)
      if (rmax >= 1.0):
          ir += 1
      if (rmax <= dv**(ir-1)):
         ir -= 1
      fmin = rmin
      fmax = rmin + dv**ir
      if (fmax==fmin):
          fmax = fmin + 1.0
   return fmax, fmin

#-----------------------------------------------------------------------
def dscaler1(f,label,isc,ist,nx,chr):
   """ displays 1d scalar field """
   xmin = 0.0; xmax = float(nx-1)
   disps(f,label,xmin,xmax,isc,ist,nx,chr)

#-----------------------------------------------------------------------
def displayfv1(fv,fvm,label,nmv):
   """
   displays velocity distribution functions
   fv = velocity distribution
   fvm = velocity moments
   label = long character string label for plot
   nmv = number of velocity intervals
   """
# internal parameters
# isc = 999 = use default display scale
# ist = 1 = display positive values only
# mks = 0 = cycle through line styles
   isc = 999; ist = 1; mks = 0
   nmvf = numpy.size(fv[:])
   nmv21 = 2*nmv + 1
# chrs = short array of characters to label individual line samples
   chrs = ['VX']
   vmax = fv[nmv21]
   vmin = -vmax
   chr = 'VD = ' + str(round(fvm[0],5)) + ', VTH = ' + \
                   str(round(fvm[1],5))
   disps(fv[:],label,vmin,vmax,isc,ist,nmv21,chr)

#-----------------------------------------------------------------------
def grasp1(part,fvms,label,isc,nx,iyp,ixp,npx,chr):
   """
   for 1 code, this subroutine displays (iyp-ixp) phase space
   plots background particles in blue and beam particles in red
   part(1,n) = position x of particle n
   part(2,n) = velocity vx of particle n
   fvms = scaled maximum of absolute values of velocities
   label = species label
   isc = power of 2 scale of range of values of velocity
   nx = system length in x direction
   iyp/ixp = phase space coordinates to be displayed
   npx = number of background particles
   chr = additional character string comment for plot
   """
   global fig, canvas, toolbar
# Check size of array
   idimp = numpy.size(part[:,0])
   np = numpy.size(part[0,:])
# Find window and frame number
   i = comms.plot_loc[label]
   j = int(i/4)
# dv = scale will be set in powers of this parameter
   dv = 2.0
# find scales for plot
   algdvi = 1.0/math.log(dv)
# find y scale for plot
   if (iyp <= 1):
      if (iyp < 1):
         iyp = 1
      ymin = 0.0
      ymax = float(nx)
# find x scale for plot
   if (ixp <= 1):
      if (ixp < 1):
         ixp = 1
      xmin = 0.0
      xmax = float(nx)
   iss = isc
# find velocity scales for plot
   jxp = ixp - 1
   jyp = iyp - 1
   if (abs(iss) > 116):
      if (iyp > 1):
         ymax = fvms[jyp-1]
         ymin = -ymax
      if (ixp > 1):
         xmax = fvms[jxp-1]
         xmin = -xmax
   else:
      ymax = dv**iss
      ymin = -ymax
      xmax = dv**iss
      xmin = -xmax
# initiate plot
   fig[i].clf()
   axes = fig[i].add_subplot(111)
# draw grid and labels
   fglabel = label + ', Time = ' + comms.stime
   axes.set_title(fglabel,fontsize=10)
   axes.set_xlabel(chr)
   axes.set_xlim(xmin=xmin,xmax=xmax)
   axes.set_ylim(ymin=ymin,ymax=ymax)
# plot data
   axes.plot(part[jxp,:npx],part[jyp,:npx],marker=',',lw=0,color='b')
   if (npx < np):
      axes.plot(part[jxp,npx:np],part[jyp,npx:np],marker=',',lw=0,
                color='r')

# update plot
   canvas[i].draw()
   canvas[i].get_tk_widget().update_idletasks()
   toolbar[i].update_idletasks()

# set status to done
   comms.set_plot_status(label)

#-----------------------------------------------------------------------
def displayw1(wt,label,t0,dtw,nt,movion):
   """
   displays time history of electric field, kinetic, and total energies
   wt = time history array for energies
   label = long character string label for plot
   t0 = initial energy time
   dtw = time between energy values
   nt = number of energy values to be displayed
   movion = number of moving ions
   """
# internal parameters
# isc = 999 = use default display scale
# ist = 2 = display minimum range
# mks = 0 = cycle through line styles
   isc = 999; ist = 2; mks = 0
# quit if array is empty or incorrect
   if (nt <= 0):
# set status to done
      comms.set_plot_status(label)
   ns = min(numpy.size(wt[0,:]),7)
# chrs = short array of characters to label individual line samples
# electrostatic labels
   if (movion==0):
     ns = ns - 1
     chrs = ['EL','ELECT','TOTAL']
   else:
      chrs = ['EL','ELECT','ION','TOTAL']
# location of total energy
   if (movion==0):
      ne = 2
   else:
      ne = 3
   chr = 'FIELD ENERGY='+str(round(wt[nt-1,0],2))+', TOTAL ENERGY='+\
         str(round(wt[nt-1,ne],0))
# tmin/tmax = range of time values in plot
   tmin = t0
   tmax = t0 + dtw*(nt - 1)
   if (nt==1):
      tmax = tmin + dtw
# all energies on common plot
   dispr(wt[:nt,0:ns]-wt[0,:ns],label,tmin,tmax,isc,ist,mks,nt,ns,chr,
         chrs)

#----------------------------------------------------------------------
def disps(f,label,xmin,xmax,isc,ist,nx,chr):
   """
   this subroutine plots an array f versus a linear function in x,
   where xmin < x < xmax.  It is plotted in solid line style.
   f = array to be plotted
   label = long character string label for plot
   xmin/xmax = range of x values in plot
   isc = power of 2 scale of y coordinate for plot
   ist = flag for choosing positive and/or negative values
   The plot has a scale in y given by ymax and ymin.
   if ist = 0, then ymax = 2**isc and ymin = -2**isc.
   if ist = 1, then ymax = 2**isc and ymin = 0.
   if ist = -1, then ymax = 0 and ymin = -2**isc.
   if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
   where fmin/fmax are the function minimum/maximum, 
   and ir = power of 2 scale for (fmax - fmin)
   if abs(isc) < 116, then the isc value passed is used for scale.
   if abs(isc) > 116, then the program finds the minimum value of isc
   which will contain the plot, determined by the absolute value of f
   nx = the number of points to be plotted, <= dimension of array f
   chr = additional long character string comment for plot
   """
   global fig, canvas, toolbar
# Check size of array
   nxv = numpy.size(f[:])
   if (nxv < nx):
      print('disps:invalid nx > nxv:nx,nxv=',nx,nxv)
      return
   if (xmax==xmin):
      xmax = xmin + 1.0e-35
# Find window and frame number
   i = comms.plot_loc[label]
   j = int(i/4)
# dv = scale will be set in powers of this parameter
   dv = 2.0
# find scales for plot
   algdvi = 1.0/math.log(dv)
   iss = isc
# find maximum/minimum of function
   if (abs(iss) > 116):
      fmax = numpy.amax(f[:nx])
      fmin = numpy.amin(f[:nx])
      if (fmax==0.0):
         fmax = 1.0e-35
# find ymax and ymin for display
      ymax, ymin = find_scales(fmax,fmin,ist)
   else:
      ymax = dv**iss
      ymin = -ymax
# check ist options
   if (ist==1):
      ymin = 0.0
   elif (ist==(-1)):
      ymax = 0.0
# initiate plot
   fig[i].clf()
   axes = fig[i].add_subplot(111)
# draw grid and labels
   fglabel = label + ', Time = ' + comms.stime
   axes.set_title(fglabel,fontsize=10)
   axes.set_xlabel(chr)
   axes.set_xlim(xmin=xmin,xmax=xmax)
   axes.set_ylim(ymin=ymin,ymax=ymax)
# plot data
   axes.plot(numpy.linspace(xmin,xmax,nx),f[:nx])

# update plot
   canvas[i].draw()
   canvas[i].get_tk_widget().update_idletasks()
   toolbar[i].update_idletasks()

# set status to done
   comms.set_plot_status(label)

#----------------------------------------------------------------------
def dispr(f,label,xmin,xmax,isc,ist,mks,nx,ngs,chr,chrs):
   """
   This subroutine plots ngs subarrays of the array f, on a common graph,
   each plot with nx points, versus a linear function in x,
   where xmin < x < xmax.
   f = array containing subarrays to be plotted
   label = long character string label for plot
   xmin/xmax = range of x values in plot
   isc = power of 2 scale of y coordinate for plot
   ist = flag for choosing positive and/or negative values
   the plots have a common scale in y given by ymax and ymin.
   if ist = 0, then ymax = 2**isc and ymin = -2**isc.
   if ist = 1, then ymax = 2**isc and ymin = 0.
   if ist = -1, then ymax = 0 and ymin = -2**isc.
   if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
   where fmin/fmax are the minimum/maximum of f, 
   and ir = power of 2 scale for (fmax - fmin)
   if abs(isc) < 116, then the isc value passed is used for scale.
   if abs(isc) > 116, then the program finds the minimum value of isc
   which will contain the plots, determined by the absolute value of f
   mks = flag to determine whether lines or markers are used,
   mks=0 means cycle through lines styles, mks > 0 means cycle through
   marker styles, using the value of mks as the initial marker,
   mks < 0 means draw the first subarray with a line, then subsequent
   subarrays with markers, using the value of abs(mks) as the initial
   marker.
   nx = number of points plotted in each subarray,
   must be <= first dimension of array f.
   ngs = number of subarrays to be plotted, 
   must be  <= second dimension of array f
   chr = additional long character string comment for plot
   chrs = array of ngs short character labels to label individual line
   or marker samples
   """
   global fig, canvas, toolbar
# Check size of array
   nxv = numpy.size(f[:,0])
   if (nxv < nx):
      print('dispr:invalid nx > nxv:nx,nxv=',nx,nxv)
      return
   ngv = numpy.size(f[0,:])
   if (ngv < ngs):
      print('dispr:invalid ngs > ngv:ngs,ngv=',ngs,ngv)
      return
   if (xmax==xmin):
      xmax = xmin + 1.0e-35
# Find window and frame number
   i = comms.plot_loc[label]
   j = int(i/4)
# dv = scale will be set in powers of this parameter
   dv = 2.0
# find scales for plot
   algdvi = 1.0/math.log(dv)
   iss = isc
# find maximum/minimum of function
   if (abs(iss) > 116):
      fmax = numpy.amax(f[:nx,:ngs])
      fmin = numpy.amin(f[:nx,:ngs])
      if (fmax==0.0):
         fmax = 1.0e-35
# find ymax and ymin for display
      ymax, ymin = find_scales(fmax,fmin,ist)
   else:
      ymax = dv**iss
      ymin = -ymax
# check ist options
   if (ist==1):
      ymin = 0.0
   elif (ist==(-1)):
      ymax = 0.0
# initiate plot
   fig[i].clf()
   axes = fig[i].add_subplot(111)
# draw grid and labels
   fglabel = label + ', Time = ' + comms.stime
   axes.set_title(fglabel,fontsize=10)
   axes.set_xlabel(chr)
   axes.set_xlim(xmin=xmin,xmax=xmax)
   axes.set_ylim(ymin=ymin,ymax=ymax)
# plot data
   x = numpy.linspace(xmin,xmax,nx)
   for k in range(0,ngs):
# use line types
      if (mks==0) or (mks < 0 and k==0):
         axes.plot(x,f[:nx,k],label=chrs[k])
# use markers
      else:
         axes.scatter(x,f[:nx,k],label=chrs[k])
# display sample lines or markers
      if (ngs <= 6):
         axes.legend(loc='upper right',fontsize=8)
      else:
          axes.legend(loc='upper right',fontsize=6)

# update plot
   canvas[i].draw()
   canvas[i].get_tk_widget().update_idletasks()
   toolbar[i].update_idletasks()

# set status to done
   comms.set_plot_status(label)

#-----------------------------------------------------------------------
def destroy_windows():
   """ callback for Quit """
   global win, num_wins
#Destroy windows
   for j in range(0,num_wins):
      win[j].destroy()

if (__name__=="__main__"):
   main()
