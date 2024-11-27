#-----------------------------------------------------------------------
from __future__ import print_function
# Basic 1D Electrostatic PIC code
# written by Viktor K. Decyk, UCLA
import sys
import math
import numpy
import matplotlib.pyplot as plt

try:
   from pic1lib import *
except ModuleNotFoundError:
   print("Execute 'make python' to compile required libraries")
   exit()

from fomplib import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64

if (wpush1.fprecision()==0):
   float_type = numpy.float32
   complex_type = numpy.complex64
   print("single precision floating point library required")
elif (wpush1.fprecision()==1):
   float_type = numpy.float64
   complex_type = numpy.complex128

#-----------------------------------------------------------------------
def main():

# import GUI if present
   sys.path.append('./Basic_Gui')
   try:
      import comms
   except ModuleNotFoundError:
      GUI = False
   else:
      GUI = True
# indx = exponent which determines grid points in x direction:
# nx = 2**indx.
   indx =   9
# npx = number of electrons distributed in x direction.
   npx = 184320
# npxb = number of beam electrons in x direction
   npxb = 18432

# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
   tend = 100.0; dt = 0.2; qme = -1.0
# vtx = thermal velocity of electrons in x direction
# vx0 = drift velocity of electrons in x direction.
   vtx = 1.0; vx0 = 0.0
# vtdx = thermal velocity of beam electrons in x direction
# vdx = drift velocity of beam electrons in x direction
   vtdx = 0.5; vdx = 5.0
# ax = smoothed particle size in x direction
   ax = .912871
# idimp = number of particle coordinates = 2
   idimp = 2

# Multi-tasking Parameter:
# nvp = number of shared memory threads (0=default)
   nvp = 0

# External Electrostatic Traveling Wave Driver
# amodex = wave number which determines k0 = 2*pi*amodex/NX
# freq = frequency of external wave driver
# trmp = ramp-up time for external wave driver
# toff = shut-off time for external wave driver
   amodex = 0.0; freq = 0.0; trmp = 0.0; toff = 0.0
# el0/er0 = external pump amplitude for left-going/right-going wave
   el0 = 0.0; er0 = 0.0

# diagnostics parameters
# nte = number of time steps between electric field display
   nte = 10
# ntv = number of time steps between velocity-space diagnostic
# nmv = number of segments in v for velocity distribution
   ntv = 10; nmv = 40
# nts = number of time steps between phase space display
   nts = 10
# ntw = number of time steps between energy display
   ntw = 10

# wke/we/wtt = particle kinetic/electric field/total energy
   wke = numpy.zeros((1),float_type)
   we = numpy.zeros((1),float_type)
   wtt = numpy.zeros((1),float_type)

# declare and initialize timing data
   tinit = 0.0; tloop = 0.0
   itime = numpy.empty((4),numpy.int32)
   ltime = numpy.empty((4),numpy.int32)
   dtime = numpy.empty((1),double_type)
   tdpost = numpy.zeros((1),float_type)
   tpush = numpy.zeros((1),float_type)
   ws = numpy.zeros((1),float_type)

# start timing initialization
   dtimer(dtime,itime,-1)

# nvp = number of shared memory threads (0=default)
#  nvp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
   omplib.init_omp(nvp)
   nvp = omplib.getnthsize()

# initialize scalars for standard code
# np = total number of particles in simulation
# nx = number of grid points in x direction
   np = npx + npxb; nx = int(math.pow(2,indx)); nxh = int(nx/2)
   nxe = nx + 2; nxeh = int(nxe/2)
# nloop = number of time steps in simulation
# ntime = current time step
   nloop = int(tend/dt + .0001) + 1; ntime = 0
   qbme = qme
   affp = float(nx)/float(np)

# allocate data for standard code
# part = particle array
   part = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
   qe = numpy.empty((nxe),float_type,'F')
# fxe = smoothed electric field with guard cells
   fxe = numpy.empty((nxe),float_type,'F')
# ffc = form factor array for poisson solver
   ffc = numpy.empty((nxh),complex_type,'F')
# mixup = bit reverse table for FFT
   mixup = numpy.empty((nxh),int_type,'F')
# sct = sine/cosine table for FFT
   sct = numpy.empty((nxh),complex_type,'F')

# allocate data for diagnostics
   if (ntv > 0):
# fv = global electron velocity distribution functions
      fv = numpy.empty((2*nmv+2),float_type,'F')
# fvm = electron vdrift, vth, entropy for global distribution
      fvm = numpy.zeros((3),float_type,'F')
# fvmx = scaled maximum of absolute values of electron velocities
      fvmx = numpy.zeros((2),float_type,'F')
   if (nts > 0):
# fvms = scaled maximum of absolute values of electron velocities
      fvms = numpy.zeros((2),float_type,'F')
# wt = energy time history array
   if (ntw > 0):
      mtw = int((nloop - 1)/ntw + 1); itw = 0; dtw = dt*float(ntw)
      wt = numpy.zeros((mtw,4),float_type,'F')

# prepare fft tables
   wpush1.fft1_init(mixup,sct,indx)
# calculate form factors
   wpush1.pois1_init(ffc,ax,affp,nx)
# initialize electrons
# uniform density in space, maxwellian in velocity space: updates part
   if (npx > 0) :
      wpush1.wdistr1(part,vtx,vx0,npx,nx)
# initialize beam electrons
# uniform density in space, maxwellian in velocity space: updates part
   if (npxb > 0):
      wpush1.wdistr1(part[:,npx:np],vtdx,vdx,npxb,nx)

# initialization time
   dtimer(dtime,itime,1)
   tinit += float(dtime[0])
# start timing loop
   dtimer(dtime,ltime,-1)

# start GUI if active
   if (GUI):
# send plot_loc, sim_data dictionaries and custom event to initialize
# plots
      import comms
      try:
         plot_loc = {}; plotnum = 0
         if (nte > 0):
            plot_loc['SMOOTHED ELECTRIC FIELD VS X'] = plotnum
            plotnum += 1
         if (ntv > 0):
            plot_loc['ELECT VELOCITY DISTR VS VX'] = plotnum
            plotnum += 1
         if (nts > 0):
            plot_loc['ELECTRON VX VERSUS X'] = plotnum
            plotnum += 1
         if (ntw > 0):
            plot_loc['ENERGY DIFFS VS TIME'] = plotnum
            plotnum += 1
         comms.update_gui(plot_loc,{"DT": dt,"TEND": tend,"NX": nx,
                          "NMV": nmv,"NPX": npx,"DTW": dtw})
# wait for plot intialization
         gui_err = comms.check_run_status()
         if (gui_err=='QUIT'):
            exit()
# Run without GUI
      except AttributeError:
         GUI = False

# * * * start main iteration loop * * *

   for ntime in range(0,nloop):
#     print("ntime = ", ntime)
# send current time to GUI, if active
      if GUI:
        comms.update_time(abs(dt*float(ntime)))

# deposit charge with standard procedure: updates qe
      dtimer(dtime,itime,-1)
      qe.fill(0.0)
      wpush1.wgpost1(part,qe,qme)
      dtimer(dtime,itime,1)
      tdpost[0] += float(dtime[0])

# add guard cells with standard procedure: updates qe
      wpush1.waguard1(qe,nx)

# transform charge to fourier space with standard procedure:
# updates qe, fxe
      isign = -1
      wpush1.wfft1r(qe,isign,mixup,sct,indx)

# calculate force/charge in fourier space with standard procedure:
# updates fxe, we
      wpush1.wpois1(qe,fxe,ffc,we,nx)

# transform force to real space with standard procedure: updates fxe, qe
      isign = 1
      wpush1.wfft1r(fxe,isign,mixup,sct,indx)

# add external traveling wave field: updates fxe
      ts = abs(dt)*float(ntime)
      wpush1.meaddext1(fxe,amodex,freq,ts,trmp,toff,el0,er0,nx)

# copy guard cells with standard procedure: updates fxe
      wpush1.wcguard1(fxe,nx)

# display efield in real space
      if (GUI):
         if (nte > 0):
            if (ntime==nte*int(ntime/nte)):
               comms.update_plot('SMOOTHED ELECTRIC FIELD VS X',fxe)
               gui_err = comms.check_plot_status()
               if (gui_err=='QUIT'):
                  break

# velocity diagnostic
      if (GUI):
         if (ntv > 0):
            if (ntime==ntv*int(ntime/ntv)):
# find maximum of electron velocity and energy
               wdiag1.wvmxt1(part,fvmx)
               fv[2*nmv+1] = fvmx[0]
               wdiag1.wvdist1(part,fv,fvm,nmv)
               comms.update_plot('ELECT VELOCITY DISTR VS VX',fv,fvm)
               gui_err = comms.check_plot_status()
               if (gui_err=='QUIT'):
                  break

# phase space diagnostic
# plot electrons vx versus x
      if (GUI):
         if (nts > 0):
            if (ntime==nts*int(ntime/nts)):
# find maximum of electron velocity and energy
               wdiag1.wvmxt1(part,fvms)
               comms.update_plot('ELECTRON VX VERSUS X',part,fvms)
               gui_err = comms.check_plot_status()
               if (gui_err=='QUIT'):
                  break

# push particles: updates part, wke
      wke[0] = 0.0
      dtimer(dtime,itime,-1)
# use parallel version
      if (nvp > 1):
         wpush1.wgmpush1(part,fxe,qbme,dt,wke,nvp,nx)
# use serial version
      else:
         wpush1.wgpush1(part,fxe,qbme,dt,wke,nx)
      dtimer(dtime,itime,1)
      tpush[0] += float(dtime[0])

      if (ntime==0):
         print("Initial Field, Kinetic and Total Energies:")
         wt0 = wke[0] + we[0]
         print("%14.7e %14.7e %14.7e" % (we[0], wke[0], wt0))

# energy diagnostic: updates wt
      if (ntw > 0):
         if (ntime==ntw*int(ntime/ntw)):
            wt[itw,:] = [we[0],wke[0],0.0,we[0]+wke[0]]
            itw += 1
            if (GUI):
               comms.update_plot('ENERGY DIFFS VS TIME',wt,itw)
               gui_err = comms.check_plot_status()
               if (gui_err=='QUIT'):
                  break

# wait for GUI each iteration step if active
      if (GUI):
# skip last interation
         if ((ntime+1) < nloop):
            gui_err = comms.check_run_status()
            if (gui_err=='QUIT'):
               break

   ntime = ntime + 1

# check if GUI needs to quit
   if (GUI):
      if (gui_err=='QUIT'):
         exit()

# loop time
   dtimer(dtime,ltime,1)
   tloop += float(dtime[0])

# * * * end main iteration loop * * *

   print("ntime = ", ntime)
   print("Final Field, Kinetic and Total Energies:")
   wtt[0] = wke[0] + we[0]
   print("%14.7e %14.7e %14.7e" % (we[0], wke[0], wtt[0]))
   print("")

   print("Energy Conservation(%) = ",100.0*(wtt[0] - wt0)/wt0)

   print("")
   ws[0] = tinit
   print("initialization time = ", ws[0])
   print("deposit time = ", tdpost[0])
   print("push time = ", tpush[0])
   wtt[0] = tdpost[0] + tpush[0]
   print("total particle time = ", wtt[0])
   ws[0] = tloop
   print("total time = ", ws[0])
   print("")

   ws[0] = 1.0e+09/(float(nloop)*float(np))
   print("Push Time (nsec) = ", tpush[0]*ws[0])
   print("Deposit Time (nsec) = ", tdpost[0]*ws[0])
   print("Total Particle Time (nsec) = ", wtt[0]*ws[0])

if (__name__=="__main__"):
   main()
