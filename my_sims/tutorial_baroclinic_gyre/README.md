Tutorial Example: "Baroclinic gyre"
(Baroclinic Ocean Gyre In Spherical Coordinates)
============================================================

Modifications done:
---------------------------------------------------------------------------------------------------------
Run with mpi (2 processors), output mnc, file/folder names, changed data diagnostics, time step, run time, grid size. 


SIZE.h file
---------------------------------------------------------------------------------------------------------

     &           sNx =  31,
     &           sNy =  62,
     &           OLx =   2,
     &           OLy =   2,
     &           nSx =   1,
     &           nSy =   1,
     &           nPx =   2,
     &           nPy =   1,
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  =   15)

packages.conf file
---------------------------------------------------------------------------------------------------------

      gfd
      diagnostics
      mnc

data.pkg file:
---------------------------------------------------------------------------------------------------------

      # Packages (lines beginning "#" are comments)
      &PACKAGES
      useMNC=.TRUE.,
      useDiagnostics=.TRUE.,
      &

data.diagnostics file
---------------------------------------------------------------------------------------------------------

      &DIAGNOSTICS_LIST
         dumpAtLast  = .TRUE.,
      #--
         fields(1:4,1) = 'UVEL','VVEL','WVEL','THETA',
         fileName(1) = 'outs_3D',
         frequency(1) = -86400.,
         timePhase(1) = -86400.,
         fileflags(1) = 'R',

      #--
         fields(1:2,2) = 'ETAN','PHIBOT',
         fileName(2) = 'outs_2D',
         frequency(2) = -86400.,
         timePhase(2) = -86400.,
         fileflags(2) = 'R',

      #--
         fields(1:3,3)  = 'ETAN    ','TRELAX  ','MXLDEPTH',
         fileName(3) = 'surfDiag',
         frequency(3) = -86400.,
         timePhase(3) = -86400.,
         fileflags(3) = 'R',

         fields(1:5,4)  = 'THETA   ','PHIHYD  ',
                        'UVEL    ','VVEL    ','WVEL    ',
      # did not specify levels => all levels are selected
         fileName(4) = 'dynDiag',
         frequency(4) = -86400.,
         timePhase(4) = -86400.,
         fileflags(4) = 'R',
      &

data.mnc file
---------------------------------------------------------------------------------------------------------

      # Example "data.mnc" file
      &MNC_01
      monitor_mnc=.FALSE.,
      # mnc_use_outdir=.TRUE.,
      # mnc_outdir_str='mnc_file_',
      &

data file
---------------------------------------------------------------------------------------------------------

      endTime = 31104000.,
      monitorFreq=25920.,



Configure and compile the code:
---------------------------------------------------------------------------------------------------------

      cd build
      ../../../tools/genmake2 -mods ../code -mpi -enable=mnc -optfile ../../../tools/build_options/darwin_amd64_gfortran
      make depend
      make
      cd ../



To run:
---------------------------------------------------------------------------------------------------------

      mkdir run01
      cd run01
      cp -r ../input/* .
      cp -r ../input .
      cp ../build/mitgcmuv .
      mpirun -np 2 ./mitgcmuv &


Comments:
---------------------------------------------------------------------------------------------------------

Run01: original tau as in tutorial 

      tauMax = 0.1;
      x = (xo-dx) : dx : xeast;
      y = (yo-dy/2) : dy : (ynorth+dy/2); 
      [X,Y] = ndgrid(x, y);  % zonal wind-stress on (XG,YC) points
      tau = -tauMax * cos(2*pi*((Y-yo)/(ny-2)/dy)); % ny-2 accounts for walls at N,S boundaries

Run02: different tau and x,y (took configuration from barotropic gyre)

      tauMax = 0.1;
      x = (-1:nx-2) / (nx-2);       % non-dim x-coordinate, located at XG points
      y = ((0:ny-1)-0.5) / (ny-2);  % non-dim y-coordinate, located at YC points
      [X,Y] = ndgrid(x, y);
      tau = tauMax * cos(pi*Y);

Run03: different tauMax

      tauMax = 0.5;
      x = (-1:nx-2) / (nx-2);       % non-dim x-coordinate, located at XG points
      y = ((0:ny-1)-0.5) / (ny-2);  % non-dim y-coordinate, located at YC points
      [X,Y] = ndgrid(x, y);
      tau = tauMax * cos(pi*Y);

Run04: different (negative) tau

      tauMax = 0.5;
      x = (-1:nx-2) / (nx-2);       % non-dim x-coordinate, located at XG points
      y = ((0:ny-1)-0.5) / (ny-2);  % non-dim y-coordinate, located at YC points
      [X,Y] = ndgrid(x, y);
      tau = -tauMax * cos(pi*Y);

Run05: different (X) tau

      tauMax = 0.5;
      x = (-1:nx-2) / (nx-2);       % non-dim x-coordinate, located at XG points
      y = ((0:ny-1)-0.5) / (ny-2);  % non-dim y-coordinate, located at YC points
      [X,Y] = ndgrid(x, y);
      tau = -tauMax * cos(pi*X);

Run06: tau same as run04, but  10yrs runtime instead of 1 yr

      tauMax = 0.5;
      x = (-1:nx-2) / (nx-2);       % non-dim x-coordinate, located at XG points
      y = ((0:ny-1)-0.5) / (ny-2);  % non-dim y-coordinate, located at YC points
      [X,Y] = ndgrid(x, y);
      tau = -tauMax * cos(pi*Y);

In data file:

      endTime = 311040000.,
      monitorFreq=259200.,

 In data.diagnostics file:
 
      frequency(n) = -31104000.,
      timePhase(n) = -31104000.,
              
Run06 took about 50 minutes!


Run07: tau same as run04, but  changed no_slips 

      tauMax = 0.5;
      x = (-1:nx-2) / (nx-2);       % non-dim x-coordinate, located at XG points
      y = ((0:ny-1)-0.5) / (ny-2);  % non-dim y-coordinate, located at YC points
      [X,Y] = ndgrid(x, y);
      tau = -tauMax * cos(pi*Y);

In data file:

      no_slip_sides=.FALSE.,
      no_slip_bottom=.TRUE.,
              


The different runs correspond to different tau and in the case of Run06 a longer runtime. MATLAB script gendata.m was used to change wind stress (tau). After changing tau the "run" steps were repeated.

Eta, U, V, and W were plotted using MATLAB script Baroclinic_Gyre_contour_plots.m. The script and plots are in /analysis folder. Run04, 06 and 07 also have hovmoller diagrams of Eta, U, V, and W to show the progression over a year, or 10 years at X=10. It looks like about 5 years are needed to reach steady state in Eta. 
