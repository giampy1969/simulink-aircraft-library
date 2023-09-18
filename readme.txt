Airlib is a library of aircraft models to be used with Simulink 3.1 (Matlab 5.3) or later.

The library consists in just two nonlinear general aircraft models, a continuous-time one and its discrete-time approximation. 

Both blocks are largely based on blocks from the FDC toolbox, written by Marc Rauw and available at http://dutchroll.sourceforge.net/fdc.html

The type of aircraft is entirely specified by the inertial and aerodynamic parameters typed on the mask. The mask has a very detaliled help that illustrates the 
meaning of inputs, states, outputs and mask parameters.

The library comes with an example, the simulink file airlibex.mdl, that shows 13 different aircraft models, including among others, a Boeing 747, an F4 Phantom, an F104 Starfighter, an X15, an IAI Pioneer (unmanned), a Convair 880 and four different Cessna models. Each model is a link to one of the two library blocks, but of course with different parameters on the mask.

Many parameters were taken by this website http://www.aae.uiuc.edu/m-selig/apasim/Aircraft-uiuc.html mantained by Michael S. Selig, Rob Deters, and Glen Dimock at the university of Urbana-Champaign. 

An example showing a Boeing 747 in a (closed loop) straight and level flight is also given. Finally, the file airtrim.mdl could be easily invoked to find a trim point for a given aircraft. 

Altough the trim points were computed using the entire trimmod utility, in May 2003 a very effective matlab function (air3m.m) was added to allow the user to easily trim any given aircraft for any combination of desired speed, altitude, and flight path angle. Moreover, when those desired conditions are supplied as matrices, then the input and state vectors for the resulting trim points are given as 5D matrices that are ready to be used with the block "interpolate matrix" from the Aerospace Blockset. Modifying the function to give as output also the 5D matrices A,B,C,D to be used with the "3D controller" block, (also from the Aerospace Blockset) should be immediate.

In February 2004, the function ab2dv.m and the model fdcwind.mdl were added. The first is a function that recovers the aerodynamic derivatives of an aircraft given its inertial data, some trim point values, and the longitudinal and lateral A and B matrices. The second is a file that compares two different but equivalent ways to handle the wind in the equations. This also proves that the way in which FDC handles the wind is fully correct provided that the right initial conditions are chosen for V,alpha,beta whenever the initial value of the wind is nonzero.

In July 2006 a full guidance and control system based on the feedback linearization of the aircraft kinematic and dynamics was added. The scheme includes a waypoint generation system, and can be easily adapted to any aicraft with known aerodynamics coefficients. Note that baing based on feedback linearization (a.k.a dynamic inversion) this controller is general in theory but is in practice often fragile especially when the dynamics is close to being nonminimum phase.

NOTE: in Aug 2017 the old xyz scope graph (sfunxyz.m) used in airgk.mdl was replaced with the new function sfun3d.m whichworks better but only for 2014b and later.
If you need to use version prior to 2014b you can either just delete the the 3D graph block, or if you need 3D trajectory visualization you can use the model airgkr11, after installing the 3DScope package from here: https://www.mathworks.com/matlabcentral/fileexchange/4915

This library should be an useful complement to the Aerospace Blockset, from The Mathworks, (http://www.mathworks.com/products/aeroblks/). Even for those who use the Aerosim Blockset, (http://www.u-dynamics.com/aerosim/default.htm) the library could serve as an useful collection of drag-and-drop ready-to-use nonlinear aircraft models.

This is an OPEN project, so feel free to contribute in any way, in particular, if you have inertial and aerodynamic data of an aircraft that is not in the library, you can send it to me and i will be pleased to put it in the library along with your name.

Giampiero Campa, Feb 2004  
Copyright 2018 The MathWorks, Inc.