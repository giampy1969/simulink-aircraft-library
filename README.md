# simulink-aircraft-library
Simulink&reg; library of aircraft models

[![View simulink-aircraft-library on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/3019-simulink-aircraft-library)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=giampy1969/simulink-aircraft-library)

Airlib is a library of aircraft models to be used with Simulink 3.1 (MATLAB 5.3) or later (tested until 2018a).

It is based on two basic blocks that implement a continuous-time and a discrete-time nonlinear general aircraft model.
The initial version, dated Feb 2003, contains 13 different nonlinear aircraft models including, among others, a Boeing 747, an F4 Phantom, an F104 Starfighter, an X15, an IAI Pioneer (unmanned), a Convair 880 and four different Cessna models.

Each model is a link to one of the two library blocks, but of course with different parameters on the mask. Indeed, the type of aircraft is entirely specified by the inertial and aerodynamic parameters typed on the mask, which also includes a very detailed help that describes the meaning of inputs, states, outputs and mask parameters.

An example showing a Boeing 747 in a (closed loop) straight and level flight is also given.

In May 2003 a very effective MATLAB&reg; function (air3m.m) was added to allow the user to easily trim any given aircraft for any combination of desired speed, altitude, and flight path angle. This is based on the jj_3m , jj_lin, and other functions from J. J. Buchholz and Daniel Kiehn (thank you J.J and Daniel).

In March 04, a function that computes the aerodynamic derivatives from the linear model matrices was added, together with a Simulink model that compares two different ways of handling the wind in the equations of motions.

In July 2006 a full guidance and control system based on the feedback linearization of the aircraft kinematic and dynamics was added. The scheme includes a waypoint generation system, and can be easily adapted to any aircraft with known aerodynamics coefficients.

In November 2024 Paolo Massioni replaced the older air2m helper functions with new ones (thank you Paolo).

Please have a look to the readme.txt file for more detailed info.
Giampiero Campa, March 2024
