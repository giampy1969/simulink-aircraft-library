function [CD,CL,CY,Cl,Cm,Cn]=ab2dv(Alg,Blg,Alt,Blt,GM,deltaE_0,alpha_0,V0,h0,T0)

% [CD,CL,CY,Cl,Cm,Cn]=ab2dv(Alg,Blg,Alt,Blt,GM1,deltaE_0,alpha_0,V0,h0,T0)
% Recovers the aerodynamic derivatives given linear model and inertial data.
% For the longitudinal matrices Alg and Blg, states are v (m/s), alpha (rad) 
% and q (rad/s), the input is deltaE (rad).
% For the lateral matrices Alt and Blt, states are beta (rad), p (rad/s),
% and r (rad/s), the inputs are deltaA (rad) and deltaR (rad).
% Therefore, Alg and Alt are 3 by 3 , Blg is 3 by 1 and Blt is 3 by 2.
% The 10 by 1 vector GM is the vector containing all the geometric data 
% that is cbar, b, S, Ix, Iz, Jxy, Jxz, Jyz and mass (see the FDC manual also).
% The last 5 parameters are deltaE (rad), alpha (rad), airspeed (m/s), 
% altitude (m) and Trust (N), relative to the aircraft trim condition.

c_bar = GM(1);
b     = GM(2);
S     = GM(3);
Ix    = GM(4);
Iy    = GM(5);
Iz    = GM(6);
Jxy   = GM(7);
Jxz   = GM(8);
Jyz   = GM(9);
m     = GM(10);

% inertia moment coefficients calculation : took directly from fdc
detI = Ix*Iy*Iz - 2*Jxy*Jxz*Jyz - Ix*Jyz^2 - Iy*Jxz^2 - Iz*Jxy^2;
I1   = Iy*Iz - Jyz^2;I2   = Jxy*Iz + Jyz*Jxz;I3   = Jxy*Jyz + Iy*Jxz;
I4   = Ix*Iz - Jxz^2;I5   = Ix*Jyz + Jxy*Jxz;I6   = Ix*Iy - Jxy^2;
Pl  = I1/detI; Pm = I2/detI; Pn = I3/detI;
Ppp = -(Jxz*I2 - Jxy*I3)/detI;Ppq = (Jxz*I1 - Jyz*I2 - (Iy-Ix)*I3)/detI;
Ppr = -(Jxy*I1 + (Ix-Iz)*I2 - Jyz*I3)/detI;Pqq = (Jyz*I1 - Jxy*I3)/detI;
Pqr = -((Iz-Iy)*I1 - Jxy*I2 + Jxz*I3)/detI;Prr = -(Jyz*I1 - Jxz*I2)/detI;
Ql  = I2/detI; Qm = I4/detI; Qn = I5/detI;
Qpp = -(Jxz*I4 - Jxy*I5)/detI;Qpq = (Jxz*I2 - Jyz*I4 - (Iy-Ix)*I5)/detI;
Qpr = -(Jxy*I2 + (Ix-Iz)*I4 - Jyz*I5)/detI;Qqq = (Jyz*I2 - Jxy*I5)/detI;
Qqr = -((Iz-Iy)*I2 - Jxy*I4 + Jxz*I5)/detI;Qrr = -(Jyz*I2 - Jxz*I4)/detI;
Rl  = I3/detI; Rm = I5/detI; Rn = I6/detI;
Rpp = -(Jxz*I5 - Jxy*I6)/detI;
Rpq = (Jxz*I3 - Jyz*I5 - (Iy-Ix)*I6)/detI;
Rpr = -(Jxy*I3 + (Ix-Iz)*I5 - Jyz*I6)/detI;
Rqq = (Jyz*I3 - Jxy*I6)/detI;
Rqr = -((Iz-Iy)*I3 - Jxy*I5 + Jxz*I6)/detI;
Rrr = -(Jyz*I3 - Jxz*I5)/detI;

g0 = 9.80665;                 % m/s^2,     gravitational acceleration at sea level
Tt0 = 288.15+0;               % K,         temperature
R_earth = 6371020;            % m,         radius of the Earth

g = g0*(R_earth/(R_earth+h0))^2;    %   m/s^2
Tt = Tt0-0.0065*h0;
rho0 = 101325*(Tt/Tt0)^(g/1.86584)/(287.05*Tt);
q0 = 1/2*rho0*V0^2;

C_D_alpha = (g-Alg(1,2))*m/(q0*S);
C_D_q = -Alg(1,3)*(2*m*V0)/(q0*S*c_bar);
C_D_deltaE = -Blg(1)*m/(q0*S);
C_D_0 = T0*cos(alpha_0)/(q0*S) - C_D_alpha*alpha_0 - C_D_deltaE*deltaE_0;

C_L_alpha  = -(m*V0*Alg(2,2)+T0)/(q0*S);
C_L_q      = (1-Alg(2,3))*(2*m*V0^2)/(q0*S*c_bar);
C_L_deltaE = -Blg(2)*(m*V0)/(q0*S);
C_L_0 = (m*g-T0*sin(alpha_0))/(q0*S) - C_L_alpha*alpha_0 - C_L_deltaE*deltaE_0;

C_m_alpha  = Alg(3,2)/(Qm*q0*S*c_bar);
C_m_q        = Alg(3,3)*(2*V0)/(Qm*q0*S*c_bar^2);
C_m_deltaE = Blg(3)/(Qm*q0*S*c_bar);
C_m_0        = -C_m_alpha*alpha_0-C_m_deltaE*deltaE_0;

C_y_beta    = Alt(1,1)*(m*V0)/(q0*S);
C_y_p        = (Alt(1,2)-sin(alpha_0))*(2*m*V0^2)/(q0*S*b);
C_y_r        = (Alt(1,3)+cos(alpha_0))*(2*m*V0^2)/(q0*S*b);
C_y_deltaA = Blt(1,1)*(m*V0)/(q0*S);
C_y_deltaR = Blt(1,2)*(m*V0)/(q0*S);

InvPRln = inv([Pl, Pn; Rl, Rn]);

temp = 1/(q0*S*b)*InvPRln*[Alt(2,1); Alt(3,1)];
C_l_beta = temp(1);     C_n_beta = temp(2);
temp = 2*V0/(q0*S*b^2)*InvPRln*[Alt(2,2); Alt(3,2)];
C_l_p = temp(1);         C_n_p = temp(2);
temp = 2*V0/(q0*S*b^2)*InvPRln*[Alt(2,3); Alt(3,3)];
C_l_r = temp(1);          C_n_r = temp(2);
temp = 1/(q0*S*b)*InvPRln*[Blt(2,1); Blt(3,1)];
C_l_deltaA = temp(1);  C_n_deltaA = temp(2);
temp = 1/(q0*S*b)*InvPRln*[Blt(2,2); Blt(3,2)];
C_l_deltaR = temp(1);  C_n_deltaR = temp(2);

CD = [C_D_0  C_D_alpha  C_D_q  C_D_deltaE 0.0];
CL = [C_L_0  C_L_alpha  C_L_q  C_L_deltaE 0.0];
Cm = [C_m_0  C_m_alpha  C_m_q  C_m_deltaE 0];

Cl = [0   C_l_beta  C_l_p  C_l_r  C_l_deltaA  C_l_deltaR];
Cn = [0   C_n_beta  C_n_p  C_n_r  C_n_deltaA  C_n_deltaR];
CY = [0   C_y_beta  C_y_p  C_y_r  C_y_deltaA  C_y_deltaR];

% Copyright 2018 The MathWorks, Inc.