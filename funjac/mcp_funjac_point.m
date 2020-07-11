function [F, J, domerr] = mcp_funjac_point(z, jacflag)

%% initialize
z = z(:);
F = [];
J = [];
domerr = 0;

%% obtain value of global variables
global h;

global q_old ;
q_xo = q_old(1);
q_yo = q_old(2);
q_zo = q_old(3);

global nu_old;
v_xo =nu_old(1);
v_yo =nu_old(2);
v_zo =nu_old(3);

global m mu e_o e_t g;

global p_x p_y p_z;

%% unknown variables

v_x = z(1);
v_y = z(2);
v_z = z(3);

p_t = z(4);
p_o = z(5);

sig = z(6);
p_n = z(7);


%% intermediate variables appear in the chain rule
% z - vector of unknown variables 

% The intermediate variables are:
% 1. [qo(z) q1(z) q2(z) q3(z)] - Quaternions, which depends on z
% 2. Rij(qo,q1,q2,q3) - elements in rotation matrix, which depends on the quaternion
% 3. Iij(qo,q1,q2,q3) - elements in RIR', which depends on the quaternion
% 4. q_x(z),q_y(z),q_z(z) - position vector in l+1, which depends on q_old and z

% Thus system of equations F can be described by the unknown variables z or
% with intermediate variables:
% 1. F(z,Rij,Iij,q_x,q_y,q_z) - F descriabed by z with intermediate variables
% 2. F(z) - F described by unknown variables z

% Based on the chain rule, Jocobain J = J1*J2*J3 where:
% J - dF(z)/dz Jacobian matrix 
% J1 = dF(z,Rij,Iij,q_x,q_y,q_z)/d[z,Rij,Iij,q_x,q_y,q_z] 
% J2 = d[z,Rij(qo,q1,q2,q3),Iij(qo,q1,q2,q3),q_x(z),q_y(z),q_z(z)]/d[z,qo,q1,q2,q3] 
% J3 = d[z,qo(z),q1(z),q2(z),q3(z)]/dz 

% q0 = -(2*((h*q1_o*w_x)/2 - q0_o + (h*q2_o*w_y)/2 + (h*q3_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
% q1 = (2*(q1_o + (h*q0_o*w_x)/2 + (h*q3_o*w_y)/2 - (h*q2_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
% q2 = (2*(q2_o - (h*q3_o*w_x)/2 + (h*q0_o*w_y)/2 + (h*q1_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
% q3 = (2*(q3_o + (h*q2_o*w_x)/2 - (h*q1_o*w_y)/2 + (h*q0_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);

q_x = q_xo+h*v_x;
q_y = q_yo+h*v_y;
q_z = q_zo+h*v_z;

% R11 = q0^2 + q1^2 - q2^2 - q3^2;
% R12 = 2*q1*q2 - 2*q0*q3;
% R13 = 2*q0*q2 + 2*q1*q3;
% R21 = 2*q0*q3 + 2*q1*q2;
% R22 = q0^2 - q1^2 + q2^2 - q3^2;
% R23 = 2*q2*q3 - 2*q0*q1;
% R31 = 2*q1*q3 - 2*q0*q2;
% R32 = 2*q0*q1 + 2*q2*q3;
% R33 = q0^2 - q1^2 - q2^2 + q3^2;
% 
% I11 = I_xx*R11^2 + I_yy*R12^2 + I_zz*R13^2;
% I12 = I_xx*R21*R11 + I_yy* R22*R12 + I_zz*R23*R13;
% I13 = I_xx*R31*R11 + I_yy* R32*R12 + I_zz*R33*R13;
% I21 = I12;
% I22 = I_xx*R21^2 + I_yy*R22^2 + I_zz*R23^2;
% I23 = I_xx*R31*R21 + I_yy*R32*R22 + I_zz*R33*R23;
% I31 = I13;
% I32 = I23;
% I33 = I_xx*R31^2 + I_yy*R32^2 + I_zz*R33^2;

%% system of equations f(x) == 0
% Newton_Euler equations
F(1) = p_t + p_x - m*(v_x - v_xo);
F(2) = p_o + p_y - m*(v_y - v_yo);
F(3) = p_n + p_z - m*(v_z - v_zo) - g*h*m;


% Friction Model without complementarity equation

F(4) = p_t*sig + e_t^2*mu*p_n*v_x;
F(5) = p_o*sig + e_o^2*mu*p_n*v_y;

% contact constraints without complementarity equation


% complementarity equation in Friction Model
F(6) = mu^2*p_n^2 - p_t^2/e_t^2 - p_o^2/e_o^2;

% complementarity equations in contact constraints

% non-penetration 
F(7) = q_z;


%% Jacobian matrix (Using chain rules) dF/d[z,q_x,q_y,q_z]*d[z,q_x,q_y,q_z]/dz

if (jacflag)
    J1 = [-m 0 0 1 0 0 0 0 0 0;
          0 -m 0 0 1 0 0 0 0 0;
          0  0 -m 0 0 0 1 0 0 0;
          e_t^2*mu*p_n 0 0 sig 0 p_t e_t^2*mu*v_x 0 0 0;
          0 e_o^2*mu*p_n 0 0 sig p_o e_o^2*mu*v_y 0 0 0;
          0 0 0 -2*p_t/e_t^2 -2*p_o/e_o^2 0 2*mu^2*p_n 0 0 0;
          0 0 h 0 0 0 0 0 0 1];
      
     J2 = [eye(7);
          h 0 0 0 0 0 0;
          0 h 0 0 0 0 0;
          0 0 h 0 0 0 0];
      
      J = J1*J2;
      J = sparse(J);
    
end
end