function [F, J, domerr] = mcp_funjac_single_convex_contact_patch_ellipsoid(z, jacflag)

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
q0_o = q_old(4);
q1_o = q_old(5);
q2_o = q_old(6);
q3_o = q_old(7);

global nu_old;
v_xo =nu_old(1);
v_yo =nu_old(2);
v_zo =nu_old(3);
w_xo =nu_old(4);
w_yo =nu_old(5);
w_zo =nu_old(6);

global m mu I_xx I_yy I_zz e_o e_r e_t g;

global r_a r_b r_c;
et = r_a; eo = r_b; er = r_c;
global p_x p_y p_z p_xt p_yt p_zt;

%% unknown variables

v_x = z(1);
v_y = z(2);
v_z = z(3);
w_x = z(4);
w_y = z(5);
w_z = z(6);

a1_x = z(7);
a1_y = z(8);
a1_z = z(9);

a2_x = z(10);
a2_y = z(11);
a2_z = z(12);

p_t = z(13);
p_o = z(14);
p_r = z(15);

sig = z(16);

l1 = z(17);
l2 = z(18);

p_n = z(19);


%% intermediate variables appear in the chain rule
% z - vector of unknown variables 

% The intermediate variables are:
% 1. [qo(z) q1(z) q2(z) q3(z)] - Quaternions, which depends on z
% 2. Rij(qo,q1,q2,q3) - elements in rotation matrix, which depends on the quaternion
% 3. Iij(qo,q1,q2,q3) - elements in RIR', which depends on the quaternion
% 3. q_x(z),q_y(z),q_z(z) - position vector in l+1, which depends on q_old and z

% Thus system of equations F can be described by the unknown variables z or
% with intermediate variables:
% 1. F(z,Rij,Iij,q_x,q_y,q_z) - F descriabed by z with intermediate variables
% 2. F(z) - F described by unknown variables z

% Based on the chain rule, Jocobain J = J1*J2*J3 where:
% J - dF(z)/dz Jacobian matrix 
% J1 = dF(z,Rij,Iij,q_x,q_y,q_z)/d[z,Rij,q_x,q_y,q_z] 
% J2 = d[z,Rij(qo,q1,q2,q3),Iij(qo,q1,q2,q3),q_x(z),q_y(z),q_z(z)]/d[z,qo,q1,q2,q3] 
% J3 = d[z,qo(z),q1(z),q2(z),q3(z)]/dz 

%% intermediate variables
q0 = -(2*((h*q1_o*w_x)/2 - q0_o + (h*q2_o*w_y)/2 + (h*q3_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
q1 = (2*(q1_o + (h*q0_o*w_x)/2 + (h*q3_o*w_y)/2 - (h*q2_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
q2 = (2*(q2_o - (h*q3_o*w_x)/2 + (h*q0_o*w_y)/2 + (h*q1_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
q3 = (2*(q3_o + (h*q2_o*w_x)/2 - (h*q1_o*w_y)/2 + (h*q0_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);

q_x = q_xo+h*v_x;
q_y = q_yo+h*v_y;
q_z = q_zo+h*v_z;

R11 = q0^2 + q1^2 - q2^2 - q3^2;
R12 = 2*q1*q2 - 2*q0*q3;
R13 = 2*q0*q2 + 2*q1*q3;
R21 = 2*q0*q3 + 2*q1*q2;
R22 = q0^2 - q1^2 + q2^2 - q3^2;
R23 = 2*q2*q3 - 2*q0*q1;
R31 = 2*q1*q3 - 2*q0*q2;
R32 = 2*q0*q1 + 2*q2*q3;
R33 = q0^2 - q1^2 - q2^2 + q3^2;

I11 = I_xx*R11^2 + I_yy*R12^2 + I_zz*R13^2;
I12 = I_xx*R21*R11 + I_yy* R22*R12 + I_zz*R23*R13;
I13 = I_xx*R31*R11 + I_yy* R32*R12 + I_zz*R33*R13;
I21 = I12;
I22 = I_xx*R21^2 + I_yy*R22^2 + I_zz*R23^2;
I23 = I_xx*R31*R21 + I_yy*R32*R22 + I_zz*R33*R23;
I31 = I13;
I32 = I23;
I33 = I_xx*R31^2 + I_yy*R32^2 + I_zz*R33^2;

%% simplified variables
rx = R11*a1_x + R21*a1_y + R31*a1_z - R11*q_x - R21*q_y - R31*q_z;
ry = R12*a1_x + R22*a1_y + R32*a1_z - R12*q_x - R22*q_y - R32*q_z;
rz = R13*a1_x + R23*a1_y + R33*a1_z - R13*q_x - R23*q_y - R33*q_z;

%% Newton_Euler equations
F(1) = p_t + p_x - m*(v_x - v_xo);
F(2) = p_o + p_y - m*(v_y - v_yo);
F(3) = p_n + p_z - m*(v_z - v_zo) - g*h*m;

F(4) = p_xt - h*(w_y*(I31*w_x + I32*w_y + I33*w_z) - w_z*(I21*w_x + I22*w_y + I23*w_z)) - I11*(w_x - w_xo) - I12*(w_y - w_yo) - I13*(w_z - w_zo) + p_n*(a1_y - q_y) - p_o*(a1_z - q_z);
F(5) = p_yt + h*(w_x*(I31*w_x + I32*w_y + I33*w_z) - w_z*(I11*w_x + I12*w_y + I13*w_z)) - I21*(w_x - w_xo) - I22*(w_y - w_yo) - I23*(w_z - w_zo) - p_n*(a1_x - q_x) + p_t*(a1_z - q_z);
F(6) = p_r + p_zt - h*(w_x*(I21*w_x + I22*w_y + I23*w_z) - w_y*(I11*w_x + I12*w_y + I13*w_z)) - I31*(w_x - w_xo) - I32*(w_y - w_yo) - I33*(w_z - w_zo) + p_o*(a1_x - q_x) - p_t*(a1_y - q_y);

%% Friction Model without complementarity equation

F(7) = p_t*sig + e_t^2*mu*p_n*v_x - e_t^2*mu*p_n*w_z*(a1_y - q_y) + e_t^2*mu*p_n*w_y*(a1_z - q_z);
F(8) = p_o*sig + e_o^2*mu*p_n*v_y + e_o^2*mu*p_n*w_z*(a1_x - q_x) - e_o^2*mu*p_n*w_x*(a1_z - q_z);
F(9) = mu*p_n*w_z*e_r^2 + p_r*sig;

%% contact constraints without complementarity equation

F(10) = a2_x - a1_x;
F(11) = a2_y - a1_y;
F(12) = a2_z - a1_z + l2;

F(13) = l1*((2*R12*ry)/eo^2 + (2*R11*rx)/et^2 + (2*R13*rz)/er^2);
F(14) = l1*((2*R22*ry)/eo^2 + (2*R21*rx)/et^2 + (2*R23*rz)/er^2);
F(15) = l1*((2*R32*ry)/eo^2 + (2*R31*rx)/et^2 + (2*R33*rz)/er^2) + 1;

%% complementarity equation in Friction Model
F(16) = mu^2*p_n^2 - p_r^2/e_r^2 - p_t^2/e_t^2 - p_o^2/e_o^2;

%% complementarity equations in contact constraints
F(17) =  1 - rz^2/er^2 - rx^2/et^2 - ry^2/eo^2;
F(18) = -a2_z;
%% non-penetration 
F(19) = a1_z;

%% Jacobian matrix (Using chain rules)

if (jacflag)
    
    % Newton_Euler equations
    J1 = [           -m,            0,  0,                                                   0,                                                   0,                                                   0,                                                              0,                                                              0,                                                              0, 0, 0,  0,              1,              0,              0,   0,                                                      0, 0,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,           -m,  0,                                                   0,                                                   0,                                                   0,                                                              0,                                                              0,                                                              0, 0, 0,  0,              0,              1,              0,   0,                                                      0, 0,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0, -m,                                                   0,                                                   0,                                                   0,                                                              0,                                                              0,                                                              0, 0, 0,  0,              0,              0,              0,   0,                                                      0, 0,                                                                    1,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0,  0,                       - I11 - h*(I31*w_y - I21*w_z), - I12 - h*(I31*w_x + 2*I32*w_y - I22*w_z + I33*w_z),   h*(I21*w_x + I22*w_y - I33*w_y + 2*I23*w_z) - I13,                                                              0,                                                            p_n,                                                           -p_o, 0, 0,  0,              0,     q_z - a1_z,              0,   0,                                                      0, 0,                                                           a1_y - q_y,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0, w_xo - w_x, w_yo - w_y, w_zo - w_z,  h*w_x*w_z,  h*w_y*w_z,    h*w_z^2, -h*w_x*w_y,   -h*w_y^2, -h*w_y*w_z,                                                               0,                                                            -p_n,                                                             p_o
                      0,            0,  0,   h*(2*I31*w_x + I32*w_y - I11*w_z + I33*w_z) - I21,                         h*(I32*w_x - I12*w_z) - I22, - I23 - h*(I11*w_x - I33*w_x + I12*w_y + 2*I13*w_z),                                                           -p_n,                                                              0,                                                            p_t, 0, 0,  0,     a1_z - q_z,              0,              0,   0,                                                      0, 0,                                                           q_x - a1_x,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0, -h*w_x*w_z, -h*w_y*w_z,   -h*w_z^2, w_xo - w_x, w_yo - w_y, w_zo - w_z,    h*w_x^2,  h*w_x*w_y,  h*w_x*w_z,                                                             p_n,                                                               0,                                                            -p_t
                      0,            0,  0, - I31 - h*(2*I21*w_x - I11*w_y + I22*w_y + I23*w_z),   h*(I11*w_x - I22*w_x + 2*I12*w_y + I13*w_z) - I32,                       - I33 - h*(I23*w_x - I13*w_y),                                                            p_o,                                                           -p_t,                                                              0, 0, 0,  0,     q_y - a1_y,     a1_x - q_x,              1,   0,                                                      0, 0,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,  h*w_x*w_y,    h*w_y^2,  h*w_y*w_z,   -h*w_x^2, -h*w_x*w_y, -h*w_x*w_z, w_xo - w_x, w_yo - w_y, w_zo - w_z,                                                            -p_o,                                                             p_t,                                                               0
           e_t^2*mu*p_n,            0,  0,                                                   0,                           e_t^2*mu*p_n*(a1_z - q_z),                          -e_t^2*mu*p_n*(a1_y - q_y),                                                              0,                                              -e_t^2*mu*p_n*w_z,                                               e_t^2*mu*p_n*w_y, 0, 0,  0,            sig,              0,              0, p_t,                                                      0, 0, e_t^2*mu*v_x - e_t^2*mu*w_z*(a1_y - q_y) + e_t^2*mu*w_y*(a1_z - q_z),                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                e_t^2*mu*p_n*w_z,                                               -e_t^2*mu*p_n*w_y
                      0, e_o^2*mu*p_n,  0,                          -e_o^2*mu*p_n*(a1_z - q_z),                                                   0,                           e_o^2*mu*p_n*(a1_x - q_x),                                               e_o^2*mu*p_n*w_z,                                                              0,                                              -e_o^2*mu*p_n*w_x, 0, 0,  0,              0,            sig,              0, p_o,                                                      0, 0, e_o^2*mu*v_y + e_o^2*mu*w_z*(a1_x - q_x) - e_o^2*mu*w_x*(a1_z - q_z),                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                               -e_o^2*mu*p_n*w_z,                                                               0,                                                e_o^2*mu*p_n*w_x
                      0,            0,  0,                                                   0,                                                   0,                                        e_r^2*mu*p_n,                                                              0,                                                              0,                                                              0, 0, 0,  0,              0,              0,            sig, p_r,                                                      0, 0,                                                         e_r^2*mu*w_z,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0,  0,                                                   0,                                                   0,                                                   0,                                                             -1,                                                              0,                                                              0, 1, 0,  0,              0,              0,              0,   0,                                                      0, 0,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0,  0,                                                   0,                                                   0,                                                   0,                                                              0,                                                             -1,                                                              0, 0, 1,  0,              0,              0,              0,   0,                                                      0, 0,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0,  0,                                                   0,                                                   0,                                                   0,                                                              0,                                                              0,                                                             -1, 0, 0,  1,              0,              0,              0,   0,                                                      0, 1,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0,  0,                                                   0,                                                   0,                                                   0,       l1*((2*R11^2)/r_a^2 + (2*R12^2)/r_b^2 + (2*R13^2)/r_c^2), l1*((2*R11*R21)/r_a^2 + (2*R12*R22)/r_b^2 + (2*R13*R23)/r_c^2), l1*((2*R11*R31)/r_a^2 + (2*R12*R32)/r_b^2 + (2*R13*R33)/r_c^2), 0, 0,  0,              0,              0,              0,   0, (2*R11*rx)/r_a^2 + (2*R12*ry)/r_b^2 + (2*R13*rz)/r_c^2, 0,                                                                    0, (2*l1*(rx + R11*a1_x - R11*q_x))/r_a^2, (2*l1*(ry + R12*a1_x - R12*q_x))/r_b^2, (2*l1*(rz + R13*a1_x - R13*q_x))/r_c^2,          (2*R11*l1*(a1_y - q_y))/r_a^2,          (2*R12*l1*(a1_y - q_y))/r_b^2,          (2*R13*l1*(a1_y - q_y))/r_c^2,          (2*R11*l1*(a1_z - q_z))/r_a^2,          (2*R12*l1*(a1_z - q_z))/r_b^2,          (2*R13*l1*(a1_z - q_z))/r_c^2,          0,          0,          0,          0,          0,          0,          0,          0,          0,       -l1*((2*R11^2)/r_a^2 + (2*R12^2)/r_b^2 + (2*R13^2)/r_c^2), -l1*((2*R11*R21)/r_a^2 + (2*R12*R22)/r_b^2 + (2*R13*R23)/r_c^2), -l1*((2*R11*R31)/r_a^2 + (2*R12*R32)/r_b^2 + (2*R13*R33)/r_c^2)
                      0,            0,  0,                                                   0,                                                   0,                                                   0, l1*((2*R11*R21)/r_a^2 + (2*R12*R22)/r_b^2 + (2*R13*R23)/r_c^2),       l1*((2*R21^2)/r_a^2 + (2*R22^2)/r_b^2 + (2*R23^2)/r_c^2), l1*((2*R21*R31)/r_a^2 + (2*R22*R32)/r_b^2 + (2*R23*R33)/r_c^2), 0, 0,  0,              0,              0,              0,   0, (2*R21*rx)/r_a^2 + (2*R22*ry)/r_b^2 + (2*R23*rz)/r_c^2, 0,                                                                    0,          (2*R21*l1*(a1_x - q_x))/r_a^2,          (2*R22*l1*(a1_x - q_x))/r_b^2,          (2*R23*l1*(a1_x - q_x))/r_c^2, (2*l1*(rx + R21*a1_y - R21*q_y))/r_a^2, (2*l1*(ry + R22*a1_y - R22*q_y))/r_b^2, (2*l1*(rz + R23*a1_y - R23*q_y))/r_c^2,          (2*R21*l1*(a1_z - q_z))/r_a^2,          (2*R22*l1*(a1_z - q_z))/r_b^2,          (2*R23*l1*(a1_z - q_z))/r_c^2,          0,          0,          0,          0,          0,          0,          0,          0,          0, -l1*((2*R11*R21)/r_a^2 + (2*R12*R22)/r_b^2 + (2*R13*R23)/r_c^2),       -l1*((2*R21^2)/r_a^2 + (2*R22^2)/r_b^2 + (2*R23^2)/r_c^2), -l1*((2*R21*R31)/r_a^2 + (2*R22*R32)/r_b^2 + (2*R23*R33)/r_c^2)
                      0,            0,  0,                                                   0,                                                   0,                                                   0, l1*((2*R11*R31)/r_a^2 + (2*R12*R32)/r_b^2 + (2*R13*R33)/r_c^2), l1*((2*R21*R31)/r_a^2 + (2*R22*R32)/r_b^2 + (2*R23*R33)/r_c^2),       l1*((2*R31^2)/r_a^2 + (2*R32^2)/r_b^2 + (2*R33^2)/r_c^2), 0, 0,  0,              0,              0,              0,   0, (2*R31*rx)/r_a^2 + (2*R32*ry)/r_b^2 + (2*R33*rz)/r_c^2, 0,                                                                    0,          (2*R31*l1*(a1_x - q_x))/r_a^2,          (2*R32*l1*(a1_x - q_x))/r_b^2,          (2*R33*l1*(a1_x - q_x))/r_c^2,          (2*R31*l1*(a1_y - q_y))/r_a^2,          (2*R32*l1*(a1_y - q_y))/r_b^2,          (2*R33*l1*(a1_y - q_y))/r_c^2, (2*l1*(rx + R31*a1_z - R31*q_z))/r_a^2, (2*l1*(ry + R32*a1_z - R32*q_z))/r_b^2, (2*l1*(rz + R33*a1_z - R33*q_z))/r_c^2,          0,          0,          0,          0,          0,          0,          0,          0,          0, -l1*((2*R11*R31)/r_a^2 + (2*R12*R32)/r_b^2 + (2*R13*R33)/r_c^2), -l1*((2*R21*R31)/r_a^2 + (2*R22*R32)/r_b^2 + (2*R23*R33)/r_c^2),       -l1*((2*R31^2)/r_a^2 + (2*R32^2)/r_b^2 + (2*R33^2)/r_c^2)
                      0,            0,  0,                                                   0,                                                   0,                                                   0,                                                              0,                                                              0,                                                              0, 0, 0,  0, -(2*p_t)/e_t^2, -(2*p_o)/e_o^2, -(2*p_r)/e_r^2,   0,                                                      0, 0,                                                           2*mu^2*p_n,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0,  0,                                                   0,                                                   0,                                                   0,       - (2*R11*rx)/r_a^2 - (2*R12*ry)/r_b^2 - (2*R13*rz)/r_c^2,       - (2*R21*rx)/r_a^2 - (2*R22*ry)/r_b^2 - (2*R23*rz)/r_c^2,       - (2*R31*rx)/r_a^2 - (2*R32*ry)/r_b^2 - (2*R33*rz)/r_c^2, 0, 0,  0,              0,              0,              0,   0,                                                      0, 0,                                                                    0,             -(2*rx*(a1_x - q_x))/r_a^2,             -(2*ry*(a1_x - q_x))/r_b^2,             -(2*rz*(a1_x - q_x))/r_c^2,             -(2*rx*(a1_y - q_y))/r_a^2,             -(2*ry*(a1_y - q_y))/r_b^2,             -(2*rz*(a1_y - q_y))/r_c^2,             -(2*rx*(a1_z - q_z))/r_a^2,             -(2*ry*(a1_z - q_z))/r_b^2,             -(2*rz*(a1_z - q_z))/r_c^2,          0,          0,          0,          0,          0,          0,          0,          0,          0,          (2*R11*rx)/r_a^2 + (2*R12*ry)/r_b^2 + (2*R13*rz)/r_c^2,          (2*R21*rx)/r_a^2 + (2*R22*ry)/r_b^2 + (2*R23*rz)/r_c^2,          (2*R31*rx)/r_a^2 + (2*R32*ry)/r_b^2 + (2*R33*rz)/r_c^2
                      0,            0,  0,                                                   0,                                                   0,                                                   0,                                                              0,                                                              0,                                                              0, 0, 0, -1,              0,              0,              0,   0,                                                      0, 0,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0
                      0,            0,  0,                                                   0,                                                   0,                                                   0,                                                              0,                                                              0,                                                              1, 0, 0,  0,              0,              0,              0,   0,                                                      0, 0,                                                                    0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,          0,          0,          0,          0,          0,          0,          0,          0,          0,                                                               0,                                                               0,                                                               0];

    
    J2 = zeros(40,23);
    J2(1:19,1:19) = eye(19); %dz/dz
    J2(20:28,end-3:end) = [2*q0,  2*q1, -2*q2, -2*q3
                          -2*q3,  2*q2,  2*q1, -2*q0
                           2*q2,  2*q3,  2*q0,  2*q1
                           2*q3,  2*q2,  2*q1,  2*q0
                           2*q0, -2*q1,  2*q2, -2*q3
                          -2*q1, -2*q0,  2*q3,  2*q2
                          -2*q2,  2*q3, -2*q0,  2*q1
                           2*q1,  2*q0,  2*q3,  2*q2
                           2*q0, -2*q1, -2*q2,  2*q3]; %dRij/d[q0,q1,q2,q3]
                       
    J2(29:37,end-3:end) = [                                                 4*I_xx*R11*q0 - 4*I_yy*R12*q3 + 4*I_zz*R13*q2,                                                 4*I_xx*R11*q1 + 4*I_yy*R12*q2 + 4*I_zz*R13*q3,                                                 4*I_yy*R12*q1 - 4*I_xx*R11*q2 + 4*I_zz*R13*q0,                                                 4*I_zz*R13*q1 - 4*I_yy*R12*q0 - 4*I_xx*R11*q3
                            2*I_xx*R11*q3 + 2*I_xx*R21*q0 + 2*I_yy*R12*q0 - 2*I_yy*R22*q3 - 2*I_zz*R13*q1 + 2*I_zz*R23*q2, 2*I_xx*R11*q2 + 2*I_xx*R21*q1 - 2*I_yy*R12*q1 + 2*I_yy*R22*q2 - 2*I_zz*R13*q0 + 2*I_zz*R23*q3, 2*I_xx*R11*q1 - 2*I_xx*R21*q2 + 2*I_yy*R12*q2 + 2*I_yy*R22*q1 + 2*I_zz*R13*q3 + 2*I_zz*R23*q0, 2*I_xx*R11*q0 - 2*I_xx*R21*q3 - 2*I_yy*R12*q3 - 2*I_yy*R22*q0 + 2*I_zz*R13*q2 + 2*I_zz*R23*q1
                            2*I_xx*R31*q0 - 2*I_xx*R11*q2 + 2*I_yy*R12*q1 - 2*I_yy*R32*q3 + 2*I_zz*R13*q0 + 2*I_zz*R33*q2, 2*I_xx*R11*q3 + 2*I_xx*R31*q1 + 2*I_yy*R12*q0 + 2*I_yy*R32*q2 - 2*I_zz*R13*q1 + 2*I_zz*R33*q3, 2*I_yy*R12*q3 - 2*I_xx*R31*q2 - 2*I_xx*R11*q0 + 2*I_yy*R32*q1 - 2*I_zz*R13*q2 + 2*I_zz*R33*q0, 2*I_xx*R11*q1 - 2*I_xx*R31*q3 + 2*I_yy*R12*q2 - 2*I_yy*R32*q0 + 2*I_zz*R13*q3 + 2*I_zz*R33*q1
                            2*I_xx*R11*q3 + 2*I_xx*R21*q0 + 2*I_yy*R12*q0 - 2*I_yy*R22*q3 - 2*I_zz*R13*q1 + 2*I_zz*R23*q2, 2*I_xx*R11*q2 + 2*I_xx*R21*q1 - 2*I_yy*R12*q1 + 2*I_yy*R22*q2 - 2*I_zz*R13*q0 + 2*I_zz*R23*q3, 2*I_xx*R11*q1 - 2*I_xx*R21*q2 + 2*I_yy*R12*q2 + 2*I_yy*R22*q1 + 2*I_zz*R13*q3 + 2*I_zz*R23*q0, 2*I_xx*R11*q0 - 2*I_xx*R21*q3 - 2*I_yy*R12*q3 - 2*I_yy*R22*q0 + 2*I_zz*R13*q2 + 2*I_zz*R23*q1
                                                                            4*I_xx*R21*q3 + 4*I_yy*R22*q0 - 4*I_zz*R23*q1,                                                 4*I_xx*R21*q2 - 4*I_yy*R22*q1 - 4*I_zz*R23*q0,                                                 4*I_xx*R21*q1 + 4*I_yy*R22*q2 + 4*I_zz*R23*q3,                                                 4*I_xx*R21*q0 - 4*I_yy*R22*q3 + 4*I_zz*R23*q2
                            2*I_xx*R31*q3 - 2*I_xx*R21*q2 + 2*I_yy*R22*q1 + 2*I_yy*R32*q0 + 2*I_zz*R23*q0 - 2*I_zz*R33*q1, 2*I_xx*R21*q3 + 2*I_xx*R31*q2 + 2*I_yy*R22*q0 - 2*I_yy*R32*q1 - 2*I_zz*R23*q1 - 2*I_zz*R33*q0, 2*I_xx*R31*q1 - 2*I_xx*R21*q0 + 2*I_yy*R22*q3 + 2*I_yy*R32*q2 - 2*I_zz*R23*q2 + 2*I_zz*R33*q3, 2*I_xx*R21*q1 + 2*I_xx*R31*q0 + 2*I_yy*R22*q2 - 2*I_yy*R32*q3 + 2*I_zz*R23*q3 + 2*I_zz*R33*q2
                            2*I_xx*R31*q0 - 2*I_xx*R11*q2 + 2*I_yy*R12*q1 - 2*I_yy*R32*q3 + 2*I_zz*R13*q0 + 2*I_zz*R33*q2, 2*I_xx*R11*q3 + 2*I_xx*R31*q1 + 2*I_yy*R12*q0 + 2*I_yy*R32*q2 - 2*I_zz*R13*q1 + 2*I_zz*R33*q3, 2*I_yy*R12*q3 - 2*I_xx*R31*q2 - 2*I_xx*R11*q0 + 2*I_yy*R32*q1 - 2*I_zz*R13*q2 + 2*I_zz*R33*q0, 2*I_xx*R11*q1 - 2*I_xx*R31*q3 + 2*I_yy*R12*q2 - 2*I_yy*R32*q0 + 2*I_zz*R13*q3 + 2*I_zz*R33*q1
                            2*I_xx*R31*q3 - 2*I_xx*R21*q2 + 2*I_yy*R22*q1 + 2*I_yy*R32*q0 + 2*I_zz*R23*q0 - 2*I_zz*R33*q1, 2*I_xx*R21*q3 + 2*I_xx*R31*q2 + 2*I_yy*R22*q0 - 2*I_yy*R32*q1 - 2*I_zz*R23*q1 - 2*I_zz*R33*q0, 2*I_xx*R31*q1 - 2*I_xx*R21*q0 + 2*I_yy*R22*q3 + 2*I_yy*R32*q2 - 2*I_zz*R23*q2 + 2*I_zz*R33*q3, 2*I_xx*R21*q1 + 2*I_xx*R31*q0 + 2*I_yy*R22*q2 - 2*I_yy*R32*q3 + 2*I_zz*R23*q3 + 2*I_zz*R33*q2
                                                                            4*I_yy*R32*q1 - 4*I_xx*R31*q2 + 4*I_zz*R33*q0,                                                 4*I_xx*R31*q3 + 4*I_yy*R32*q0 - 4*I_zz*R33*q1,                                                 4*I_yy*R32*q3 - 4*I_xx*R31*q0 - 4*I_zz*R33*q2,                                                 4*I_xx*R31*q1 + 4*I_yy*R32*q2 + 4*I_zz*R33*q3];
 
    %dIij/d[q0,q1,q2,q3]                                                                    
    J2(38:40,1:3) = h*eye(3);                   
    
    J3 = zeros(23,19);  
    J3(1:19,1:19) = eye(19);  %dz/dz
    N_J3 = (h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2);
    J3(20:23,4:6) = [ -(h*(q1_o*h^2*w_y^2 - q2_o*w_x*h^2*w_y + q1_o*h^2*w_z^2 - q3_o*w_x*h^2*w_z + 2*q0_o*w_x*h + 4*q1_o))/N_J3, -(h*(q2_o*h^2*w_x^2 - q1_o*w_y*h^2*w_x + q2_o*h^2*w_z^2 - q3_o*w_y*h^2*w_z + 2*q0_o*w_y*h + 4*q2_o))/N_J3, -(h*(q3_o*h^2*w_x^2 - q1_o*w_z*h^2*w_x + q3_o*h^2*w_y^2 - q2_o*w_z*h^2*w_y + 2*q0_o*w_z*h + 4*q3_o))/N_J3
                       (h*(q0_o*h^2*w_y^2 - q3_o*w_x*h^2*w_y + q0_o*h^2*w_z^2 + q2_o*w_x*h^2*w_z - 2*q1_o*w_x*h + 4*q0_o))/N_J3,  (h*(q3_o*h^2*w_x^2 - q0_o*w_y*h^2*w_x + q3_o*h^2*w_z^2 + q2_o*w_y*h^2*w_z - 2*q1_o*w_y*h + 4*q3_o))/N_J3, -(h*(q2_o*h^2*w_x^2 + q0_o*w_z*h^2*w_x + q2_o*h^2*w_y^2 + q3_o*w_z*h^2*w_y + 2*q1_o*w_z*h + 4*q2_o))/N_J3
                      -(h*(q3_o*h^2*w_y^2 + q0_o*w_x*h^2*w_y + q3_o*h^2*w_z^2 + q1_o*w_x*h^2*w_z + 2*q2_o*w_x*h + 4*q3_o))/N_J3,  (h*(q0_o*h^2*w_x^2 + q3_o*w_y*h^2*w_x + q0_o*h^2*w_z^2 - q1_o*w_y*h^2*w_z - 2*q2_o*w_y*h + 4*q0_o))/N_J3,  (h*(q1_o*h^2*w_x^2 + q3_o*w_z*h^2*w_x + q1_o*h^2*w_y^2 - q0_o*w_z*h^2*w_y - 2*q2_o*w_z*h + 4*q1_o))/N_J3
                       (h*(q2_o*h^2*w_y^2 + q1_o*w_x*h^2*w_y + q2_o*h^2*w_z^2 - q0_o*w_x*h^2*w_z - 2*q3_o*w_x*h + 4*q2_o))/N_J3, -(h*(q1_o*h^2*w_x^2 + q2_o*w_y*h^2*w_x + q1_o*h^2*w_z^2 + q0_o*w_y*h^2*w_z + 2*q3_o*w_y*h + 4*q1_o))/N_J3,  (h*(q0_o*h^2*w_x^2 - q2_o*w_z*h^2*w_x + q0_o*h^2*w_y^2 + q1_o*w_z*h^2*w_y - 2*q3_o*w_z*h + 4*q0_o))/N_J3];%d(q0 q1 q2 q3)/d(wx wy wz)
 
    J = J1*J2*J3;
    J = sparse(J);
    
end
end
