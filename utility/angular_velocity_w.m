function w = angular_velocity_w(q,w_b)
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);


E = [-q1 q0 -q3 q2;
 -q2 q3 q0 -q1;
 -q3 -q2 q1 q0];

G = [-q1 q0 q3 -q2;
     -q2 -q3 q0 q1;
     -q3 q2 -q1 q0];
 
R = E*G';

w = R*w_b;

end