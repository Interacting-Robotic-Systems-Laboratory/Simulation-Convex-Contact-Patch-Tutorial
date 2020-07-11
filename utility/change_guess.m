function Z = change_guess(A,Z,Q)
unit = A.unit;


   
e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;
mu =A.cof;


m = A.mass;
g = A.gravity*unit;
h = A.h;

R = rand;


nu = Z(1:6); 
v_x = nu(1);
v_y = nu(2);
w_x = nu(4);
w_y = nu(5);
w_z = nu(6);

q_x = Q(1);
q_y = Q(2);
q_z = Q(3);

if length(A.initial_q) == 3
    Z = R*Z(7);
    return
end

a1 = Z(7:9);
a2 = a1;
ECP = [a1;a2];
  
p_n = R*m*g*h; % from 0 to mgh


v_t = v_x - w_z*(p_n*(a1(1) - q_x))/p_n;
v_o = v_y + w_z*(p_n*(a1(2) - q_y))/p_n;
v_r = w_z;
    
    
    
sig = sqrt(e_t^2*v_t^2 + e_o^2*v_o^2 + e_r^2*v_r^2);
    

    
    
p_t = -e_t^2*mu*p_n*v_t/sig;
p_o = -e_o^2*mu*p_n*v_o/sig;
p_r = -e_r^2*mu*p_n*v_r/sig;
Con_wrench = [p_t;p_o;p_r];      
    
    
    

    
  
  
if size(A.dim,2) == 2   
    La = [0;0;0;1];
     
    Z = [nu;ECP;Con_wrench;sig;La;p_n];
elseif size(A.dim,2) == 3
    La = [0;0;0;1;0;0;0];
    Z = [nu;ECP;Con_wrench;sig;La;p_n];
end



end