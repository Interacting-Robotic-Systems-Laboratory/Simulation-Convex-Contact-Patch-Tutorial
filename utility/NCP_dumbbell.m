function A = NCP_dumbbell(A)
unit = A.unit;
N = A.N;


global e_t e_o e_r mu;

e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;
mu =A.cof;


global m g;
m = A.mass;
g = A.gravity*unit;

global radius heg I_xx I_yy I_zz;

radius = A.dim(1)*unit;  % in fixed body frame's x direction
heg = A.dim(2)*unit;  % in fixed body frame's z direction

I_xx = (m/12)*(3*radius^2+heg^2);
I_yy = (m/12)*(3*radius^2+heg^2);
I_zz = (m/2)*radius^2;



global h;
h = 0.01;

% applied wrenches
% applied wrenches
global p_x p_y p_z p_xt p_yt p_zt;
p_x = 0;
p_y = 0;
p_z = 0;
p_xt = 0;
p_yt = 0;
p_zt = 0;

%% initial configuration and state
% q_old - position and orientation vector at l, q_old=[q_xo;q_yo;q_zo;q0_o;q1_o;q2_o;q3_o]
global q_old;

q_old(1:3,1) = A.initial_q(1:3)*unit;
q_old(4:7,1) = A.initial_q(4:7);
% nu_old - generalized velocity vector at l, nu_old=[v_xo;v_yo;v_zo;w_xo;w_yo;w_zo]
global nu_old;

nu_old(1:3,1) = A.initial_v(1:3)*unit;
nu_old(4:6,1) = A.initial_v(4:6); 


%% define the infinity and initial guess and fun to use
 A = initial_guess(A);
 l = A.l;
 u = A.u;
 Z = A.Z;
 
 fun = A.fun;
 Q = q_old;
 theta_z = 0;
for i=1:N
        p_x = 0;
        p_y = 0;
        p_z = 0;
        p_xt = 0;
        p_yt = 0;
        p_zt = 0;

    if i <= 3
        p_y = A.Total(i);
        p_xt = -p_y*heg/2;
    elseif (i> 100)&&(i<106)
        
        p_x = A.Total(i);
       
    elseif i>=106
       if i == 106
            w_y = nu_old(5);
       end
       
        w_z = nu_old(6);
        R = [cos(theta_z) sin(theta_z) 0;
            -sin(theta_z) cos(theta_z) 0;
            0 0 1];
        theta_z = theta_z+w_z*h;

        w_r = R*nu_old(4:6);
        A.w_r(:,i) = w_r;
       

        if abs(w_r(3))<0.3
            p_x = A.Total(i)*cos(theta_z);
            p_y = A.Total(i)*sin(theta_z);
            p_zt = A.Total(i)*heg/2;
        end
       
        
      
    end
    A.input(:,i) = [p_x;p_y;p_z;p_xt;p_yt;p_zt];
    
    tic;
    [A.z(:,i),~,~,~,status] = pathmcp(Z,l,u,fun);
    time_NCP = toc;
    if status == 1 
        A.time_NCP(i) = time_NCP;
    else
        A.time_NCP(i) = 0;
    end
    j = 1;
    if i == 1
        while status == 0
        j = j+1;
        Z = change_initial_guess(A,Z,Q);
        tic
        [A.z(:,i),~,~,~,status] = pathmcp(Z,l,u,fun);
        time_NCP = toc;
        if status == 1 
            A.time_NCP(i) = time_NCP;
        else
            A.time_NCP(i) = 0;
        end
            if j>=30
                error('Path can not found the solution, change your initial guess');
            end
        end
        
    else
        while status == 0
            j = j+1;
            Z = change_guess(A,Z,Q);
            tic
            [A.z(:,i),~,~,~,status] = pathmcp(Z,l,u,fun);
            time_NCP = toc;
            if status == 1 
                A.time_NCP(i) = time_NCP;
            else
                A.time_NCP(i) = 0;
            end
            if j>=30
                error('Path can not found the solution, change your initial guess');
            end
        end
    end
    
   

   [Q,Nu] = kinematic_map(q_old,A.z(:,i),h); % the function which returns the state vectors
   
   Z = A.z(:,i); % updating the initial guess for each iteration
   A.q(:,i) = Q; 
   
  

   q_old = Q; % updating the beginning value for next time step
   nu_old = Nu; 
   
   
end

end