function  movie_ellipse(A)
N = A.N;
unit = A.unit;

range =1*[-1 1 -1 1 -1 1];
r_a = A.dim(1)*unit;  % in fixed body frame's x direction
r_b = A.dim(2)*unit;  % in fixed body frame's y direction
r_c = A.dim(3)*unit;  % in fixed body frame's z direction

N_edge = 10;
[x, y, z] = ellipsoid(0,0,0,r_a,r_b,r_c,N_edge);  
x = reshape(x,1,[]);
y = reshape(y,1,[]);
z = reshape(z,1,[]);
set(gca,'nextplot','replacechildren'); 
v = VideoWriter('test.avi');
open(v);
for i = 1:1:N
 
    
    q0 = A.q(4,i);
    q1 = A.q(5,i);
    q2 = A.q(6,i);
    q3 = A.q(7,i);
    
    q_x = A.q(1,i);
    q_y = A.q(2,i);
    q_z = A.q(3,i);
    
    E = [-q1 q0 -q3 q2;
     -q2 q3 q0 -q1;
     -q3 -q2 q1 q0];

    G = [-q1 q0 q3 -q2;
         -q2 -q3 q0 q1;
         -q3 q2 -q1 q0];
 
    R = E*G';
    
    H = [R,[q_x;q_y;q_z];zeros(1,3),1];
    
    W = H*[x;y;z;ones(1,(N_edge+1)*(N_edge+1))];
    x_w = reshape(W(1,:),(N_edge+1),(N_edge+1));
    y_w = reshape(W(2,:),(N_edge+1),(N_edge+1));
    z_w = reshape(W(3,:),(N_edge+1),(N_edge+1));
    surf(x_w,y_w,z_w,'FaceColor',[1 1 0],'Linestyle','-'); 
    hold on
    plot3(A.z(7,i),A.z(8,i),A.z(9,i),'r.');
    hold on
    [X,Y] = meshgrid(range(1):range(2)-range(1):range(2),range(3):range(4)-range(3):range(4));
	Z = X.*0+Y.*0;
    surf(X,Y,Z,'FaceColor','k','FaceAlpha',0.1,'Linestyle','-');  
    
    axis('equal');
    axis(range);
    xlabel('x (meter)');
    ylabel('y (meter)');
    zlabel('z (meter)');
    
    view([9,3]);
    hold off
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);
end




