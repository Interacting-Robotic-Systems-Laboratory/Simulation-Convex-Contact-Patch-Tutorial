function  movie_cylinder(A)
N = A.N;


r = A.dim(1)*A.unit;
h = A.dim(2)*A.unit;



    
    
    set(gca,'nextplot','replacechildren'); 
    v = VideoWriter('test.avi');
    open(v);
    for i = 1:N
    q0 = A.q(4,i);
    q1 = A.q(5,i);
    q2 = A.q(6,i);
    q3 = A.q(7,i);
    
    q_x = A.q(1,i);
    q_y = A.q(2,i);
    q_z = A.q(3,i);
    
    
    a1_x = A.z(7,i);
    a1_y = A.z(8,i);
    a1_z = A.z(9,i);
    
    
    
    E = [-q1 q0 -q3 q2;
     -q2 q3 q0 -q1;
     -q3 -q2 q1 q0];

    G = [-q1 q0 q3 -q2;
         -q2 -q3 q0 q1;
         -q3 q2 -q1 q0];
 
    R = E*G';
    
    H = [R,[q_x;q_y;q_z];zeros(1,3),1];
    
    
    
    
    [x,y,z] = cylinder_mod( H,r,h );
    
    
     a1=surf(x,y,z,'FaceColor',[1 1 0]);
     hold on
     
     
    axis('equal');
    axis([-0.5 1.5 -1.5 0.5 0 2]);
    xlabel('x (meter)');
    ylabel('y (meter)');
    zlabel('z (meter)');
    
    
    hold off
    frame = getframe(gcf);
    writeVideo(v,frame);
    end
close(v);
end




