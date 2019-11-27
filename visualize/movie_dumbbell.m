function movie_dumbbell(A)
N = A.N;


r = A.dim(1)*A.unit; %radius of the bell
height = A.dim(2)*A.unit;
r1 = 0.5*r; %radius of the bar
heightr = 0.2*height;
q_t = [0;0;height/2-heightr/2;1];
q_b = [0;0;-height/2+heightr/2;1];
    
q_t1 = [0;0;height/2;1];
q_1 = [0;0;height/2-heightr/2;1];    
q_b1 = [0;0;height/2-heightr;1]; 

q_t2 = [0;0;-height/2+heightr;1];
q_2 = [0;0;-height/2+heightr/2;1];
q_b2 = [0;0;-height/2;1];

a1_c  = [0;r;-height/2;1];
a2_c  = [0;-r;-height/2;1];
range =0.5*[-0.3 2.7 -0.3 2.7 0 1.5];



set(gca,'nextplot','replacechildren'); 
v = VideoWriter('test.avi');
open(v);
for i = 1:5:250
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
    
    q_w1 = H*q_1;
    q_w2 = H*q_2;
    
    q_wb1 = H*q_b1;
    q_wb2 = H*q_b2;
  
    q_wt1 = H*q_t1;
    q_wt2 = H*q_t2;
    
    a1_cw = H*a1_c;
    a2_cw = H*a2_c;
    %Cylinder
    [x1,y1,z1]=cylinder_mod([q_w1(1),q_w1(2),q_w1(3)],r,heightr);
    [x2,y2,z2]=cylinder_mod([q_w2(1),q_w2(2),q_w2(3)],r,heightr);
    [x3,y3,z3]=cylinder_mod([q_x,q_y,q_z],r1,height-heightr);
    
    %ellipsoid
    [x_b1, y_b1, z_b1] = ellipsoid(q_wb1(1),q_wb1(2),q_wb1(3),r,r,0,30);
    [x_b2, y_b2, z_b2] = ellipsoid(q_wb2(1),q_wb2(2),q_wb2(3),r,r,0,30);
    
    [x_t1, y_t1, z_t1] = ellipsoid(q_wt1(1),q_wt1(2),q_wt1(3),r,r,0,30);
    [x_t2, y_t2, z_t2] = ellipsoid(q_wt2(1),q_wt2(2),q_wt2(3),r,r,0,30);
    
    
    
    
    a1=surf(x1,y1,z1,'FaceColor',[1 1 0]);
    hold on
    a2=surf(x2,y2,z2,'FaceColor',[1 1 0]);
    hold on
    a3=surf(x3,y3,z3,'FaceColor',[1 1 0]);
    hold on
    a4=surf(x_b1,y_b1,z_b1,'FaceColor',[1 1 0],'EdgeColor','none');
    hold on
    a5=surf(x_b2,y_b2,z_b2,'FaceColor',[1 1 0],'EdgeColor','none');
    hold on
    a6=surf(x_t1,y_t1,z_t1,'FaceColor',[1 1 0],'EdgeColor','none');
    hold on
    a7=surf(x_t2,y_t2,z_t2,'FaceColor',[1 1 0],'EdgeColor','none');
    hold on
    plot3([a1_cw(1),a2_cw(1)],[a1_cw(2),a2_cw(2)],[a1_cw(3),a2_cw(3)],'r','LineWidth',3);
    input = A.q(4:7,i);
    output = quater2rotate(input);
    direction = output(2:4);
    theta = output(1);
    
    rotate(a1,direction,(theta/pi)*180,[q_w1(1),q_w1(2),q_w1(3)]);
    rotate(a2,direction,(theta/pi)*180,[q_w2(1),q_w2(2),q_w2(3)]);
    rotate(a3,direction,(theta/pi)*180,[q_x,q_y,q_z]);
      
    rotate(a4,direction,(theta/pi)*180,[q_wb1(1),q_wb1(2),q_wb1(3)]);
    rotate(a5,direction,(theta/pi)*180,[q_wb2(1),q_wb2(2),q_wb2(3)]);
    
    rotate(a6,direction,(theta/pi)*180,[q_wt1(1),q_wt1(2),q_wt1(3)]);
    rotate(a7,direction,(theta/pi)*180,[q_wt2(1),q_wt2(2),q_wt2(3)]);
 
    axis('equal');
    
    
    
    
    [X,Y] = meshgrid(range(1):range(2)-range(1):range(2),range(3):range(4)-range(3):range(4));
	Z = X.*0+Y.*0;
    surf(X,Y,Z,'FaceColor',[0 1 0],'Linestyle','-');  

    axis(range);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');


    hold off
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);





end