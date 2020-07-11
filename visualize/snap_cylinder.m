function  snap_cylinder(A)

r = A.dim(1)*A.unit;
h = A.dim(2)*A.unit;
q_t1 = [0;0;h/2;1];
q_t2 = [0;0;-h/2;1];
a1_c  = [0;r;-h/2;1];
a2_c  = [0;-r;-h/2;1];
a3_c  = [0;r;h/2;1];
a4_c  = [0;-r;h/2;1];
range =[-0.2 1.9 -1.35 0.75 0 1.5];

    
    
    for i = 500:500
        
        axis(range);
        
         % frame
        
        
        xr = 0.2;
        yr = -0.6;
        zr = 0.2;
        L = 0.3;
        p1 = [xr,yr,zr;xr,yr,zr;xr,yr,zr];
        p2 = [xr+L,yr,zr;xr,yr+L,zr;xr,yr,zr+L];
        w=10*L;
        arrow3(p1,p2,'k',w,3*w);
        axis equal;
        hold on
    % applied force and torque
   
        if (i>1)&&(i<=30)
            
            p1 = [xr,yr,zr];
            p2 = [xr,yr+L,zr]';
            arrow3(p1,p2','r',w,3*w);

            

            rr = L/3;
            th = 0:pi/50:3*pi/2;
            yunit = rr*cos(th) + yr;
            zunit = rr*sin(th) + zr;
            xunit = 0*sin(th) + xr+L;
            plot3(xunit,yunit,zunit,'r');

            p1 = [xr+L,yr+r,zr];
            p2 = p1-[0,0,r];
            arrow3(p1,p2,'r',w,3*w);

        elseif (71<=i)&&(i<=90)
            p1 = [xr,yr,zr];
            p2 = [xr+L,yr,zr]';
            arrow3(p1,p2','r',w,3*w);   
        elseif (150<=i)&&(i<=160)
            p1 = [xr,yr,zr];
            p2 = [xr,yr,zr+L]';
            arrow3(p1,p2','r',w,3*w);
            
            rr = L/3;
            th = pi/2:pi/50:2*pi;
            yunit = rr*cos(th) + yr;
            zunit = rr*sin(th) + zr;
            xunit = 0*sin(th) + xr+L;
            plot3(xunit,yunit,zunit,'r');

            p1 = [xr+L,yr+r,zr];
            p2 = p1+[0,0,r];
            arrow3(p1,p2,'r',w,3*w);
        elseif (161<=i)&&(i<=210)
            p1 = [xr,yr,zr];
            p2 = [xr,yr,zr+L]';
            arrow3(p1,p2','r',w,3*w);
            
            rr = L/3;
            th = pi/2:pi/50:2*pi;
            yunit = rr*cos(th) + yr;
            zunit = rr*sin(th) + zr;
            xunit = 0*sin(th) + xr+L;
            plot3(xunit,yunit,zunit,'r');

            p1 = [xr+L,yr+r,zr];
            p2 = p1+[0,0,r];
            arrow3(p1,p2,'r',w,3*w);
        elseif (221<=i)&&(i<=231)
            rr = L/3;
            th = pi/2:pi/50:2*pi;
            yunit = rr*cos(th) + yr;
            zunit = 0*sin(th) + zr + L;
            xunit = rr*sin(th) + xr;
            plot3(xunit,yunit,zunit,'r');

            p1 = [xr,yr+r,zr+L];
            p2 = p1+[r,0,0];
            arrow3(p1,p2,'r',w,3*w);
        elseif (276<=i)&&(i<=320)
            p1 = [xr,yr,zr];
            p2 = [xr,yr+L,zr]';
            arrow3(p1,p2','r',w,3*w);
        elseif (325<=i)&&(i<=342)
            p1 = [xr,yr,zr];
            p2 = [xr,yr,zr+L]';
            arrow3(p1,p2','r',w,3*w);
            
            rr = L/3;
            th = pi/2:pi/50:2*pi;
            yunit = 0*cos(th) + yr+L;
            zunit = rr*sin(th) + zr;
            xunit = rr*cos(th) + xr;
            plot3(xunit,yunit,zunit,'r');

            p1 = [xr+r,yr+L,zr];
            p2 = p1+[0,0,r];
            arrow3(p1,p2,'r',w,3*w);
        end
    
    
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
    
    q_wt1 = H*q_t1;
    q_wt2 = H*q_t2;
    a1_cw = H*a1_c;
    a2_cw = H*a2_c;
    a3_cw = H*a3_c;
    a4_cw = H*a4_c;
    
    [x,y,z] = cylinder_mod([q_x,q_y,q_z],r,h );
    [x_t1, y_t1, z_t1] = ellipsoid(q_wt1(1),q_wt1(2),q_wt1(3),r,r,0,30);
    [x_t2, y_t2, z_t2] = ellipsoid(q_wt2(1),q_wt2(2),q_wt2(3),r,r,0,30);
    
    a = surf(x,y,z,'FaceColor',[1 1 0]);
    hold on
    a6=surf(x_t1,y_t1,z_t1,'FaceColor',[1 1 0],'EdgeColor','none');
    hold on
    a7=surf(x_t2,y_t2,z_t2,'FaceColor',[1 1 0],'EdgeColor','none');
    hold on
    plot3([a1_cw(1),a2_cw(1)],[a1_cw(2),a2_cw(2)],[a1_cw(3),a2_cw(3)],'r','LineWidth',1);
    hold on
    plot3([a3_cw(1),a4_cw(1)],[a3_cw(2),a4_cw(2)],[a3_cw(3),a4_cw(3)],'r','LineWidth',1);
    input = A.q(4:7,i);
    output = quater2rotate(input);
    direction = output(2:4);
    theta = output(1);
    rotate(a,direction,(theta/pi)*180,[q_x,q_y,q_z]); 
     
    rotate(a6,direction,(theta/pi)*180,[q_wt1(1),q_wt1(2),q_wt1(3)]);
    rotate(a7,direction,(theta/pi)*180,[q_wt2(1),q_wt2(2),q_wt2(3)]);
    hold on
    % start and goal
    rr = 0.15;
    th = 0:pi/50:2*pi;
    yunit = rr*cos(th);
    xunit = rr*sin(th);
    zunit = 0*th;
    plot3(xunit,yunit,zunit,':k','LineWidth',1);
    axis('equal');
    axis(range);
    hold on
    xunit = rr*sin(th)+0.65;
    yunit = rr*cos(th)-1.2;
    plot3(xunit,yunit,zunit,':k','LineWidth',1);
    
    % the obstacle
    H = eye(4);
    [xo1,yo1,zo1] = mashgrid_cuboid(0.5,0.6,0.05,H,[0.9;0.5;0.025]);
    surf(xo1,yo1,zo1,'FaceColor',[0 1 0]);
    
    [xo2,yo2,zo2] = mashgrid_cuboid(0.75,0.75,0.05,H,[0.25;-1.25/2;0.025]);
    surf(xo2,yo2,zo2,'FaceColor',[0 1 0]);
    axis('equal');
    axis(range);
     
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    grid on
    view([-140,30]);
    hold off    
    export_fig manipulation_8.png -m3.5 -transparent
    end

end




