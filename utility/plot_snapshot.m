function  plot_snapshot(A)

N = A.N;
r = A.dim(1)*A.unit; %radius of the bell
height = A.dim(2)*A.unit;
r1 = 0.5*r; %radius of the bar
heightr = 0.2*height;

a11 = [0;0;height/2;1];
a12 = [0;0;height/2-heightr;1];

a21 = [0;0;-height/2+heightr;1];
a22 = [0;0;-height/2;1];
figure 
for i = 100:10:N
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
    
    a11_w = H*a11;
    a12_w = H*a12;
    a21_w = H*a21;
    a22_w = H*a22;
    
    
    
    plot(A.q(1,i),A.q(2,i),'ko','MarkerSize',5,'LineWidth',1);
    hold on
    plot(A.z(7,i),A.z(8,i),'kx','MarkerSize',6,'LineWidth',2);
    hold on
    plot([a12_w(1),a21_w(1)],[a12_w(2),a21_w(2)],'r-.','LineWidth',1.5);
    hold on
    plot([a11_w(1),a12_w(1)],[a11_w(2),a12_w(2)],'k','LineWidth',1.5);
    hold on
    plot([a21_w(1),a22_w(1)],[a21_w(2),a22_w(2)],'k','LineWidth',1.5);
    hold on
    plot(A.q(1,i),A.q(2,i),'ko','MarkerSize',5,'LineWidth',1);
    hold on
    plot(A.z(7,i),A.z(8,i),'kx','MarkerSize',6,'LineWidth',2);
    hold on
    %plot(A.z(7,100:N),A.z(8,100:N),'k');
    axis square
    axis([-0.2 1.3 0 1.5]);
end
legend({'CM','ECP'},'FontSize',12);
xlabel('x coordinate (m)','FontSize',12);
ylabel('y coordinate (m)','FontSize',12);