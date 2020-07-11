function plot_configuration(z,q)
% v_x = z(1,:);
% v_y = z(2,:);
q_z = q(3,:);
v_z = z(3,:);
T =0.01:0.01:5;
figure
subplot(2,1,1)
% plot(T',q_x,'r','LineWidth',1);
% hold on
% plot(T',q_y,'b-.','LineWidth',1);
% hold on
plot(T',v_z,'k','LineWidth',1);
hold on
plot(0.65,0,'r.','MarkerSize',10);
plot(2.15,0,'r.','MarkerSize',10);
plot(4.71,0,'r.','MarkerSize',10);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
hold on
pgon1 = polyshape([0 0 0.65 0.65],[0.5 -1 -1 0.5]);
plot(pgon1,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon2 = polyshape([1.5 1.5 2.15 2.15],[0.5 -1 -1 0.5]);
plot(pgon2,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon3 = polyshape([3.25 3.25 4.7 4.7],[0.5 -1 -1 0.5]);
plot(pgon3,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
text(0.25,-0.75,'T1')
text(1,-0.75,'T2')
text(1.75,-0.75,'T3')
text(2.6,-0.75,'T4')
text(3.9,-0.75,'T5')
text(4.75,-0.75,'T6')
hold on
plot(0.65,0,'r.','MarkerSize',10);
plot(2.15,0,'r.','MarkerSize',10);
plot(4.71,0,'r.','MarkerSize',10);
grid on
axis([0 5 -1 0.5])
legend('v_z(t)','Impact events','location','Northwest');


subplot(2,1,2)
plot(T',q_z,'k','LineWidth',1);
hold on
plot(T',0.15*ones(500,1),'k--','LineWidth',1);
hold on
T =0:0.5:5;
plot(T',0.1*ones(11,1),'k:','LineWidth',1);
hold on
pgon1 = polyshape([0 0 0.65 0.65],[0.5 -1 -1 0.5]);
plot(pgon1,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon2 = polyshape([1.5 1.5 2.15 2.15],[0.5 -1 -1 0.5]);
plot(pgon2,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon3 = polyshape([3.25 3.25 4.7 4.7],[0.5 -1 -1 0.5]);
plot(pgon3,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
text(0.25,0.12,'T1')
text(1,0.12,'T2')
text(1.75,0.12,'T3')
text(2.6,0.12,'T4')
text(3.9,0.12,'T5')
text(4.75,0.12,'T6')
hold on
legend('q_z(t)','Surface','Line','location','Northwest');
axis([0 5 0.09 0.21])
xlabel('Time (s)');
ylabel('Position (m)');
yticks([0.1,0.12,0.14,0.15,0.16,0.18,0.2]);
yticklabels({'\color{red} 0.1','0.12','0.14','\color{red} 0.15','0.16','0.18','0.2'});

grid on
export_fig manipulation_trajectory.png -m3.5 -transparent
end