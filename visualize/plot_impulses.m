function plot_impulses(Impulses)
P_x = Impulses(:,1);
P_y = Impulses(:,2);
P_z = Impulses(:,3);

P_xt = Impulses(:,4);
P_yt = Impulses(:,5);
P_zt = Impulses(:,6);

T =0.01:0.01:5;
figure
subplot(2,1,1)
plot(T',P_x,'r','LineWidth',1);
hold on
plot(T',P_y,'b-.','LineWidth',1);
hold on
plot(T',P_z,'k--','LineWidth',1);
hold on
pgon1 = polyshape([0 0 0.65 0.65],[10 -10 -10 10]);
plot(pgon1,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon2 = polyshape([1.5 1.5 2.15 2.15],[10 -10 -10 10]);
plot(pgon2,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon3 = polyshape([3.25 3.25 4.7 4.7],[10 -10 -10 10]);
plot(pgon3,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
text(0.25,5,'T1')
text(1,5,'T2')
text(1.75,5,'T3')
text(2.6,5,'T4')
text(3.9,5,'T5')
text(4.75,5,'T6')
axis([0 5 0 7]);
legend('p_x','p_y','p_z');
xlabel('Time (s)');
ylabel('Impulse (Ns)');
grid on
subplot(2,1,2)
plot(T',P_xt,'r','LineWidth',1);
hold on
plot(T',P_yt,'b-.','LineWidth',1);
hold on
plot(T',P_zt,'k--','LineWidth',1);
hold on
pgon1 = polyshape([0 0 0.65 0.65],[10 -10 -10 10]);
plot(pgon1,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon2 = polyshape([1.5 1.5 2.15 2.15],[10 -10 -10 10]);
plot(pgon2,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
pgon3 = polyshape([3.25 3.25 4.7 4.7],[10 -10 -10 10]);
plot(pgon3,'FaceColor','k','FaceAlpha',0.1,'LineStyle','--')
hold on
text(0.25,-0.6,'T1')
text(1,-0.6,'T2')
text(1.75,-0.6,'T3')
text(2.6,-0.6,'T4')
text(3.9,-0.6,'T5')
text(4.75,-0.6,'T6')
legend('p_{x\tau}','p_{y\tau}','p_{z\tau}');
xlabel('Time (s)');
ylabel('Angular mpulse (Nms)');
axis([0 5 -1 1]);
grid on
export_fig manipulation_impulses.png -m3.5 -transparent
end