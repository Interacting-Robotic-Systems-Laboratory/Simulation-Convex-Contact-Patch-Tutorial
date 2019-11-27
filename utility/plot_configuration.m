function plot_configuration(A)
q_x = A.q(1,:);
q_y = A.q(2,:);
q_z = A.q(3,:);

T = 1:A.N;
plot(T*A.h,q_x,'b','LineWidth',2);
hold on
plot(T*A.h,q_y,'r--','LineWidth',2);
hold on
plot(T*A.h,q_z,'k-.','LineWidth',2);


legend({'q_x','q_y','q_z'},'FontSize',12);
xlabel('Time (s)','FontSize',12);
ylabel('Position (m)','FontSize',12);


end