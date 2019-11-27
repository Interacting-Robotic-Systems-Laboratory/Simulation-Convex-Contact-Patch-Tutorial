function plot_impulse(A)
Input = A.input;
px = Input(1,:);
py = Input(2,:);
pxt = Input(4,:);
pzt = Input(6,:);
N=A.N;
T=1:N;
figure
subplot(2,1,1)
plot(T*A.h,px(1:N),'k','LineWidth',1.5);
hold on
plot(T*A.h,py(1:N),'b--','LineWidth',1);
hold on
legend({'p_x','p_y'},'FontSize',12);
xlabel('Time (s)','FontSize',12);
ylabel('Impulse (Ns)','FontSize',12);
subplot(2,1,2)
plot(T*A.h,pxt(1:N),'k','LineWidth',1.5);
hold on
plot(T*A.h,pzt(1:N),'b--','LineWidth',1);
hold on
legend({'p_{x\tau}','p_{z\tau}'},'FontSize',12);
xlabel('Time (s)','FontSize',12);
ylabel('Angular Impulse (Nms)','FontSize',12);
end

