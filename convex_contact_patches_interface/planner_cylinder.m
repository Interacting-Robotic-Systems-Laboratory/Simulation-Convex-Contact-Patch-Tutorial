function A = planner_cylinder(A)
N =250;
A.N =N; % number of time steps in simulation
A.Impulses = zeros(N,6);
P=150*0.01;
A.Total = zeros(N,1);
for i = 1:3
    A.Total(i) = P/2;
end

for i = 100:N
    A.Total(i) = P/5;
end



end