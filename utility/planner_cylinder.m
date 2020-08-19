function A = planner_cylinder(A)
N =500;
A.N =N; % number of time steps in simulation
A.Impulses = zeros(N,6);
P=300*0.01;
heg = A.dim(2);

for i = 1:30
        A.Impulses(i,2) = P;
        A.Impulses(i,4) = -P*heg/2;    
end

for i = 71:90
        A.Impulses(i,1) = P;          
end


for i = 150:160
        A.Impulses(i,3) = 2*P; 
        A.Impulses(i,4) = 2*P*heg/2; 
end

for i = 161:210
        A.Impulses(i,3) = 1*P; 
        A.Impulses(i,4) = 1*P*heg/2; 
end


for i = 221:231
        A.Impulses(i,6) = -1*P*heg/2; 
end

for i = 276:320
        A.Impulses(i,2) = P; 
end

for i = 325:342
        A.Impulses(i,3) = 2*P; 
        A.Impulses(i,5) = -2*P*heg/2; 
end

end