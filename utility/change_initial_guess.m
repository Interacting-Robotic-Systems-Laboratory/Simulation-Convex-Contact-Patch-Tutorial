function Z_new = change_initial_guess(A,Z)



  
R = rand;



Z_new = Z;

if size(A.dim,2) == 3
   Z_new(24) = R*Z(24); 
elseif size(A.dim,2) == 2
    Z_new(21) = R*Z(21); 
end


end