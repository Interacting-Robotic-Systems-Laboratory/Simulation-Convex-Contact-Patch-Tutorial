function Z_new = change_initial_guess(A,Z)



  
R = rand;



Z_new = Z;

if strcmp(A.object,'particle') == 1
    Z_new(7) = R*Z(7);
    return
end

if strcmp(A.object,'cube') == 1
   Z_new(24) = R*Z(24); 
elseif strcmp(A.object,'cylinder') == 1
   Z_new(21) = R*Z(21); 
elseif strcmp(A.object,'ellipse') == 1
   Z_new(19) = R*Z(19);
end


end