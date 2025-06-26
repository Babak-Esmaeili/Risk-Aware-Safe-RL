function P = P_matrix_calc(a,b,theta)

    P = [ ((cos(theta)^2)/(a^2))+((sin(theta)^2)/(b^2))       cos(theta)*sin(theta)*(1/a^2-1/b^2) 
               cos(theta)*sin(theta)*(1/a^2-1/b^2)        ((sin(theta)^2)/(a^2))+((cos(theta)^2)/(b^2)) ];

end