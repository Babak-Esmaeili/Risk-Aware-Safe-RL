function F = find_F(vertices)

    nStates = size(vertices,1);
    nVertices = size(vertices,2);

    m = zeros(1,nVertices);
    c = zeros(1,nVertices);
    
    F = zeros(nVertices,nStates);

    for i = 1:nVertices

        if(i ~= nVertices)

            m(i) = (vertices(2,i+1)-vertices(2,i))/(vertices(1,i+1)-vertices(1,i));
    
            c(i) = -m(i)*vertices(1,i) + vertices(2,i);
    
            F(i,:) = [-m(i)/c(i) 1/c(i)];

        else

             m(i) = (vertices(2,1)-vertices(2,i))/(vertices(1,1)-vertices(1,i));
             c(i) = -m(i)*vertices(1,i) + vertices(2,i);

             F(i,:) = [-m(i)/c(i) 1/c(i)];

        end

%         if(c(i)<=0)
% 
%             coeff = 1;
% 
%         else
% 
%             coeff = -1;
%         
%         end
% 
%         F(i,:) = coeff*F(i,:);

    end

end
