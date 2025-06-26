function [dist, idx] = dist_calc(x_curr,X_h)

    nVertices = size(X_h,2);

    all_dist = zeros(1,nVertices);

    for j = 1:nVertices
        
        if(j ~= nVertices)

            all_dist(j) = norm(x_curr-X_h(:,j)) + norm(x_curr-X_h(:,j+1));

        else

            all_dist(j) = norm(x_curr-X_h(:,j)) + norm(x_curr-X_h(:,1));

        end

    end

    [dist, idx] = min(all_dist);

end
