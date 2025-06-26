function dist = find_distance_interpolation(x_curr,vertices,idx)

%     syms x1 x2;

    vert_1 = vertices(:,idx);
    vert_2 = vertices(:,idx+1);

    m = (vert_2(2)-vert_1(2))/(vert_2(1)-vert_1(1));
    c = -m*vert_1(1) + vert_1(2);
    
    vert_1_curr = [0 0]';
    vert_2_curr = x_curr;

    m_curr = (vert_2_curr(2)-vert_1_curr(2))/(vert_2_curr(1)-vert_1_curr(1));
    c_curr = -m_curr*vert_1_curr(1) + vert_1_curr(2);

%     eq_1 = m*x1+c-x2 == 0;
%     eq_2 = m_curr*x1+c_curr-x2 == 0;
% 
%     sol = solve(eq_1,eq_2,[x1 x2]);
% 
%     dist = norm(double([sol.x1 sol.x2]'));

    x1 = det([-c -1;-c_curr -1])/det([m -1;m_curr -1]);
    x2 = det([m -c;m_curr -c_curr])/det([m -1;m_curr -1]);

    dist = norm([x1 x2]');

end