function [Q_ind, count] = deter_Q(xi,vertices_all,c,opt)
%DETER_Q xi in R(2,p), vertices_all in R(2,m), c in cell{n,1}
%   Q-ind(i,1) being j indicates that the ith point lies in the jth polytope.  
%  Q_ind(i,1) = 0 means xi does not belong to any polytope.   
if nargin >4
    error('deter_Q function: too many inputs');
end


px = xi(1,:); py = xi(2,:);

p = size(xi,2); Q_ind = zeros(p,1); count = zeros(p,1);
for i = 1:size(c,1)
    vertices = vertices_all(:,c{i});
    vx = vertices(1,:);
    vy = vertices(2,:);
    [in,on] = inpolygon(px,py,vx,vy);
    Q_ind(in == 1) = i;
    Q_ind(on == 1) = i;
    count(in == 1) = count(in == 1) + 1;
end

switch nargin
    case 3
        for i = 1: p
            if Q_ind(i,1)  == 0
                error('deter_Q error: one point does not belong to any polytope');
            end

            if count(i,1) > 1
                warning('deter_Q error: one point seems belonging to at least 2 polytopes');
            end
        end
    case 4
        for i = 1: p
            if count(i,1) > 1
                warning('deter_Q error: one point seems belonging to at least 2 polytopes');
            end
        end
end

end