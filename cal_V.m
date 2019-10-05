function V = cal_V(point,f_ex,f_ey,vertices,opt)
% cal V given point in R(2,m),f_ex,f_ey in R(2,1),vertices in R(2,n);

if nargin >5
    error('cal_V function: too many inputs');
end

n = size(vertices,2); % number of edges

px = point(1,:);
py = point(2,:);
vx = vertices(1,:); vx(1,n+1) = vx(1,1); vx = vx';
vy = vertices(2,:); vy(1,n+1) = vy(1,1); vy = vy';

m = size(point,2);
V = zeros(2,m);

switch nargin
    case 4
        V_c = cal_Vc(point,f_ex,f_ey,vertices);
        V_f = cal_Vf(point,f_ex,f_ey,vertices); 
        [b,~] = cal_b_s_p(px,py,vx,vy);
        for i = 1:m
            V(:,i) = (b(1,i)*V_c(:,i) + (1-b(1,i))*V_f(:,i))/norm(b(1,i)*V_c(:,i) + (1-b(1,i))*V_f(:,i));
        end
    case 5
        xi_g = opt;
        V_c = cal_Vc(point,f_ex,f_ey,vertices, xi_g); % not normalize, maximize norm = 1
        V_f = cal_Vf(point,f_ex,f_ey,vertices, xi_g); % normalize
        [b,~] = cal_b_s_p(px,py,vx,vy,xi_g); %change
        for i = 1:m
            V(:,i) = (b(1,i)*V_c(:,i) + (1-b(1,i))*V_f(:,i));
        end
end

end
