function [b,s] = cal_b_s_p(px,py,vx,vy,opt)
%cal_b_p calculates b(s(p)) given px, py in R(1,m), vx, vy in R(n,1)
% note vx(1,1) = vx(n,1); vy(1,1) = vy(n,1)

if nargin >5
    error('cal_V function: too many inputs');
end

samp_size = size(px,2);
b = zeros(samp_size,1); s = zeros(samp_size,1);
vect_size = size(vx,1)-1;

switch nargin
    case 4
        shift_vx = circshift(vx,1);
        lx = [vx(2:end,1) shift_vx(2:end,1)]';
        shift_vy = circshift(vy,1);
        ly = [vy(2:end,1) shift_vy(2:end,1)]';
        for i = 1:samp_size
            one_point = [px(1,i); py(1,i)]; % this is one point
            d = dis_p2l(one_point,lx,ly,1);
            [di, ind]= min(d); prod = 1;
            for j = 1:vect_size
                if j ~= ind
                    prod = prod*(d(j)-di)/d(j);
                end
            end 
            s(i) = 1 - prod;
            if (s(i) < 1- 10e-6) && (s(i)>10e-6)
                b(i) = lamda_func(s(i))/(lamda_func(s(i))+lamda_func(1-s(i)));
            else
                if (s(i) >= 1- 10e-6)
                b(i) = 1;
                else
                   b(i) = 0; 
                end
            end

        end
    case 5
        xi_g = opt;
        point = [px; py];
        vertices = [vx(2:end,1) vy(2:end,1)]'; % vertices in R(2,m)
        % determin with region the point is in.
        c = cell(vect_size,1);
        vertices_all = [vertices, xi_g];
        n = vect_size;
        for i = 1:n-1
            c(i,1) = {[i,i+1,n+1]};
        end
        c(n,1) = {[1,n,n+1]};
        ii = deter_Q(point,vertices_all,c); 
        
        for i = 1:samp_size
            one_point = [px(1,i); py(1,i)]; % this is one point
            poly_ver = vertices_all(:,c{ii(i)}); % should be in R(2,3).
            vx_poly_ver = poly_ver(1,:)';
            vy_poly_ver = poly_ver(2,:)';
            shift_vx_poly_ver = circshift(vx_poly_ver,1);
            lx_vx_poly_ver = [vx_poly_ver shift_vx_poly_ver]';
            shift_vy_poly_ver = circshift(vy_poly_ver,1);
            ly_vy_poly_ver = [vy_poly_ver shift_vy_poly_ver]';            
           
            d = dis_p2l(one_point,lx_vx_poly_ver,ly_vy_poly_ver,1);
            if (d(2) < 10e-6)
                b(i) = 0;
            else
                if (d(1) < 10e-6 || d(3) < 10e-6 )
                b(i) = 1;
                else
                   b(i) = tanh(5*d(2)*d(2)/(d(1)*d(3))); 
                end
            end

        end            
end

