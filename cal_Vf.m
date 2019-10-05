function V_f = cal_Vf(point,f_ex,f_ey,vertices,opt)
%calculate V_f, given m ps in R(2,m),f_ex in R(2,1),f_ey in R(2,1),  
% n vertices in R(2,n)


if nargin >5
    error('cal_V function: too many inputs');
end

%% get the index of the exit face if it is an intermediate cell
f_v1 = [f_ex(1,1); f_ey(1,1) ]; % point
f_v2 = [f_ex(2,1); f_ey(2,1) ]; % point
[exis_flag1,exis_ind1] = ismembertol(f_v1',vertices','Byrows',true);
[exis_flag2,exis_ind2] = ismembertol(f_v2',vertices','Byrows',true);
m = size(point,2);
n = size(vertices,2); % number of edges

V_c = zeros(2,m); % needs to be deleted

% calculate n_fi of segment AB in R(2,n); CAB
n_fi = zeros(2,n);
for i = 1:n
    if (i ~= n)&&(i ~= 1)
        A = vertices(:,i);
        B = vertices(:,i+1);
        C = vertices(:,i-1);
        v_AB = B-A;
        v_AC = C -A;
        n_fi_temp = [-v_AB(2,1); v_AB(1,1)];
        n_fi(:,i) = sign(v_AC'*n_fi_temp)*n_fi_temp/norm(n_fi_temp);
    else
        if i == 1
            A = vertices(:,i);
            B = vertices(:,i+1);
            C = vertices(:,end);
            v_AB = B-A;
            v_AC = C -A;
            n_fi_temp = [-v_AB(2,1); v_AB(1,1)];
            n_fi(:,i) = sign(v_AC'*n_fi_temp)*n_fi_temp/norm(n_fi_temp);
        else
            A = vertices(:,i);
            B = vertices(:,1);
            C = vertices(:,i-1);
            v_AB = B-A;
            v_AC = C -A;
            n_fi_temp = [-v_AB(2,1); v_AB(1,1)];
            n_fi(:,i) = sign(v_AC'*n_fi_temp)*n_fi_temp/norm(n_fi_temp);
        end
    end
end


V_f = zeros(2,m);
        
switch nargin
    case 4
        if ~(exis_flag1&&exis_flag2)
            error('cal_Vc function: did not find the vertice of the exit edge');
        end

        if ~(( abs(exis_ind1 - exis_ind2)== 1)||((exis_ind2 == 1)&&(exis_ind1 == n)||((exis_ind1 == 1)&&(exis_ind2 == n))))
                error('cal_Vc function: the vertice of the exit edge do not match!')
        end
        
        % calculate alpha_i in R(m,n): the mth point in region n
        exis_ind_min = min(exis_ind1,exis_ind2);
        exis_ind_max = max(exis_ind1,exis_ind2);
        if (exis_ind_min == 1)&&(exis_ind_max == n)
            n_fi(:,n) = - n_fi(:,n);
        else
            n_fi(:,exis_ind_min) = - n_fi(:,exis_ind_min);
        end

        % determin which region the point is in.
        d = zeros(m,n);
        vx = vertices(1,:); vx(1,n+1) = vx(1,1); vx = vx';
        vy = vertices(2,:); vy(1,n+1) = vy(1,1); vy = vy';
        shift_vx = circshift(vx,1);
        lx = [shift_vx(2:end,1) vx(2:end,1)]';
        shift_vy = circshift(vy,1);
        ly = [shift_vy(2:end,1) vy(2:end,1)]';
        for i =1:m
            d(i,:) =  dis_p2l(point(:,i),lx,ly,1);
        end
    
        [~,ii] = min(d,[],2); % minimize each column, with index ii in R(m,1)
          
        for i = 1:m
             V_f(:,i) = n_fi(:,ii(i));
        end
    case 5
        xi_g = opt;
        % determin with region the point is in
        
        c = cell(n,1);
        vertices_all = [vertices, xi_g];
        for i = 1:n-1
            c(i,1) = {[i,i+1,n+1]};
        end
        c(n,1) = {[1,n,n+1]};
        
        ii = deter_Q(point,vertices_all,c);
        for i = 1:m
             V_f(:,i) = n_fi(:,ii(i));
        end
%         V_f = zeros(size(V_f));

end

end
