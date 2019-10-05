function [u,dot_sp] = fb_cons_one_chart(x,a,svertices_all,c,seg_ind, x_g)  
%% for one state x on the sphere, calculate the control input  u(x)

    xi = phi_a(x,a);
    vertices_all = phi_a(svertices_all,a);
    Q_ind = deter_Q(xi,vertices_all,c);
    svertices = svertices_all(:,c{Q_ind});
    vertices = phi_a(svertices,a);
    Q_ind_g = 0;
    if a'*x_g >0
        xi_g = phi_a(x_g,a);
        Q_ind_g = deter_Q(xi_g,vertices_all,c);
    end

    if Q_ind ~= Q_ind_g % intermediate polytope case
        f_ex = vertices_all(1,seg_ind{Q_ind,1})';
        f_ey = vertices_all(2,seg_ind{Q_ind,1})';
        v = cal_V(xi,f_ex, f_ey, vertices);
        [u, ~] = cal_v2u(v,x,a);
        Pi = cross_vec(x); % dynamics
        u = u/norm(Pi*u);
        dot_sp = Pi*u;
    else % goal polytope case
        f_ex = vertices_all(1, seg_ind{Q_ind,1})';
        f_ey = vertices_all(2, seg_ind{Q_ind,1})';
        xi_g = phi_a(x_g,a);
        v = cal_V(xi,f_ex, f_ey, vertices,xi_g); 
        [u,~] = cal_v2u(v,x,a);
        Pi = cross_vec(x); % dynamics        
        u_temp = 10*norm(x_g - x);
        b = tanh(u_temp); % Ranging from 0-1
       if (b >1 || b < 0)
           warning('something wrong with b');
       end
        u = b* u/norm(Pi*u);
        dot_sp = Pi*u; 
    end
end