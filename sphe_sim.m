function  [x_t,t] = sphe_sim(x_0,x_g,a,svertices_all,c,seg_ind,time,opt)
%SPHE_SIM simulation of the spherical dynamics
% 

if nargin >8
    error('sphe_sim function: too many inputs');
end

dt = 0.001;
steps = time/dt;
x_t = ones(3, steps);
x_t(:,1) = x_0;

switch nargin
    case 7
        for i = 1:steps-1
            u = fb_cons_one_chart(x_t(:,i),a,svertices_all,c,seg_ind, x_g);
            Pi = cross_vec(x_t(:,i));
            dot_x_t = Pi*u;
            x_t(:,i +1) = x_t(:,i) + dt*dot_x_t;
            x_t(:,i +1) = x_t(:,i +1)/norm(x_t(:,i +1));
        end
    case 8
        for i = 1:steps-1
            svertices_chart = opt;
            u = fb_cons_two_charts(x_t(:,i),a,svertices_all,svertices_chart,c,seg_ind, x_g);
            Pi = cross_vec(x_t(:,i));
            dot_x_t = Pi*u;
            x_t(:,i +1) = x_t(:,i) + dt*dot_x_t;
            x_t(:,i +1) = x_t(:,i +1)/norm(x_t(:,i +1));
        end
end
t = steps*dt;
end

