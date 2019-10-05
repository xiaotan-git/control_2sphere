function [u, dot_sp] = cal_trans_chart( x,a1,a2,f_ex_share,f_ey_share, f_ex, f_ey, vertices_a, vertices_b, svertices_share)
%CAL_TRANS_CHART Calculate u for one state x that lies in Reg1,a or Reg2,b
% x the state; a1 in R(3,1) the phi_a1; a2 in R(3,1) the phi_a2; 
% f_ex in R(1,1), f_ey in R(),vertices,seg,svertices_a2
% vertices_a represents the coordination of polytope 1 in chart a
% vertices_b represents the coordination of polytope 2 in chart b

% from cell 2 TO cell 1, illustrated in Fig.7. or Algorithm 3

xi1 = phi_a(x,a1); 

v_fb = phi_a(svertices_share,a2);
xi_ex_b = (v_fb(:,1) + v_fb(:,2))/2; % be a central point on the exit face

f_ex_b = v_fb(1,:)';
f_ey_b = v_fb(2,:)';

V_f_2b = cal_Vf(xi_ex_b,f_ex_b,f_ey_b,vertices_b); 
V_f_1b = V_f_2b; %cal V_f_1b;


V_c_1a = cal_Vc(xi1,f_ex,f_ey,vertices_a); % NB. the input f_ex, f_ey and f_ex_share, f_ey_share are different.
u = cal_v2u(V_c_1a,x,a1); % from phi_a_inv(xx,a1) --> 2-sphere
V_c_1b = cal_u2v(u,x,a2); %cal V_c_1b; from 2-sphere  -->  phi_a(xx,a2)

px = xi1(1,:); py = xi1(2,:);
vx = vertices_a(1,:);  vx(1,size(vertices_a,2)+1) = vx(1,1); vx = vx';
vy = vertices_a(2,:);  vy(1,size(vertices_a,2)+1) = vy(1,1); vy = vy';
b = cal_b_s_p(px,py,vx,vy); %cal b

v = (b*V_c_1b + (1-b)*V_f_1b)/norm(b*V_c_1b + (1-b)*V_f_1b); %cal V_1b


X = ['b is  ', num2str(b),'   \n V_c_1b is  ', num2str(reshape(V_c_1b,1,2)), '    \n V_f_1b is   ',num2str(reshape(V_f_1b,1,2)) ];
disp(X)

u = cal_v2u(v,x,a2);
Pi = cross_vec(x);
u = u/norm(Pi*u);
dot_sp = Pi*u; 
end
% sv_1 = svertices_a1(:,Q_ind1_share);
% sv_2 = svertices_a2(:,Q_ind2_share);
% v_1a = phi_a(sv_1,a1);
% v_2b = phi_a(sv_2,a2);
% % v_GVD_1a = GVD(v_1a,share_edge); % work out the voronoi points in 2,b
% sv_GVD_1 = phi_a_inv(v_GVD_1a,a1);
% v_Reg1b = phi_a(sv_GVD_1,a2); 
% % p = f(v_Reg1b, v_2b);  % such that the three conditions on the face vector field are satisfied


% if flag == 1      % cell 1, the latter cell, suppose cell 1 is not the goal cell
%     V_c_1a = cal_Vc(xi1,f_ex,f_ey,vertices); 
%     u = cal_v2u(V_c_1a,x,a2);
%     V_c_1b = cal_u2v(u,x,a1); %cal V_c_1b; 
%     V_f_1b = (p - xi2)/norm(p - xi2); %cal V_f_1b;
%     px = xi1(1,:); py = xi1(2,:);
%     vx = v_1a(1,:);  vx(1,size(v_1a,2)+1) = vx(1,1); vx = vx';
%     vy = v_1a(2,:);  vy(1,size(v_1a,2)+1) = vy(1,1); vy = vy';
%     b = cal_b_s_p(px,py,vx,vy);
%     v = (b*V_c_1b + (1-b)*V_f_1b)/norm(b*V_c_1b + (1-b)*V_f_1b); %cal V_1b
%     u = cal_v2u(v,x,a2);
%     Pi = cross_vec(x);
%     u = u/norm(Pi*u);
%     dot_sp = Pi*u; 
% 
% else             % cell 2, the former cell
%     V_c_2b = cal_Vc(xi2,f_ex,f_ey,vertices); %cal V_c_2b; 
%     px = xi2(1,:); py = xi2(2,:);
%     vx = v_2b(1,:);  vx(1,size(v_2b,2)+1) = vx(1,1); vx = vx';
%     vy = v_2b(2,:);  vy(1,size(v_2b,2)+1) = vy(1,1); vy = vy';
%     b = cal_b_s_p(px,py,vx,vy);
%     v = (b*V_c_2b + (1-b)*V_f_2b)/norm(b*V_c_2b + (1-b)*V_f_2b); %cal V_2b
%     u = cal_v2u(v,x,a2);
%     Pi = cross_vec(x);
%     u = u/norm(Pi*u);
%     dot_sp = Pi*u;
% end

