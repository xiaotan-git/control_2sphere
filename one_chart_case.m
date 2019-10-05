%%%%%%%%%%%%%%  spherical polytopes in one chart.
%%%%%%%%%%%%%%  Both state trajectories and vector field
%%%%%%%%%%%%%% illustrations are included.
clc; clear;
%% initial data
% suppose we just know the feasible spherical polytopes expressed by
% svertices_all and c
% the sample points should be strictly inside the polygons. Needs
% change here. Right now it is working but should be impoved.

% v = [-4 -4 -4 -3 -2.5 -2.5 -2 0 1  1 2  2    3   3  4  4    4  1   2;
%       4  1 -4  1 -1.5 -2.5  3 0 2 -1 2 -1 -1.5 -2.5 4 -1.5 -4 -1.5 -1.5];% 2*n
% v = 0.3*v;
% c = {[1,7,9,11,15],[1,2,4,7],[2,3,6,5,4],[4,5,18,10,8],[7,8,10,9],[10,18,19,12],...
%      [11,12,19,13],[15,11,13,16],[6,3,17,14],[13,14,17,16]}';
% a = [1, -1, 0]'; a = a/norm(a);
% svertices_all = phi_a_inv(v,a); 
% 
% 
% initial_state = [0, 0.58, 0.67 , 0.4;
%                  -1, -0.8, -0.6 -0.34;
%                  -0.3, -0.15, -0.48,-0.8];
% initial_state = normc(initial_state); % normalize the columns
% x_g = [1 0 0.3]'/norm([1 0 0.3]);
% 
% 
% seg_ind ={[1,15],[1,7],[2,4],[10,8],[7,9],[10,18],...
%      [11,13],[15,11],[6,3],[13,16]}';
% time = 4;
% 
% sphe_sample_spoint_temp = [];
% for i = 1:size(c,1)
%     spherical_v = svertices_all(:,c{i,1});
%     v = phi_a(spherical_v,a);
%     vx = [v(1,:) v(1,1)];
%     vy = [v(2,:) v(2,1)];
%     [px, py] = sampl_poly(vx,vy);
%     sample_point = [px; py]; 
%     sphe_sample_spoint_temp= [ sphe_sample_spoint_temp, phi_a_inv(sample_point,a)]; % 3*m
% end
% 
% sample_spoint_temp = phi_a(sphe_sample_spoint_temp,a);
% vertices_all = phi_a(svertices_all,a);
% [~,count] = deter_Q(sample_spoint_temp,vertices_all,c);
% sphe_sample_spoint = sphe_sample_spoint_temp(:,count == 1);

load('initial_data.mat');
% save('initial_data.mat', 'a','c','initial_state','seg_ind','sphe_sample_spoint',...
%     'svertices_all','vertices_all','x_g','time');

%%  state trajectories
% only one point x in R(3,1)is allowed for now.
x_0 = initial_state(:,1);
x_1 = initial_state(:,2);
x_2 = initial_state(:,3);
x_3 = initial_state(:,4);
u = fb_cons_one_chart(x_0,a,svertices_all,c,seg_ind, x_g); 

% simulation
[x_t,t] = sphe_sim(x_3,x_g,a,svertices_all,c,seg_ind,time);

% plotting
figure
sphere
axis equal
for i = 1:length(c) 
svertices = svertices_all(:,c{i});
draw_spolygon(svertices);
end
axis off
plot3( x_0(1,1), x_0(2,1), x_0(3,1),'mo','linewidth',6);
plot3( x_1(1,1), x_1(2,1), x_1(3,1),'mo','linewidth',6);
plot3( x_2(1,1), x_2(2,1), x_2(3,1),'mo','linewidth',6);
plot3( x_3(1,1), x_3(2,1), x_3(3,1),'mo','linewidth',6);
plot3( x_g(1,1), x_g(2,1), x_g(3,1),'rx','linewidth',12);
plot3(x_t(1,:), x_t(2,:), x_t(3,:),'b-');


%% vector field illustration
% calculate dot_x for all sphe_sample_spoint
spherical_vec = zeros(3, size(sphe_sample_spoint,2));
for i = 1 : size(sphe_sample_spoint,2)
    spherical_one_point =  sphe_sample_spoint(:,i);
    [ ~ , dot_sp] = fb_cons_one_chart(spherical_one_point,a,svertices_all,c,seg_ind, x_g);
    spherical_vec(:,i) = dot_sp;
end

% plotting
figure
sphere
axis equal
for i = 1:length(c) 
svertices = svertices_all(:,c{i});
draw_spolygon(svertices);
end
axis off
q2 = quiver3(sphe_sample_spoint(1,:),sphe_sample_spoint(2,:),...
    sphe_sample_spoint(3,:),spherical_vec(1,:),spherical_vec(2,:),spherical_vec(3,:));
q2.AutoScale = 'on'; q2.AutoScaleFactor = 0.3;
hold on;
plot3( x_g(1,1), x_g(2,1), x_g(3,1),'r*')

% view(1,1,0);
