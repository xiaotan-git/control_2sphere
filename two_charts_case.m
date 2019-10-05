%%%%%%%%%%%%%%  spherical polytopes in one chart.
%%%%%%%%%%%%%%  Both state trajectories and vector field
%%%%%%%%%%%%%% illustrations are included.
clc; clear;
%% initial data
% suppose we just know the feasible spherical polytopes expressed by
% svertices_all and c
% the sample points should be strictly inside the polygons. Needs
% change here. Right now it is working but should be impoved.

% v1 = [-4 -4 -4 -3 -2.5 -2.5 -2 0 1  1 2  2    3   3  4  4    4  1   2;
%       4  1 -4  1 -1.5 -2.5  3 0 2 -1 2 -1 -1.5 -2.5 4 -1.5 -4 -1.5 -1.5];% 2*n
% v1 = 0.3*v1;
% c1 = {[1,7,9,11,15],[1,2,4,7],[2,3,6,5,4],[4,5,18,10,8],[7,8,10,9],[10,18,19,12],...
%      [11,12,19,13],[15,11,13,16],[6,3,17,14],[13,14,17,16]}';
% a1 = [1, -1, 0]'; a1 = a1/norm(a1);
% 
% v2 = [-1  -1     0 0   0 0.75 1.5  0 1.01;
%        1  -0.75  0 0 1.5    1 0.25 0 1.5];% 2*n
% c2 = {[1,3,6,5],[1,2,3],[2,3,4],[3,4,8,7],[7,8,9],[6,7,9],...
%      [5,6,9]}';
% a2 = [-1, -1, 0]'; a2 = a2/norm(a2);
% svertices_all_temp = phi_a_inv(v1,a1); vertices_all_temp = phi_a(svertices_all_temp,a2);
% v2(:,4) = vertices_all_temp(:,3);
% v2(:,8) = vertices_all_temp(:,2);
% c2_temp = cell(size(c2));
% for i = 1: size(c2,1)
%     c2_temp(i,1) = {c2{i,1} + size(v1,2)};
% end
% 
% a = [a1 a2];
% c = [c1'  c2_temp']';
% 
% svertices_a1 = phi_a_inv(v1,a1); 
% svertices_a2 = phi_a_inv(v2,a2);
% svertices_all = [svertices_a1, svertices_a2];
% svertices_chart = ones(size(svertices_all,2),1);
% for i = size(svertices_a1,2)+1: size(svertices_chart,1)
%     svertices_chart(i,1) = 2;
% end
% % vertices_all = phi_a(svertices_all,a2);
% 
% initial_state = [0, 0.58, 0.67 , 0.4, -0.9, -0.56 -0.18;
%                  -1, -0.8, -0.6 -0.34, -0.3, - 0.18 -0.95;
%                  -0.3, -0.15, -0.48,-0.8, 0.3, 0.6 0.25];
% initial_state = normc(initial_state); % normalize the columns
% x_g = [1 0 0.3]'/norm([1 0 0.3]);
% 
% 
% 
% seg_ind ={[1,15],[1,7],[2,4],[10,8],[7,9],[10,18],...
%      [11,13],[15,11],[6,3],[13,16],[20,22],[21,22],[22,23],[23,27],...
%      [26,27],[26,28],[25,28]}';
% time = 4;
% 
% sphe_sample_spoint_temp = []; sphe_sample_spoint = [];
% for i = 1:size(c,1)
%     spherical_v = svertices_all(:,c{i,1});
%     chart = min(svertices_chart(c{i,1},1));
%     if chart == 1
%         v = phi_a(spherical_v,a1);
%         vx = [v(1,:) v(1,1)];
%         vy = [v(2,:) v(2,1)];
%         [px, py] = sampl_poly(vx,vy);
%         sample_point = [px; py]; 
%         sphe_sample_spoint_temp=  phi_a_inv(sample_point,a1); % 3*m
%         sample_spoint_temp = phi_a(sphe_sample_spoint_temp,a1);
%         vertices_a1 = phi_a(svertices_a1,a1);
%         [~,count] = deter_Q(sample_spoint_temp,vertices_a1,c1);
%         sphe_sample_spoint = [ sphe_sample_spoint sphe_sample_spoint_temp(:,count == 1)];
%     else
%         v = phi_a(spherical_v,a2);
%         vx = [v(1,:) v(1,1)];
%         vy = [v(2,:) v(2,1)];
%         [px, py] = sampl_poly(vx,vy);
%         sample_point = [px; py]; 
%         sphe_sample_spoint_temp=  phi_a_inv(sample_point,a2); % 3*m  
%         sample_spoint_temp = phi_a(sphe_sample_spoint_temp,a2);
%         vertices_a2 = phi_a(svertices_a2,a2);
%         [~,count] = deter_Q(sample_spoint_temp,vertices_a2,c2);
%         sphe_sample_spoint = [ sphe_sample_spoint sphe_sample_spoint_temp(:,count == 1)];
%     end
% end


% sample_spoint_temp = phi_a(sphe_sample_spoint_temp,a);
% vertices_all = phi_a(svertices_all,a);
% [~,count] = deter_Q(sample_spoint_temp,vertices_all,c);
% sphe_sample_spoint = sphe_sample_spoint_temp(:,count == 1);

% svertices_all = round(svertices_all,7);
load('initial_data_2charts.mat');
% save('initial_data_2charts.mat', 'a','svertices_chart','initial_state','seg_ind','sphe_sample_spoint',...
%     'svertices_all','c','x_g','time');

%%  state trajectories
% only one point x in R(3,1)is allowed for now.
x_0 = initial_state(:,1);
x_1 = initial_state(:,2);
x_2 = initial_state(:,3);
x_3 = initial_state(:,4);
x_4 = initial_state(:,5);
x_5 = initial_state(:,6);
x_6 = initial_state(:,7);
x_7 = [ 0 -0.95 0.3]'; x_7 = x_7/norm(x_7); 
u = fb_cons_two_charts(x,a,svertices_all,svertices_chart,c,seg_ind, x_g); 
% where we need to use the accociate chart.

% simulation, something wrong with x_4 and x_5
[x_t,t] = sphe_sim(x_5,x_g,a,svertices_all,c,seg_ind,10,svertices_chart);

% plotting
figure
sphere
axis equal
for i = 1:length(c) 
svertices = svertices_all(:,c{i});
draw_spolygon(svertices);
end
axis off
plot3( -x_g(1,1), -x_g(2,1), -x_g(3,1),'go');
plot3( x_1(1,1), x_1(2,1), x_1(3,1),'go');
plot3( x_2(1,1), x_2(2,1), x_2(3,1),'go');
plot3( x_3(1,1), x_3(2,1), x_3(3,1),'go');
plot3( x_4(1,1), x_4(2,1), x_4(3,1),'go');
plot3( x_5(1,1), x_5(2,1), x_5(3,1),'go');
plot3( x_6(1,1), x_6(2,1), x_6(3,1),'go');
plot3( x_7(1,1), x_7(2,1), x_7(3,1),'go');
plot3( x_g(1,1), x_g(2,1), x_g(3,1),'r*');
plot3(x_t(1,:), x_t(2,:), x_t(3,:),'b-');


%% vector field illustration
% calculate dot_x for all sphe_sample_spoint
spherical_vec = zeros(3, size(sphe_sample_spoint,2));
for i = 1 : size(sphe_sample_spoint,2)
    spherical_one_point =  sphe_sample_spoint(:,i);
    [ ~ , dot_sp] = fb_cons_two_charts(spherical_one_point,a,svertices_all,svertices_chart,...
                                    c,seg_ind, x_g);
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
q2.AutoScale = 'on'; q2.AutoScaleFactor = 0.5; q2.LineWidth = 0.2;
hold on;
plot3( x_g(1,1), x_g(2,1), x_g(3,1),'r*')

% view(1,1,0);
