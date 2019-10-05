function d = dis_p2l(p, lx,ly,opt1)
%% return the distance of a point to a line segment;
% when flag = 1, we need just one point  p in R(2,1)
% p in R(2,1) or R(2,m); lx in R(2,m); ly in R(2,m); return d in R(1,m)
% p = rand(2,1); lx = rand(2,5);ly = rand(2,5);
if nargin >4
    error('dis_p2l function: too many inputs');
end
switch nargin
    case 3
        flag = 0;
    case 4
        flag = opt1;
end

m = size(lx,2);
if flag == 1
    p = p*ones(1,m);
end

v_BA = [lx(1,:)-lx(2,:);ly(1,:) - ly(2,:)];
v_AB = -v_BA;
v_BP = [p(1,:)-lx(2,:);p(2,:) - ly(2,:)];
v_AP = [p(1,:)-lx(1,:);p(2,:) - ly(1,:)];
theta_ABP = acos(diag(v_BA'*v_BP)./(vecnorm(v_BA).*vecnorm(v_BP))');


d_orth = vecnorm(v_BP)'.*sin(theta_ABP);
d_AP = vecnorm(v_AP)';
d_BP = vecnorm(v_BP)';

for i = 1:m
    if (v_BA(:,i)'*v_BP(:,i))*(v_AB(:,i)'*v_AP(:,i)) >0
        d(i) = d_orth(i);
    else
        d(i) = min(d_AP(i),d_BP(i));
    end
end

% d
% d_orth'
% min(d_AP,d_BP)
% 
% plot(lx,ly); legend;
% hold on;
% plot(p(1,5),p(2,5),'r+')
end


