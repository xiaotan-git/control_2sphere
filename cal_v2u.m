function [u, dot_sp] = cal_v2u(v,sp,opt)
% cal_v2u, given v in R(2,m), sp in R(3,m), opt in R(3,1)
%  no normaliziation in this file
J2 = [ 1 0 0; 0 1 0]; e3 = [0,0,1]'; I = diag([1,1,1]);
if nargin >3
    error('dis_p2l function: too many inputs');
end
switch nargin
    case 2
        a = e3;
    case 3
        a = opt;
end

theta = acos(e3'*a);
v_r = cross(a,e3);
r = [v_r; theta];
mat = vrrotvec2mat(r);

m = size(v,2);
u = zeros(3,m); dot_sp = zeros(3,m);
for i = 1:m
    grad_phi = J2*mat*(a'*sp(:,i)*I - sp(:,i)*a')/((a'*sp(:,i))^2); 
    Pi_p = cross_vec(sp(:,i)); % dynamics inbedded, maybe a problem
    Sig = grad_phi*Pi_p;
    u(:,i) = Sig'*inv(Sig*Sig')*v(:,i); % badly singular, should be changed
    dot_sp(:,i) = Pi_p*u(:,i);
end
end