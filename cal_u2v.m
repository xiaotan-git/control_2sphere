function v = cal_u2v(u,sp,opt)
%cal_u2v, given u in R(6,m), sp in R(3,m), opt in R(3,1)
%   Detailed explanation goes here
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

m = size(u,2);
v = zeros(2,m);
%  sp = spoint; i = 100;
for i = 1:m
    grad_phi = J2*mat*(a'*sp(:,i)*I - sp(:,i)*a')/((a'*sp(:,i))^2); 
    Pi_p =  cross_vec(sp(:,i));
    Sig = grad_phi*Pi_p;
    v(:,i) = (Sig *u(:,i))/norm(Sig *u(:,i));
%     u(:,i) = Sig'*inv(Sig*Sig')*v(:,i);
%     dot_sp(:,i) = Pi_p*u(:,i);
end
end