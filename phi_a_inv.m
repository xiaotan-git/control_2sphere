function spoint = phi_a_inv(point,opt)
%PHI_A, given spoint in R(3,m), opt in R(3,m)
%   Detailed explanation goes here
J2 = [ 1 0 0; 0 1 0]; e3 = [0,0,1]';
if nargin >3
    error('dis_p2l function: too many inputs');
end
switch nargin
    case 1
        a = e3;
    case 2
        a = opt;
end


m = size(point,2); spoint = zeros(3,m);     
theta = acos(e3'*a);
v = cross(a,e3);
r = [v; theta];
mat = vrrotvec2mat(r);

for i = 1:m
    spoint(:,i) = mat'*(J2'*point(:,i)+e3)/(sqrt(1+point(:,i)'*point(:,i)));
end
end
