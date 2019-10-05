function point = phi_a(spoint,opt)
%PHI_A, given spoint in R(3,m), opt in R(3,1)
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

if (max(vecnorm(spoint)) -1 > 10E-6)||(1- min(vecnorm(spoint))>10E-6)
    error('phi_a function: check the spherical points');
end

m = size(spoint,2); point = zeros(2,m);     
theta = acos(e3'*a);
v = cross(a,e3);
r = [v; theta];
mat = vrrotvec2mat(r);

for i = 1:m
   if a'*spoint(:,i) <= 10E-6
       point(:,i) = [0 0]';
%        error('phi_a function: check the spherical points');
   else
       point(:,i) = J2*mat*spoint(:,i)/(a'*spoint(:,i));       
   end
end
end
