function draw_spolygon(svertices)
%DRAW_SPOLYGON given svertices in R(3,n);
%   Detailed explanation goes here

if ~((1- min(vecnorm(svertices))<10E-6)||(max(vecnorm(svertices))-1<10E-6))
        error('draw_spolygon function: the spherical vertices are not on the sphere')
end
hold on;
n = size(svertices,2);
for i = 1:(n-1)
    A = svertices(:,i);
    B = svertices(:,i+1);

    theta = acos(A'*B);
    v = cross(A,B);
    sp_int = zeros(3,21);
    for j = 1:21
        r = [v; theta*(j-1)/20];
        mat = vrrotvec2mat(r);
        sp_int(:,j) = mat*A;
    end
    plot3(sp_int(1,:),sp_int(2,:),sp_int(3,:),'k','LineWidth',2);
end
A = svertices(:,n);
B = svertices(:,1);

theta = acos(A'*B);
v = cross(A,B);
sp_int = zeros(3,21);
for j = 1:21
    r = [v; theta*(j-1)/20];
    mat = vrrotvec2mat(r);
    sp_int(:,j) = mat*A;
end
plot3(sp_int(1,:),sp_int(2,:),sp_int(3,:),'k','LineWidth',2);
end

