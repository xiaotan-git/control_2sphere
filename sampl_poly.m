function [pinx,piny] = sampl_poly(lx,ly)
% sampling the polygon with uniform samping points
xmin = min(min(lx)); xmax = max(max(lx));
ymin = min(min(ly));ymax = max(max(ly));

pointx = xmin + 0.001: 0.05 :xmax - 0.001;
pointy = ymin + 0.001:0.05:ymax- 0.001;
[a,b] = size(pointx);
[c,d] = size(pointy);
m = a*b*c*d;
[pX, pY] = meshgrid(pointx,pointy);
px = reshape(pX,1,m)';
py = reshape(pY,1,m)';
[in,on] = inpolygon(px,py,lx,ly);
% plot(px(in),py(in),'+')

pinx = px(in)'; piny = py(in)';

end

