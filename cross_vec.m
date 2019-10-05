function M = cross_vec(v)
%CROSS_VEC asymetric matrix induced by v
%   此处显示详细说明
v1 = v(1,1); v2 = v(2,1);  v3 = v(3,1); 
M(1,1) = 0;  M(1,2) = -v3;  M(1,3) = v2; 
M(2,1) = v3;  M(2,2) = 0;  M(2,3) = -v1; 
M(3,1) = -v2;  M(3,2) = v1;  M(3,3) = 0; 

end

