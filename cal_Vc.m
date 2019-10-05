function V_c = cal_Vc(point,f_ex,f_ey,vertices,opt)
%calculate V_c, given m p in R(2,m),f_ex in R(2,1),f_ey in R(2,1),  
% n vertices in R(2,n);
if nargin >5
    error('cal_V function: too many inputs');
end


f_v1 = [f_ex(1,1); f_ey(1,1) ];
f_v2 = [f_ex(2,1); f_ey(2,1) ];
[exis_flag1,exis_ind1] = ismembertol(f_v1',vertices','Byrows',true);
[exis_flag2,exis_ind2] = ismembertol(f_v2',vertices','Byrows',true);
vect_size = size(vertices,2);

if ~(exis_flag1&&exis_flag2)
    error('cal_Vc function: did not find the vertice of the exit edge');
end

if ~(( abs(exis_ind1 - exis_ind2)== 1)||((exis_ind2 == 1)&&(exis_ind1 == vect_size)||((exis_ind1 == 1)&&(exis_ind2 == vect_size))))
        error('cal_Vc function: the vertice of the exit edge do not match!')
end



switch nargin
    case 4
        %  the four points in the order CABD 
        if (min(exis_ind1,exis_ind2) == 1)&&(max(exis_ind1,exis_ind2) == vect_size)
            B = vertices(:,vect_size);
            A = vertices(:,1);    
            C = vertices(:,2);
            D = vertices(:,vect_size - 1);
        else
            exis_ind_min = min(exis_ind1,exis_ind2);
            exis_ind_max = max(exis_ind1,exis_ind2);
            B = vertices(:,exis_ind_min);
            A = vertices(:,exis_ind_max);
            if exis_ind_min == 1
                C = vertices(:,3);
                D = vertices(:,vect_size);
            else
                if exis_ind_max == vect_size
                    C = vertices(:,1);
                    D = vertices(:,vect_size - 2);
                else
                    C = vertices(:,exis_ind_max + 1);
                    D = vertices(:,exis_ind_min - 1);  
                end
            end
        end

        if (A - C)'*(B-D)/(norm(A - C)*norm(B-D))> 1 - 10E-6
            p = (A + B)/2 + 0.5* (A-C);
        else
            Mat_temp = [C(2,1) - A(2,1), -(C(1,1) - A(1,1));
                D(2,1) - B(2,1), -(D(1,1) - B(1,1))];
            Vec_temp = [((C(2,1) - A(2,1))*A(1,1) -(C(1,1) - A(1,1))*A(2,1));
                ((D(2,1) - B(2,1))*B(1,1) -(D(1,1) - B(1,1))*B(2,1))];
            p_intersec = inv(Mat_temp)*Vec_temp;
            if ((B(2,1) - A(2,1))*(p_intersec(1,1) - A(1,1)) - (B(1,1) - A(1,1))*(p_intersec(2,1) - A(2,1)))*((B(2,1) - A(2,1))*(C(1,1) - A(1,1)) - (B(1,1) - A(1,1))*(C(2,1) - A(2,1))) >0
                p = p_intersec + 2*((A + B )/2 -  p_intersec);
            else
                p = p_intersec + 1/2*((A + B )/2 -  p_intersec);
            end
        end
        V_c = (p - point)./vecnorm(p - point);
    case 5
        xi_g = opt;
        V_c_temp = (xi_g - point);
        V_c = ((V_c_temp)./vecnorm(V_c_temp));
%         V_c = zeros(size(V_c));
end
end

