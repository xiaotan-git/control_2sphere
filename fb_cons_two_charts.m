function [u,dot_sp] = fb_cons_two_charts(x,a,svertices_all,svertices_chart,c,seg_ind, x_g)
% FB_CONS_TWO_CHARTS for one state x on the sphere,
% calculate the control input  u(x)

% determin which projection should be used and whether the 
% state belongs to the overlapping region of two charts.

% suppose that the discrete plan is from chart 2 to chart 1, known in
% advance

 a1 = a(:,1);
 a2 = a(:,2);
 n = sum(svertices_chart( svertices_chart == 1)); % no.svertices in chart 1
 m = 0;
 for i = 1:size(c,1)
     m = m + min(c{i,1} <= n); % no.spherical polytopes in chart 1
 end
svertices_a1 = svertices_all(:,1:n);
c1 = c(1:m,1);
vertices_a1 = phi_a(svertices_a1,a1);
    
svertices_a2 = svertices_all(:,n+1:end);

seg_ind1 = cell(m,1);
for i = 1: m
    c1(i,1) = {c{i,1}};
    b = seg_ind{i,1};
    seg_ind1(i,1) = {b} ;
end

seg_ind2 = cell(size(c,1)-m,1);
for i = 1: size(c,1)-m
    c2(i,1) = {c{i+m,1} - n};
    b = seg_ind{i+m,1}-n;
    seg_ind2(i,1) = {b};
end
    
vertices_a2 = phi_a(svertices_a2,a2);

[svertices_share, ia, ib] = intersect(svertices_a1',svertices_a2','rows');
svertices_share = svertices_share';

 for i = 1:m
     if min(ismember(ia,c1{i,1}))
        Q_ind1_share = i;
     end
 end
 for i = 1:size(c,1)-m
     if min(ismember(ib,c2{i,1}))
        Q_ind2_share = i;
     end
 end

u = [ 0 0 0]'; dot_sp = [ 0 0 0]';

if (x'*a1 <= 0)&&(x'*a2 <= 0) % not in both charts
    error('fb_cons_two_charts error: check the input state, which lies out side the two charts.');
else
    if  x'*a1 > 0
        xi1 = phi_a(x,a1);
        Q_ind1 = deter_Q(xi1,vertices_a1,c1,1); 
        if (Q_ind1 ~= 0) % in chart 1, and in Q_ind1
            [u,dot_sp] = fb_cons_one_chart(x,a1,svertices_a1,c1,seg_ind1, x_g);
        end        
    end
    if x'*a2 >0
            xi2 = phi_a(x,a2);
            Q_ind2 = deter_Q(xi2,vertices_a2,c2,1); % seems something wrong with my points.
            if (Q_ind2 ~= 0) % in chart 2, and in Q_ind2
                [u,dot_sp] = fb_cons_one_chart(x,a2,svertices_a2,c2,seg_ind2, x_g);
            end
%         else     %  in chart 1 and chart 2, but not in any polytope
%             error('fb_cons_two_charts error: check the input state, which lies out side the two charts.');
    end 
    if (x'*a1 > 0)&&(x'*a2 >0)&&((Q_ind1 == 0 && Q_ind2 == Q_ind2_share )||(Q_ind1 == Q_ind1_share && Q_ind2 == 0 ))
        if  Q_ind1 == Q_ind1_share
            vertices_a = vertices_a1(:,c1{Q_ind1_share});
            vertices_b = vertices_a2(:,c2{Q_ind2_share});
            
            vx = vertices_a(1,:); vx(1,size(vertices_a,2)+1) = vx(1,1); vx = vx';
            vy = vertices_a(2,:); vy(1,size(vertices_a,2)+1) = vy(1,1); vy = vy';
            shift_vx = circshift(vx,1);
            lx = [vx(2:end,1) shift_vx(2:end,1)]';
            shift_vy = circshift(vy,1);
            ly = [vy(2:end,1) shift_vy(2:end,1)]';
            d = dis_p2l(xi1,lx,ly);
            [~,ind_dmin] = min(d);
            
            v_share = phi_a(svertices_share,a1);
            L = [lx; ly];
            L_share = [v_share(1,1) v_share(1,2); v_share(1,2) v_share(1,1); ...
                       v_share(2,1) v_share(2,2); v_share(2,2) v_share(2,1)];
            L = round(L,4); L_share = round(L_share,4); % this is a trick and may not work in some circumstances
            [~, ia, ~] = intersect(L',L_share','rows');        
            if ind_dmin == ia
                 % from cell 2 TO cell 1, considering the state x in cell 1
                f_ex_share = v_share(1,:)'; f_ey_share = v_share(2,:)';
                f_ex = vertices_a1(1,seg_ind1{Q_ind1_share})';
                f_ey = vertices_a1(2,seg_ind1{Q_ind1_share})';                
                [u, dot_sp] = cal_trans_chart( x,a1,a2,f_ex_share,f_ey_share, f_ex, f_ey,vertices_a, vertices_b, svertices_share); 
                disp('now it is in region 1,a; chart transition is needed.');
            end
%         else
%             vertices = vertices_a2(:,c2{Q_ind2});
%             vx = vertices(1,:); vx(1,size(vertices,2)+1) = vx(1,1); vx = vx';
%             vy = vertices(2,:); vy(1,size(vertices,2)+1) = vy(1,1); vy = vy';
%             shift_vx = circshift(vx,1);
%             lx = [vx(2:end,1) shift_vx(2:end,1)]';
%             shift_vy = circshift(vy,1);
%             ly = [vy(2:end,1) shift_vy(2:end,1)]';
%             d = dis_p2l(xi2,lx,ly);
%             [~,ind_dmin] = min(d);
%             
%             v_share = phi_a(svertices_share,a2);
%             L = [lx; ly];
%             L_share = [v_share(1,1) v_share(1,2); v_share(1,2) v_share(1,1); ...
%                        v_share(2,1) v_share(2,2); v_share(2,2) v_share(2,1)];
%             L = round(L,4); L_share = round(L_share,4); % this is a trick and may not work in some circumstances
%             [~, ia, ~] = intersect(L',L_share','rows');        
%             if ind_dmin == ia
%                 flag = 2; % from cell 2 TO cell 1
%                 f_ex = v_share(1,:); f_ey = v_share(2,:);                
%                 [u, dot_sp] = cal_trans_chart(x,a1,a2,f_ex,f_ey,vertices,Q_ind1_share,Q_ind2_share,...
%                                     svertices_a1, svertices_a2, x_g, flag); 
%                 disp('now it is in region 2,b');
%             end
        end
    end
 end

end
