function [new_B] = estimate_LU(data_x,data_y,epsilon,new_B_wo_bound,v_p,v_q)
N = size(new_B_wo_bound,2);
L(1,N) = 0;
U(1,N) = 0;
for ii = 1:N
    t_l = new_B_wo_bound(1:3,ii);
    t_u = new_B_wo_bound(4:6,ii);  
    t_c = 0.5*(t_l+t_u);
    
    L(ii) = get_lower_bound(data_x,data_y,epsilon,t_c,v_p,v_q);
    U(ii) = get_upper_bound(data_x,data_y,epsilon,t_c,t_l,t_u,v_p,v_q);
    
end
new_B = [new_B_wo_bound;L;U];
end

function UU = get_upper_bound(data_x,data_y,epsilon,t_c,t_l,t_u,v_p,v_q)
N = size(data_x,2);
ang_alpha = acos(v_p'*(data_x+t_c)./(vecnorm(data_x+t_c)));
ang_gama = acos(v_q'*(data_y)./(vecnorm(data_y)));
delta = sqrt(3)*abs(t_u-t_l)/2;
delta = delta(1);
dis = vecnorm(data_x+t_c);
r = ones(1,N)*pi;
flag2 = delta<dis;
r(flag2) = asin(delta./dis(flag2));
UU = sprank((abs(ang_alpha'-ang_gama)-r')<=epsilon);

end

function [LL] = get_lower_bound(data_x,data_y,epsilon,t_c,v_p,v_q)
ang_alpha = acos(v_p'*(data_x+t_c)./(vecnorm(data_x+t_c)));
ang_gama = acos(v_q'*(data_y)./(vecnorm(data_y)));
LL = sprank(abs(ang_alpha'-ang_gama)<=epsilon);
end
