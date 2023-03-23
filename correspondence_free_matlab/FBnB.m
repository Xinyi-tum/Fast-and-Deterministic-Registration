function [t_opt,L_global,corr] = FBnB(data_x,data_y,epsilon,t_lb,t_ub,v_p,v_q)

t_B = [t_lb;t_ub];
best_branch = t_B;
B = [];
M = size(data_x,2);
N = size(data_y,2);
L_global = 0;
U_global = max(M,N);
iter = 1;

while(U_global-L_global>0)
    new_B_wo_bound = Divide(best_branch);
    new_B = estimate_LU(data_x,data_y,epsilon,new_B_wo_bound,v_p,v_q);
    B = [B,new_B];
    [L_now,ind_L_now] = max(B(end-1,:));
    [U_global,ind_U] = max(B(end,:));
    if(L_now > L_global)
        L_global = L_now;
        opt_branch = B(:,ind_L_now);
    end
    best_branch = B(:,ind_U);
    B(:,ind_U) = [];
    B(:,B(end,:) < L_global) = [];
end
t_opt = 0.5*(opt_branch(1:3)+opt_branch(4:6));

if(U_global==L_global) 
        data_x_tran = data_x+t_opt;
        Length = vecnorm(data_x_tran);
        ang_alpha = acos(v_p'*(data_x_tran./Length));
        ang_gama = acos(v_q'*(data_y)./(vecnorm(data_y)));
        
        matching = abs(ang_alpha'-ang_gama)<=epsilon;
        corr = dmperm(matching>0);
end
end

