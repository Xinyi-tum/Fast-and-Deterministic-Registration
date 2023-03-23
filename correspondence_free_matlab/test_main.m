clc;clear;close all

num = 100;
overlap = 0.8;
noise_level = 0.01;
M = 1;%number of trials
N = 0;
epsilon = 0.0175*0.5;%threshold of the inlier
for i = 1:M
    
    [data_x,data_y,R_theta,R_v,R_gt,t_gt,v_p,v_q,corr_gt] = gen_data_overlap(num,overlap,noise_level);%generate synthetic data
    t_lb = [-100;-100;-100];%translation domain
    t_ub = [100;100;100];
    tic
    [t_opt,L_global,corr_opt] = FBnB(data_x,data_y,epsilon,t_lb,t_ub,v_p,v_q);%translation search
    time1(i) = toc;
    tic
    [R_opt,R_vot,theta_opt,amounts] = voting(data_x,data_y,corr_opt,t_opt,v_p,v_q);%rotation voting
    time2(i) = toc;
    
    eR = transpose(R_gt)*R_opt;
    e_r(i) = acosd((trace(eR)-1)/2);%rotation error
    e_t(i) = norm(t_opt-t_gt);%translation error
end
E_r = mean(e_r);
E_t = mean(e_t);
raito = N/M;%success rate
Time1 = median(time1);
Time2 = median(time2);

disp(['Points amount: ',num2str(num)]);
disp(['Noise level: ',num2str(noise_level)]);
disp(['Rot. error: ',num2str(E_r),'(deg.)']);
disp(['Trans. error: ', num2str(E_t), '(m)']);
disp(['Median time for trans.:',num2str(Time1),'(s)']);
disp(['Median time for Rot.:',num2str(Time2),'(s)',newline]);
