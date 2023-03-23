function [data_x,data_y,R_theta,R_v,R_gt,t_gt,v_p,v_q,corr_gt] = gen_data_overlap(num,overlap,noise_level)

num_overlap = num*overlap;
R_theta=  (rand(1)*2-1)*3;
v = rand(3,1)*2-1;
R_v = v./norm(v);
R_gt = rotationVectorToMatrix(R_theta*R_v);
scale = 100;
t_gt = (rand(3,1)*2-1)*scale;

data_x_inlier = (rand(3,num)*2-1)*scale;
data_y_inlier_ = R_gt*(data_x_inlier+t_gt);

corr_gt_inlier = randperm(num);
data_y_inlier = data_y_inlier_(:,corr_gt_inlier);

v = data_x_inlier(:,1)-data_x_inlier(:,2);%random gravity direction
v_p = v./norm(v);
v_q = R_gt*v_p;

data_x = data_x_inlier+normrnd(0,noise_level,3,num);
data_y = data_y_inlier(:,1:num_overlap)+normrnd(0,noise_level,3,num_overlap);

corr_gt = corr_gt_inlier;
end
