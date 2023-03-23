function [data_x,data_y,R_gt,t_gt] = gen_data_pcd(num_inlier,num_outlier,noise_level)

num = num_inlier+num_outlier;
R_theta = (rand(1)*2-1)*pi;
R_gt = [cos(R_theta),sin(R_theta),0;-sin(R_theta),cos(R_theta),0;0,0,1];
% v = rand(3,1)*2-1;
% R_v = v./norm(v);
% R_gt = rotationVectorToMatrix(R_theta*R_v);
scale = 100;
t_gt = (rand(3,1)*2-1)*scale;
% T_gt = [R_gt,t_gt]';
% fid = fopen('gt_syn.txt','w');
% fprintf(fid,'%.10f %.10f %.10f %.10f\n',T_gt);
data_x_inlier = (rand(3,num_inlier)*2-1)*scale;
data_y_inlier = R_gt*(data_x_inlier+t_gt);
v = data_x_inlier(:,1)-data_x_inlier(:,2);
v_p = v./norm(v);
v_q = R_gt*v_p;
data_x_outlier = (rand(3,num_outlier)*2-1)*scale;
data_y_outlier = (rand(3,num_outlier)*2-1)*scale;
data_x = [data_x_inlier,data_x_outlier]+normrnd(0,noise_level,3,num_inlier+num_outlier);
data_y = [data_y_inlier,data_y_outlier]+normrnd(0,noise_level,3,num_inlier+num_outlier);%
%%save data
v=[v_p;v_q]';
fid = fopen('input_syn.txt','w');
fprintf(fid,'%d\n',num);
fid = fopen('input_syn.txt','a');
fprintf(fid,'%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n',v);
for ii = 1:num
    data_point = [data_y(:,ii);data_x(:,ii)]';
    fprintf(fid,'%.10f %.10f %.10f %.10f %.10f %.10f\n',data_point);
end
fclose('all');

data_x2 = data_x';
data_y2 = data_y';
ptCloud1 = pointCloud(single(data_x2(:,1:3)));
pcwrite(ptCloud1, 'point_p.pcd', 'Encoding', 'ascii');  
ptCloud2 = pointCloud(single(data_y2(:,1:3)));
pcwrite(ptCloud2, 'point_q.pcd', 'Encoding', 'ascii'); 
 
end





