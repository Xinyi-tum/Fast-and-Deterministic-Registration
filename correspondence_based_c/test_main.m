clc;clear;close all


%%%%check compile
if exist('./FBnB_registration','file') == 0
    command1 = 'bash ./compile.sh';
    [status1,cmdout1] = system(command1);
end

num_all = 100;
outlier_rate = 0.5;
num_inlier = round(num_all*(1-outlier_rate));
num_outlier = round(num_all-num_inlier);
noise_level = 0.1;
M = 1;
N = 0;
epsilon = 0.0175*1;
t_domain = [0;0;0;200];
sett = [epsilon, t_domain'];
fid = fopen('parameter_syn.txt','w');
fprintf(fid,'%.10f,%.10f,%.10f,%.10f,%.10f\n',sett);
for i = 1:M
    [data_x,data_y,R_gt,t_gt]= gen_data_pcd(num_inlier,num_outlier,noise_level);
    command1 = './FBnB_registration ./input_syn.txt ./parameter_syn.txt ./results_syn.txt';
    [status1,cmdout1] = system(command1);
    results = load ('results_syn.txt');
    A = results(1:9);
    B = reshape(A,[3,3]);
    R_opt = B.';
    t_opt = results(10:12)';
    e_r(i) = acosd((trace(transpose(R_gt)*R_opt)-1)/2);
    e_t(i) = norm(R_opt*t_opt-R_gt*t_gt);
    time(i) = results(13);
    
end
E_r = mean(e_r);
E_t = mean(e_t);
Tim = median(time);

disp(['Points amount: ',num2str(num_inlier+num_outlier)]);
disp(['Outlier rate: ',num2str(num_outlier/(num_inlier+num_outlier)*100),'(%)']);
disp(['Noise level: ',num2str(noise_level),newline]);
disp(['Rot. error of ours: ',num2str(E_r),'(deg.)']);
disp(['Trans. error of ours: ', num2str(E_t), '(m)']);
disp(['Running time of ours: ',num2str(Tim),'(ms)',newline]);

