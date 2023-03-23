function [R_opt,R_vot,theta_opt,amounts] = voting(data_x,data_y,corr_opt,t_opt,v_p,v_q)
a = acos(dot(v_p,v_q));
axis = cross(v_p,v_q);
axis = axis./norm(axis);
R_pq = rotationVectorToMatrix(-a*axis);%rotation with the minimum geodesic motion % v_q=R_pq*v_p;

M = sum(corr_opt>0);
data_y(:,corr_opt==0) = [];
data_x = data_x(:,corr_opt(corr_opt>0));
v1 = R_pq*(data_x+t_opt);

m_ = cross(v1,repmat(v_q,1,M));
m = m_./vecnorm(m_);
n_ = cross(data_y,repmat(v_q,1,M));
n = n_./vecnorm(n_);

theta = acos(dot(m,n));
v2 = data_y;
flag = v_q'*cross(v1,v2)<0;%check direction
theta(flag) = -theta(flag);

edges = linspace(-pi,pi,361);

% figure
% histogram(theta,edges);
% title('voting for rotation')
% xticks([-pi -pi/2 0 pi/2 pi])
% xticklabels({'-\pi' '-\pi/2' '0' '\pi/2' '\pi'})
% xlabel('rad');
% grid on

angle_hist = histcounts(theta,edges);%histogram voting[-pi:0.0175:pi]
[amounts,theta_index] = max(angle_hist);
theta_opt = 0.5*(edges(theta_index)+edges(theta_index+1));

R_vot = rotationVectorToMatrix(-theta_opt*v_q);
R_opt = R_vot*R_pq;
end

