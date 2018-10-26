%[xx,yy,w_f,theta_f,t_vec,l2_f_vec,l2_w_vec]=vorticity_solve();

[xx,yy,l2_f_vec,l2_w_vec]=vorticity_solve();

X(1,:)=xx;
X(2,:)=yy;

save('coords.mat','X')
save('omega.mat','w_f')
save('theta.mat','theta_f')

T(1,:)=t_vec;
T(2,:)=l2_f_vec;
T(3,:)=l2_w_vec;

save('norms.mat','T')
