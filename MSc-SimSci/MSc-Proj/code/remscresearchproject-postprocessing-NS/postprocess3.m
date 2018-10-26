function [xa,ya]=postprocess3(nlow,nhigh,nres)

num_str=num2str(nlow);
filename_theta=strcat('output_theta',num_str);
filename_w=strcat('output_vorticity',num_str);

w=open_octave_loc(filename_w);
t=open_octave_loc(filename_theta);

[d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad]=postprocess0(w,t);

%  [yl,xl]=hist(reshape(Xangle,512^2,1),-(pi/2):0.01:(9*pi/2));
%  [x_ss,y_ss,~,~]=gaussian_smoothing_loc(xl,yl,0.15);

[yw,xw]=hist(reshape(phi,nres^2,1),0:0.01:(pi/2));
x_ss=xw;
y_ss=yw;

sum_vec_x=x_ss;
sum_vec_y=y_ss;



for k=(nlow+1):1:(nhigh)
    k
    filectr=k;
    num_str=num2str(filectr);
    filename_theta=strcat('output_theta',num_str);
    filename_w=strcat('output_vorticity',num_str);
    
    w=open_octave_loc(filename_w);
    t=open_octave_loc(filename_theta);
    
    [d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad]=postprocess0(w,t);
    
%     [yl,xl]=hist(reshape(Xangle,512^2,1),-(pi/2):0.01:(9*pi/2));    
%     [x_ss,y_ss,~,~]=gaussian_smoothing_loc(xl,yl,0.15);

    [yw,xw]=hist(reshape(phi,nres^2,1),0:0.01:(pi/2));
    x_ss=xw;
    y_ss=yw; 
    
    sum_vec_x=sum_vec_x+x_ss;
    sum_vec_y=sum_vec_y+y_ss;
    
end

xa=sum_vec_x/(nhigh-nlow);
ya=sum_vec_y/(nhigh-nlow);








end

function C=open_octave_loc(inputfilename)

Nh=5;

fid=fopen(inputfilename);
C_temp=textscan(fid, '%f', 'HeaderLines',Nh);
fclose(fid);

C=C_temp{1};

C=reshape(C,512,512);

end

function [x_ss,y_ss,x_a,y_a]=gaussian_smoothing_loc(x,y,w)

% Gaussian filter for a function on a periodic domain.

x=x-2*pi;
L=x(length(x))-x(1);

n=length(x);
x_a=0*(1:(3*n));
y_a=x_a;

dx=abs(x(2)-x(1));

for i=1:n
    x_a(i)=x(i)-L;
    y_a(i)=y(i);
end

for i=1:n
    x_a(i+n)=x(i);
    y_a(i+n)=y(i);
end

for i=1:n
    x_a(i+2*n)=x(i)+L;
    y_a(i+2*n)=y(i);
end

y_s=y_a;

for j=1:(3*n)
    y_s(j)=sum(G(x_a-x_a(j)).*y_a)/sum(G(x_a-x_a(j))); 
end

x_ss=0*(1:n);
y_ss=x_ss;
for i=1:n
    x_ss(i)=x_a(i+n);
    y_ss(i)=y_s(i+n);
end

for i=1:floor(n/2)
    x_ss(i+floor(n/2))=-x_ss(i);
    y_ss(i+floor(n/2))=y_ss(i);
end

[x_ss1,ix]=sort(x_ss);

for i=1:length(x_ss1)
    y_ss1(i)=y_ss(ix(i));
end

x_ss=x_ss1;
y_ss=y_ss1;

nrm=sum(y_ss)*dx;
y_ss=y_ss/nrm;

    function z=G(x)
        z=exp(-x.^2/(2*w*w));
        %z=exp(-abs(x)/w);
    end

end
