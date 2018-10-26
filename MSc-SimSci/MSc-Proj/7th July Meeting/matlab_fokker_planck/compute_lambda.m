function [lambda,p_lambda]=compute_lambda(x,y,z,p_joint)

Lx=abs(2*min(x));
Ly=abs(2*min(y));
Lz=abs(2*min(z));
dz=abs(z(2)-z(1));
dy=abs(y(2)-y(1));
Ny=length(y);
dx=abs(x(2)-x(1));

p_xy=zeros(length(x),length(y));

for i=1:length(x)
    for j=1:length(y)
        p_xy(i,j)=sum(p_joint(i,j,:))*dz;
    end
end

lambda_range=-100:0.01:100;
p_lambda=0*lambda_range;

for i=1:length(lambda_range)
    sum_val=0;
    lambda_val=lambda_range(i);
    
    if(lambda_val==0)
        sum_val=0;
    else
    
        for ii=1:length(x)
            g_val=-2*sin(x(ii));
            if((abs(lambda_val)/(abs(g_val)+1e-8))>=(Ly/2))
                summand=0;
            else
                ll_map=1+(1/dy)*((lambda_val/g_val)+(Ly/2));
                ll_map=floor(ll_map);
                
                if(abs(ll_map)>Ny)
                    summand=0;
                else
                    summand=p_xy(ii,ll_map)*(dx/abs(g_val));
                end
            end
            sum_val=sum_val+summand;
            
        end

    end
    
    p_lambda(i)=sum_val;
end

[~,ix]=min(abs(lambda_range));

p_lambda(ix)=(p_lambda(ix-1)+p_lambda(ix+1))/2;

lambda=lambda_range;

end