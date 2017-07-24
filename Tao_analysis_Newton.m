%根据油井对水井的单位脉冲响应计算注入水通过不同的tao储层的流量
%给tao设定5个值，利用牛顿方法对这五个值进行优化求解
%分两步优化。第一步先优化n，第二步再优化tao；
clear tao;
clear n;
L=5; %参数tao的个数
%设置参数向量
N=length(data);  %数据的长度
for k=1:L
    tao(k)=k*60;
%     tao(k)=200;
end
% clear n;
for k=1:L
    n(k)=10;
end
loop=1;
step=0.0001;
while(loop<1500)
    %修改参数tao的值
    for k=1:L
        sum=0;        
        for i=1:N
            sum1=0;
            for j=1:L
                sum1=sum1+n(j)*(1-exp(-(i-1)/tao(j)));
            end
            sum=sum+(data(i)-sum1)*n(k)*i/(tao(k)*tao(k))*exp(-(i-1)/tao(k));
        end
        if sum>=0;
            tao(k)=tao(k)-1;
        else
            tao(k)=tao(k)+1;
        end            
        
%         if tao(k)<1
%            tao(k)=3;
%         end
    end
    %修改参数n的值
%     for k=1:L
%         sum=0;        
%         for i=1:N
%             sum1=0;
%             for j=1:L
%                 sum1=sum1+n(j)*(1-exp(-i/tao(j)));
%             end
%             sum=sum+(data(i)-sum1)*n(k)*(-1)*(1-exp(-i/tao(k)));
%         end
%         n(k)=n(k)-sum*step;
%         if n(k)<0
%             n(k)=0;
%         end
%     end
clear M
    for i=1:N
      for j=1:L
         M(i,j)=(1-exp(-(i-1)/tao(j)));
      end      
    end
    Positive_element=eye(L);
    for j=1:L
        if n(j)<0
            Positive_element(j,j)=1000;
        else
            Positive_element(j,j)=0;
        end
    end
    
    n=inv(M'*M+0.21*eye(L))*M'*data;
%     n=inv(M'*M+0.21*eye(L)+Positive_element)*M'*data;
%     for j=1:L
%         if n(j)<0;
%           n(j)=0;
%         end
%     end
    loop=loop+1;
end

for i=1:N
     sum1=0;
     for j=1:L
         sum1=sum1+n(j)*(1-exp(-(i-1)/tao(j)));
     end
    nihe(i)=sum1;
    E(i)=data(i)-sum1;
end

plot(E,'r');
hold on;
plot(data,'b');
hold on;
plot(nihe,'g')
norm(E)
% 
% 
% N_vect=inv(M'*M)*M'*data;
% Nihe=M*N_vect;
% plot(Nihe)
% hold on
% plot(data,'r')