%�����;���ˮ���ĵ�λ������Ӧ����ע��ˮͨ����ͬ��tao���������
%��tao�趨5��ֵ������ţ�ٷ����������ֵ�����Ż����
%�������Ż�����һ�����Ż�n���ڶ������Ż�tao��
clear tao;
clear n;
L=5; %����tao�ĸ���
%���ò�������
N=length(data);  %���ݵĳ���
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
    %�޸Ĳ���tao��ֵ
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
    %�޸Ĳ���n��ֵ
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