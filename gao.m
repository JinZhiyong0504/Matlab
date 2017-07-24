
N=182;% 节点数
time=8192;%时间步数
data=zeros(time,N);
data(1,1)=1;%初始边界条件；
%计算每个节点的参数值
for i=1:N-1
  if mod(i-1,20)==0
      rn(i)=46.4163;
  else
      rn(i)=8.822698;
%       rn(i)=8;
  end
end
%初始化
for i=2:N-1
    data(2,i)=data(1,i+1)*2*rn(i)/(rn(i)+rn(i-1))+data(1,i-1)*2*rn(i-1)/(rn(i)+rn(i-1));    
end
    data(2,N)=data(1,N-1);
    data(2,1)=0;
%计算
for j=2:time-1
    for i=2:N-1  
        data(j+1,i)=data(j,i+1)*2*rn(i)/(rn(i)+rn(i-1))+data(j,i-1)*2*rn(i-1)/(rn(i)+rn(i-1))-data(j-1,i);            
    end
    data(j+1,1)=0;
    data(j+1,N)=data(j,N-1);
end
%求傅立叶变换,对第s个节点的震动进行傅立叶变换,s变化范围[1,89]
% s=181;
% Y = fft(data(:,s),time);
% Pyy = Y.* conj(Y);
% % f = 1000*(0:time/2)/time;
% % plot(f,Pyy(1:time/2+1));
% 
% f = 1:time;
% plot(f,Pyy);

% figure(2);
% plot(f(800:1300),Pyy(800:1300));


s=182
Fs=1/0.000089;
x=data(:,s);
X=fftshift(fft(x));
f=(1:time)*Fs/(time)-Fs/2;
plot(f,abs(X)); 
title('Frequency content of y')
xlabel('frequency (Hz)')
figure(2)
t=1:length(data(:,s));
t=t*0.000089;
plot(t,data(:,s));



% t=0:0.001:2;
% Fs=1000;
% Fc=250;
% x=cos(2*pi*Fc*t);
% X=fftshift(fft(x));
% f=(0:2000)*Fs/2001-Fs/2;
% plot(f,abs(X)) 