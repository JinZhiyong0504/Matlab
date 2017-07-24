clear M;
clear M1;
clear err_lamd;
clear lamd;
clear tao;
clear sum;
clear vector;
clear tao_vector;
clear injection;
clear pre;
clear pre1;
clear pre_5;

err=0;
num=30  %Á£×ÓÊı

v=rand(6,num);
times=0;
w=0.5;
present=50*rand(6,num);
pbest=present;
gbest=present(:,1);
lamd=0;
% v=w*v+3*rand(0,1)*(pbest-present)+2*rand(0,1)*(gbest-present);
% present=present+v;
gerr=1000000;
perr=1000000*ones(1,num);
while (times<2000)
     for i=1:num
        for j=1:length(present(:,i))
           if present(j,i)<=0.01
               present(j,i)=1;
           end
       end
        [err,lamd]=compute_pwf_value(present(:,i),data);
        if err<perr(i)
            perr(i)=err;
            pbest(:,i)=present(:,i);
        end
        if err<gerr
            gerr=err;
            gbest=present(:,i);
            lamd_best=lamd;
        end
       v(:,i)=w*v(:,i)+10*rand(1,1)*(pbest(:,i)-present(:,i))+10*rand(1,1)*(gbest-present(:,i));
       present(:,i)=present(:,i)+v(:,i);
              
    end   
    times=times+1;
%     gbest
end

tao=gbest(1);
tao1=gbest(2);
tao2=gbest(3);
tao3=gbest(4);
tao4=gbest(5);
tao5=gbest(6);
L=data(3:end,1);
product=data(:,1);
injection1=data(:,2);
injection2=data(:,3);
injection3=data(:,4);
injection4=data(:,5);
pwf=data(:,6);
ridge=0.2;
  for i=3:length(data(:,1))           
                  sum1=0;
                  sum2=0;
                  sum3=0;
                  sum4=0;
                  sum5=0;
                  sum1=sum1+exp((1-i)/tao1)*injection1(1)/tao1/2;
                  sum2=sum2+exp((1-i)/tao2)*injection2(1)/tao2/2;
                  sum3=sum3+exp((1-i)/tao3)*injection3(1)/tao3/2;
                  sum4=sum4+exp((1-i)/tao4)*injection4(1)/tao4/2;
                  sum5=sum5+exp((1-i)/tao5)*pwf(1)/tao5/2;
                  for j=2:i-1
                      sum1=sum1+exp((j-i)/tao1)*injection1(j)/tao1;
                      sum2=sum2+exp((j-i)/tao2)*injection2(j)/tao2;
                      sum3=sum3+exp((j-i)/tao3)*injection3(j)/tao3;
                      sum4=sum4+exp((j-i)/tao4)*injection4(j)/tao4;
                      sum5=sum5+exp((j-i)/tao5)*pwf(j)/tao5;
                  end                    
                  sum1=sum1+injection1(i)/tao1/2;
                  sum2=sum2+injection2(i)/tao2/2;
                  sum3=sum3+injection3(i)/tao3/2;
                  sum4=sum4+injection4(i)/tao4/2;
                  sum5=sum5+pwf(i)/tao5/2+pwf(1)*exp((1-i)/tao5)-pwf(1);
%                   sum1=sum1+injection1(i);
%                   sum2=sum2+injection2(i);
%                   sum3=sum3+injection3(i);
%                   sum4=sum4+injection4(i);

                  
                  pre5(i-2)=product(1)*exp((1-i)/tao);
                  pre1(i-2)=sum1;
                  pre2(i-2)=sum2;
                  pre3(i-2)=sum3;
                  pre4(i-2)=sum4;    
                  pwf_c(i-2)=sum5;
               end
               M=ones(length(pre1),4);

%                M(:,5)=pre5';
               M(:,1)=pre1';
               M(:,2)=pre2';
               M(:,3)=pre3';
               M(:,4)=pre4';
               M(:,5)=pwf_c';
%                M(:,5)=pre5';
%                M(:,5)=ones(length(pre1),1);
               Mt=M'*M;
               lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M'*(L-pre5');
               
               
figure(1);

plot(M*lamd-L+pre5','red');
hold on;
plot(L,'green');
hold on;
plot(M*lamd+pre5','blue');

c1=data(:,2:length(data(1,:)));
lamd1=inv(c1'*c1)*(c1'*data(:,1));
figure(2);
plot(c1*lamd1-data(:,1),'red');
hold on;
plot(data(:,1),'green');
hold on;
plot(c1*lamd1,'blue');


vector(:,1)=gbest(2:end);
vector(:,2)=lamd_best;
vector1=lamd1
% err=compute_value(tao(1),tao(1),tao(2),tao(2),tao(2),data);