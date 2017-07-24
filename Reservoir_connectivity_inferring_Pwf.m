%井组为单位进行连通性分析，daa文件第一列为中心井油井的产油量，四周第二、三、四、五列分别为四口水井。
clear M;
clear M1;
clear lamd;
clear tao;
clear sum;
clear vector;
clear tao_vector;
clear injection;
clear pre;
clear pre1;
clear pre_5;

N=2; %注水井井数
colum=[2,3];%各注水井的列数
step=1;%步长
% tao[1),tao[2),tao3...分别是相应的水井的扩散系数，taoc是油井的扩散系数
% 四口水井的权值系数
err=1000000000;
for i=1:N
    tao(i)=10;    
end
tao_c=10;
product=ones(length(data(:,1)),1);
for i=1:N
  injection(:,i)=ones(length(data(:,1)),1);
end
L=ones(length(data(:,1))-2,1);
L=data(3:end,1);
product=data(:,1);
for i=1:N
  injection(:,i)=data(:,colum(i));
end


   for i=3:length(data(:,1))           

                  for k=1:N
                      sum(k)=0;
                  end
                  for k=1:N
                      sum(k)=sum(k)+exp((1-i)/tao(k))*injection(1,k)/tao(k)/2;
                  end
                  for j=2:i-1
                      for k=1:N
                          sum(k)=sum(k)+exp((j-i)/tao(k))*injection(j,k)/tao(k);
                      end
                  end 
                  
                  for k=1:N
                      sum(k)=sum(k)+injection(i,k)/tao(k)/2;
                  end

                  for j=1:N                  
                     pre(i-2,j)=sum(j);                  
                  end
                  pre_5(i-2)=product(1)*exp((1-i)/tao_c);
               end
               M=ones(length(pre(:,1)),N);

%                M(:,5)=pre';
               for j=1:N                  
                  M(:,j)=pre(:,j);                  
               end
%                M(:,N+1)=ones(length(pre(:,1)),1);
               
              k=M(:,1)'*M(:,1);
              unit=k*eye(length(M(1,:)));
               lamd=inv(M'*M+unit)*M'*(L-pre_5');
               err_L=norm(M*lamd-L+pre_5');       
% while(norm(d_vector1+d_vector2)>0.00001)
err1=0;
times=1;
while(times<100)    
    %            计算各个tao的变化
    for k=1:N
           tao11=tao(k)+step;
           M1=M;
            for i=3:length(data(:,1))           
               sum1=0;                  
                  sum1=sum1+exp((1-i)/tao11)*injection(1,k)/tao11/2;
                  for j=2:i-1
%                       sum1=sum1+exp((j-i)/tao11)*injection(j,k)/tao11;
                        sum1=sum1+exp((j-i)/tao11)*injection(j,k)/tao11;
                  end
                  sum1=sum1+injection(i,k)/tao11/2;                
                  pre1(i-2)=sum1;                     
             end
               M1(:,k)=pre1';                       
               lamd=inv(M1'*M1+unit)*M1'*(L-pre_5');
               err=norm(M1*lamd-L+pre_5');
               if err_L>err
                   dtao=step;
               else
                   if tao(k)>1
                      dtao=-1*step;
                    else
                       dtao=0;
                   end
               end        
               tao(k)=tao(k)+dtao;
       end       
       
  %修改生产数据的参数     
       tao11=tao_c+step;
           M1=M;
            for i=3:length(data(:,1))           
               pre_5(i-2)=product(1)*exp((1-i)/tao11);                                       
             end                                      
               lamd=inv(M1'*M1+unit)*M1'*(L-pre_5');
               err=norm(M1*lamd-L+pre_5');
               if err_L>err
                   dtao=step;
               else
                   if tao_c>1
                      dtao=-1*step;
                    else
                       dtao=0;
                   end
               end        
               tao_c=tao_c+dtao;     

 %计算修改参数后的误差 
            for i=3:length(data(:,1))
                  for j=1:N
                      sum(j)=0;
                  end
                  for k=1:N
                      sum(k)=sum(k)+exp((1-i)/tao(k))*injection(1,k)/tao(k)/2;
                  end
                  for j=2:i-1
                      for k=1:N
                          sum(k)=sum(k)+exp((j-i)/tao(k))*injection(j,k)/tao(k);
                      end
                  end                   
                  for k=1:N
                      sum(k)=sum(k)+injection(i,k)/tao(k)/2;
                  end
                  for k=1:N                  
                     pre(i-2,k)=sum(k);                  
                  end
                  pre_5(i-2)=product(1)*exp((1-i)/tao_c);
           end
              
               M=pre;                                    
%                M(:,N+1)=ones(length(pre(:,1)),1);
               k=1;
               unit=k*eye(length(M(1,:)));
               lamd=inv(M'*M+unit)*M'*(L-pre_5');
               err1=err_L;
               err_L=norm(M*lamd-L+pre_5');    
               
                if err_L>err1
                   step=step/4;
               end 
    times=times+1;
end
tao_vector=tao';
vector(:,1)=tao_vector;
vector(:,2)=lamd;
plot(M*lamd-L+pre_5','red');
hold on;
plot(L,'green');
hold on;
plot(M*lamd+pre_5','blue');

%直接用最小二乘计算
c1=data(:,2:length(data(1,:)));
lamd1=inv(c1'*c1)*c1'*data(:,1);
figure(2);
plot(c1*lamd1-data(:,1),'red');
hold on;
plot(data(:,1),'green');
hold on;
plot(c1*lamd1,'blue');