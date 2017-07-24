%井组为单位进行连通性分析，daa文件第一列为中心井油井的产油量，四周第二、三、四、五列分别为四口水井。

% 四口水井的 tao(扩散系数)
tao1=0
tao2=0
tao3=0
tao4=0
% 四口水井的权值系数
lam1=0
lam2=0
lam3=0
lam4=0

lamda=0
err=1000000000;
dt=50
maxdata1=0
maxdata2=0
product=ones(length(data(:,1)),1);
injection1=ones(length(data(:,1)),1);
injection2=ones(length(data(:,1)),1);
injection3=ones(length(data(:,1)),1);
injection4=ones(length(data(:,1)),1);
L=ones(length(data(:,1))-2,1);
L=data(3:end,1);
product=data(:,1);
injection1=data(:,2));
injection2=data(:,3);
injection3=data(:,4);
injection4=data(:,5);

for tao1=1:30
    for tao2=1:30
        for tao3=1:30
            for tao4=1:30
                
               for i=3:length(data(:,1))           
                  sum1=0;
                  sum2=0;
                  sum3=0;
                  sum4=0;
                  for j=1:i
                      sum1=sum1+exp((j-i)/tao1)*injection1(j)/tao1;
                      sum2=sum2+exp((j-i)/tao2)*injection2(j)/tao2;
                      sum3=sum3+exp((j-i)/tao3)*injection3(j)/tao3;
                      sum4=sum4+exp((j-i)/tao4)*injection4(j)/tao4;
                  end                    
                      pre1(i-2)=product(1)*exp((1-i)/tao1)+sum1;
                      pre2(i-2)=product(1)*exp((1-i)/tao2)+sum2;
                      pre3(i-2)=product(1)*exp((1-i)/tao3)+sum3;
                      pre4(i-2)=product(1)*exp((1-i)/tao4)+sum4;                       
               end
               M(:,1)=pre1;
               M(:,2)=pre2;
               M(:,3)=pre3;
               M(:,4)=pre4;              
               lamd=inv(M'*M)*M'*L;
               err_L=norm(M*lamd-L);
               if err_L<err
                   lamd_o=lamd;
                   tao=(tao1,tao2,tao3,tao4);
                   err=err_L;
               end
               
               
            end
        end
    end
end


