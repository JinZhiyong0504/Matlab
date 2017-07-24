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
clear err_lamd;
N=1;
step=0.5;
penalty=0;
% 四口水井的 tao(扩散系数)
ini=5
tao=ini;
tao1=ini;
tao2=ini;
tao3=ini;
tao4=ini;
tao5=ini;
% 四口水井的权值系数
err=1000000000;
product=ones(length(data(:,1)),1);
injection1=ones(length(data(:,1)),1);
injection2=ones(length(data(:,1)),1);
injection3=ones(length(data(:,1)),1);
injection4=ones(length(data(:,1)),1);
L=ones(length(data(:,1))-2,1);
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
               
   Mt=M'*M;
   lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M'*(L-pre5');
              for k=1:length(lamd)-1
                   if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k)+1;
                   else                      
                       err_lamd(k)=lamd(k);
                   end
               end               
               err_L=norm(M*lamd-L+pre5')+100000*norm(err_lamd);

% while(norm(d_vector1+d_vector2)>0.00001)
err1=0;

while(abs(err1-err_L)>0.1)
    
    %              计算tao1的变化方向 
           tao11=tao1+step;
           if tao11<0
               tao11=0.01;
           end
           M1=M;
            for i=3:length(data(:,1))           
               sum1=0;                  
                  sum1=sum1+exp((1-i)/tao11)*injection1(1)/tao11/2;
                  for j=2:i-1
                      sum1=sum1+exp((j-i)/tao11)*injection1(j)/tao11;
                  end
                  sum1=sum1+injection1(i)/tao11/2;                
                  pre1(i-2)=sum1;                     
               end
               M1(:,1)=pre1';                       
               Mt=M1'*M1;
               lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M1'*(L-pre5');
   
               for k=1:length(lamd)
                   if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k);
                   else               
                       err_lamd(k)=lamd(k)-1;
                   end
               end               
               err=norm(M1*lamd-L+pre5')+penalty*norm(err_lamd);
               if err_L>err
                   dtao1=step;
               else
                   if tao1>1
                      dtao1=-1*step;
                    else
                       dtao1=0;
                   end
               end               
               % 计算tao2的变化方向 
                tao21=tao2+step;
                if tao21<0
                    tao21=0.01;
                end
                M1=M;
               for i=3:length(data(:,1))           
                 sum2=0;
                 sum2=sum2+exp((1-i)/tao21)*injection2(1)/tao21/2;                  
                 for j=2:i-1
                      sum2=sum2+exp((j-i)/tao21)*injection2(j)/tao21;
                  end       
                  sum2=sum2+injection2(i)/tao21/2;
                  pre2(i-2)=sum2;                    
               end
               M1(:,2)=pre2';             
               Mt=M1'*M1;
               lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M1'*(L-pre5');
               for k=1:length(lamd)
                   if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k);
                   else
                       err_lamd(k)=lamd(k)-1;
                   end
               end               
               err=norm(M1*lamd-L+pre5')+penalty*norm(err_lamd);
               if err_L>err
                   dtao2=step;
               else
                   if tao2>1
                       dtao2=-1*step;
                   else
                       dtao2=0;
                   end
               end
               
               %计算tao3的变化方向 
               tao31=tao3+step;
               if tao31<0
                    tao31=0.01;
                end
               M1=M;
               for i=3:length(data(:,1))           
               sum3=0;
               sum3=sum3+exp((1-i)/tao31)*injection3(1)/tao31/2;                  
                 for j=2:i-1
                      sum3=sum3+exp((j-i)/tao31)*injection3(j)/tao31;
                  end   
               sum3=sum3+injection3(i)/tao31/2;
               pre3(i-2)=sum3;                 
               end

               M1(:,3)=pre3';            
               Mt=M1'*M1;
               lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M1'*(L-pre5');
               for k=1:length(lamd)
                   if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k);
                   else
                       err_lamd(k)=lamd(k)-1;
                   end
               end               
               err=norm(M1*lamd-L+pre5')+penalty*norm(err_lamd);
               if err_L>err
                   dtao3=step;
               else
                   if tao3>1
                       dtao3=-1*step;
                   else
                       dtao3=0;
                   end
               end          
               
%              计算tao4的变化方向  
               tao41=tao4+step;
               if tao41<0
                    tao41=0.01;
                end
               M1=M;
               for i=3:length(data(:,1))           
               sum4=0;
               sum4=sum4+exp((1-i)/tao41)*injection4(1)/tao41/2;
                  for j=2:i-1
                      sum4=sum4+exp((j-i)/tao41)*injection4(j)/tao41;
                  end   
               sum4=sum4+injection4(i)/tao41/2;
               pre4(i-2)=sum4;                       
               end

               M1(:,4)=pre4';              
               Mt=M1'*M1;
               lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M1'*(L-pre5');
               for k=1:length(lamd)
                   if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k);
                   else
                       err_lamd(k)=lamd(k)-1;
                   end
               end               
               err=norm(M1*lamd-L+pre5')+penalty*norm(err_lamd);
               if err_L>err
                   dtao4=step;
               else
                   if tao4>1
                       dtao4=-1*step;
                   else
                       dtao4=0;
                   end
               end
%                              计算tao5的变化方向  
               tao51=tao5+step;
               if tao51<0
                    tao51=0.01;
                end
               M1=M;
               for i=3:length(data(:,1))           
               sum5=0;
               sum5=sum5+exp((1-i)/tao41)*pwf(1)/tao51/2;
                  for j=2:i-1
                      sum5=sum5+exp((j-i)/tao51)*pwf(j)/tao51;
                  end   
               sum5=sum5+pwf(i)/tao5/2+pwf(1)*exp((1-i)/tao5)-pwf(1);
               pwf_c(i-2)=sum5;                       
               end

               M1(:,5)=pwf_c';              
               Mt=M1'*M1;
               lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M1'*(L-pre5');
               for k=1:length(lamd)
                   if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k);
                   else
                       err_lamd(k)=lamd(k)-1;
                   end
               end               
               err=norm(M1*lamd-L+pre5')+penalty*norm(err_lamd);
               if err_L>err
                   dtao4=step;
               else
                   if tao4>1
                       dtao4=-1*step;
                   else
                       dtao4=0;
                   end
               end
 %              计算tao的变化方向  
               tao_1=tao+step;
               if tao_1<0
                    tao_1=0.01;
                end
               M1=M;
               for i=3:length(data(:,1))           
                 pre_c(i-2)=product(1)*exp((1-i)/tao_1);                       
               end            
               Mt=M1'*M1;
               lamd=inv(Mt+ridge*Mt(1,1)*eye(length(Mt(:,1))))*M1'*(L-pre5');
               for k=1:length(lamd)
                  if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k);
                   else
                       err_lamd(k)=lamd(k)-1;
                   end
               end               
               err=norm(M1*lamd-L+pre5')+penalty*norm(err_lamd);
               if err_L>err
                   dtao=step;
               else
                   if tao>1
                       dtao=-1*step;
                   else
                       dtao=0;
                   end
               end
%                修正之后计算矩阵和向量
               d_vector2=[dtao1;dtao2;dtao3;dtao4];
               tao1=tao1+dtao1;
               tao2=tao2+dtao2;
               tao3=tao3+dtao3;
               tao4=tao4+dtao4;
               tao=tao+dtao;
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
               for k=1:length(lamd)
                   if (lamd(k)<1)&&(lamd(k)>0)
                       err_lamd(k)=0;
                   elseif lamd(k)<0
                       err_lamd(k)=-1*lamd(k);
                   else
                       err_lamd(k)=lamd(k)-1;
                   end
               end                           
               err1=err_L;
               err_L=norm(M*lamd-L+pre5')+penalty*norm(err_lamd);
               if err_L>err1
                   step=step/4;
               end            
               
    N=N+1;
end
figure(1);
tao_vector=[tao1,tao2,tao3,tao4,tao5];
vector(:,1)=tao_vector';

vector(:,2)=lamd(1:length(tao_vector));
figure(1);
plot(M*lamd-L+pre5','red');
hold on;
plot(L,'green');
hold on;
plot(M*lamd+pre5','blue');


%直接用最小二乘计算
figure(2);
c1=data(:,2:length(data(1,:)));
lamd1=inv(c1'*c1)*c1'*data(:,1);
figure(2);
plot(c1*lamd1-data(:,1),'red');
hold on;
plot(data(:,1),'green');
hold on;
plot(c1*lamd1,'blue');
err_lamd1=norm(c1*lamd1-data(:,1));
vector1=lamd1;
% for tao1=1:30
%     for tao2=1:30
%         for tao3=1:30
%             for tao4=1:30
%                 
%                for i=3:length(data(:,1))           
%                   sum1=0;
%                   sum2=0;
%                   sum3=0;
%                   sum4=0;
%                   for j=1:i
%                       sum1=sum1+exp((j-i)/tao1)*injection1(j)/tao1;
%                       sum2=sum2+exp((j-i)/tao2)*injection2(j)/tao2;
%                       sum3=sum3+exp((j-i)/tao3)*injection3(j)/tao3;
%                       sum4=sum4+exp((j-i)/tao4)*injection4(j)/tao4;
%                   end                    
%                       pre1(i-2)=product(1)*exp((1-i)/tao1)+sum1;
%                       pre2(i-2)=product(1)*exp((1-i)/tao2)+sum2;
%                       pre3(i-2)=product(1)*exp((1-i)/tao3)+sum3;
%                       pre4(i-2)=product(1)*exp((1-i)/tao4)+sum4;                       
%                end
%                M=ones(length(pre1),4);
%                M(:,1)=pre1';
%                M(:,2)=pre2';
%                M(:,3)=pre3';
%                M(:,4)=pre4';              
%                lamd=inv(M'*M)*M'*L;
%                err_L=norm(M*lamd-L);
%                if err_L<err
%                    lamd_o=lamd;
%                    tao=[tao1 tao2 tao3 tao4];
%                    err=err_L;
%                end
%                
%                
%             end
%         end
%     end
% end


