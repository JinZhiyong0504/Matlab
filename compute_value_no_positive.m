function [err,lamd]=compute_value_no_positive(tao_v,data)
tao=tao_v(1);
tao1=tao_v(2);
tao2=tao_v(3);
tao3=tao_v(4);
tao4=tao_v(5);
penalty=100000;
k=0.2;
product=ones(length(data(:,1)),1);
injection1=ones(length(data(:,1)),1);
injection2=ones(length(data(:,1)),1);
injection3=ones(length(data(:,1)),1);
L=ones(length(data(:,1))-2,1);
L=data(3:end,1);
product=data(:,1);
injection1=data(:,2);
injection2=data(:,3);
injection3=data(:,4);


            signal=1;
            for j=1:length(tao_v)              
               if tao_v(j)<=0.01
                 signal=0;
               end
           end
 if (signal==0)
     lamd=1000*ones(length(tao_v),1);
     err=100000000;
 else
           for i=3:length(data(:,1))           
                  sum1=0;
                  sum2=0;
                  sum3=0;

                  sum1=sum1+exp((1-i)/tao1)*injection1(1)/tao1/2;
                  sum2=sum2+exp((1-i)/tao2)*injection2(1)/tao2/2;
                  sum3=sum3+exp((1-i)/tao3)*injection3(1)/tao3/2;

                  for j=2:i-1
                      sum1=sum1+exp((j-i)/tao1)*injection1(j)/tao1;
                      sum2=sum2+exp((j-i)/tao2)*injection2(j)/tao2;
                      sum3=sum3+exp((j-i)/tao3)*injection3(j)/tao3;

                  end                    
                  sum1=sum1+injection1(i)/tao1/2;
                  sum2=sum2+injection2(i)/tao2/2;
                  sum3=sum3+injection3(i)/tao3/2;

                 
                  pre5(i-2)=product(1)*exp((1-i)/tao);
                  pre1(i-2)=sum1;
                  pre2(i-2)=sum2;
                  pre3(i-2)=sum3;
                     
               end
               M=ones(length(pre1),3);

%                M(:,5)=pre';
               M(:,1)=pre1';
               M(:,2)=pre2';
               M(:,3)=pre3';

%                M(:,5)=ones(length(pre1),1);
               
                
%                oen=eye(length(M(1,:)));
%                unit=k*oen*(M'*M);
%                M'*M
%                tao_v
                Mt=M'*M;

               lamd=inv(Mt+k*Mt(1,1)*eye(length(Mt(:,1))))*M'*(L-pre5');
           
%                for k=1:length(lamd)
%                                    
%                    if (lamd(k)<1)&&(lamd(k)>0)
%                        err_lamd(k)=0;
%                    elseif lamd(k)<0
%                        err_lamd(k)=-1*lamd(k);
%                    else                      
%                        err_lamd(k)=lamd(k)-1;
%                    end
%                end                      
%                
%                for k=1:length(tao_v)
%                    if tao_v(k)<0
%                        err_tao(k)=-1*tao_v(k)+1;
%                    else
%                        err_tao(k)=0;                       
%                    end
%                end
%                err_L=norm(M*lamd-L+pre5')+penalty*norm(err_lamd)+penalty*norm(err_tao);          
               err_L=norm(M*lamd-L+pre5');      
               err=err_L;
               
     end

