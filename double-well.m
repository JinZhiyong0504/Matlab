t=0
lamda=0
least=100000
% for dt=1:30 
% %     for tao=0.1:0.1:30
% %         err=0;
% %         for i=1:length(data(:,1))-dt
% %            data1(i)=data(i+dt,2)
% %            sum=0
% %            for j=i:i+dt
% %                sum=data(j,3)*exp((j-i+dt)/tao)
% %            end
% %            data2(i)=data(i,2)*exp(-dt/tao)+sum/tao;           
% %         end
% %         err=norm(data1-data2)
% %        if err<least 
% %            least=err;
% %            t=dt;
% %            lamda=tao;
% %        end
% %     end
% dt=dt+1
% end
