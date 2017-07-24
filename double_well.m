maxdt=0
lamda=0
least=100000;
dt=50
maxdata1=0
maxdata2=0
data1=ones(length(data(:,1)),1);
data1=data(:,2);
data2=ones(length(data(:,1)),1);
data2(1,1)=data(1,2);
    for tao=1:30
        err=0;
        for i=2:length(data(:,1))           
           sum=0;
           for j=2:i-1
               sum=sum+exp((j-i)/tao)*data(j,1);
           end
           data2(i)=0.5*data(1,1)*exp((1-i)/tao)+0.5*data(i,1)+sum; 
           data2(i)=data2(i)/tao+data(1,1)*exp((1-i)/tao);
        end
        err=norm(data1-data2);
       if err<least 
           least=err;
           maxdt=dt
           lamda=tao
           maxdata1=data1
           maxdata2=data2
       end
    end    
    erro=maxdata2-data(:,2);