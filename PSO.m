clear M;
clear M1;
clear err_lamd;
clear lamd;
clear tao;
clear sum;
clear vector;
clear tao_vector;
clear injection;
clear times;
%һ�μ��������;���עˮ������ͨ��ϵ��productionΪ�;���ʷ��injectionΪˮ����ʷ��
global data;
global production;
production=data(:,1);
global injection;
injection=data(:,2:end);
global start;
start=10;
global injection_number;
injection_number=length(injection(1,:));
global production_number;
production_number=length(production(1,:));
global Pwf_flag;
Pwf_flag=0;
global ridge;
ridge=0;
global Q;
Q=production(start:end,:);
global lamd;
k=0;
err=0;
particles_num=64;    %Ϊÿ���;��ĸ���������ʼ����������
splitup=2;
%ÿ���;���Ӧ�Ĳ���Ϊinjection_number+1(ѹ���㶨����)��injection_number+2�����Ǿ���ѹ��������
%particles_numΪ��������
%production_numberΪ�;���������
%��ʼ����ʹ�������Ӿ��ȵķֲ��������ռ��ڣ��м������������м���forѭ��������ʼ������λ�ø�ֵ��
if Pwf_flag==0
    %������ѹ���ÿ���;���Ӧ�Ĳ���Ϊ5�������ڴ�ģ����Ϊ5������
    for P_n=1:production_number
    ks=1;
    for i=1:splitup
        for j=1:splitup
            for k=1:splitup
                for m=1:splitup
                    for n=1:splitup
                        for mn=1:splitup
                        particle(1,ks,P_n)=i*10;
                        particle(2,ks,P_n)=j*10;
                        particle(3,ks,P_n)=k*10;
                        particle(4,ks,P_n)=m*10;
                        particle(5,ks,P_n)=n*10;
                        particle(6,ks,P_n)=mn*10;
                        ks=ks+1;        
                        end
                    end
                end
            end
        end
    end
    end

    
else
    for P_n=1:production_num
    ks=1;
    for i=1:splitup
        for j=1:splitup
            for k=1:splitup
                for m=1:splitup
                    for n=1:splitup  
                        for ns=1:splitup
                        particle(1,ks,P_n)=i*10;
                        particle(2,ks,P_n)=j*10;
                        particle(3,ks,P_n)=k*10;
                        particle(4,ks,P_n)=m*10;
                        particle(5,ks,P_n)=n*10;
                        particle(6,ks,P_n)=ns*10;                        
                        ks=ks+1;        
                        end
                    end
                end
            end
        end
    end
    end
end
%��ʼ��ÿ�����ӵĳ�ʼ�ٶȣ�
if Pwf_flag==0
    v=10*rand(injection_number+1,particles_num,production_number);
else
    v=10*rand(injection_number+2,particles_num,production_number);
end
times=0;
w=0.2;
present=particle;%ÿ�����ӵĵ�ǰλ�ã�
pbest=particle;%ÿ���;�ÿ�����Ӿ���������λ�ã�
gbest=particle(:,1,production_number);%ÿ���;����������ӵ�ȫ������λ�ã�
gerr=1000000*ones(production_number);%ÿ���;���ȫ������λ�õ���
perr=1000000*ones(particles_num,production_number);%ÿ���;���ÿ�����ӵ�ȫ������λ�õ���
tao_num=length(present(:,1,1));
while (times<100)      
    for P_n=1:production_number
        for i=1:particles_num  
            for ks=1:tao_num                
               if present(ks,i,P_n)<1
                   present(ks,i,P_n)=1;
               elseif present(ks,i,P_n)>1000
                   present(ks,i,P_n)=1000;
               end
           end
        end
        for k=1:particles_num
            [err(k),lamd(:,k)]=compute_value(P_n,present(:,k,P_n));%k���;�����i�����ӵ����
        end
        for k=1:particles_num
            if err(k)<perr(k,P_n)
               perr(k,P_n)=err(k);
               pbest(:,k,P_n)=present(:,k,P_n);
            end
            if err(k)<gerr(P_n)
                gerr(P_n)=err(k);
                gbest(:,P_n)=present(:,k,P_n);               
            end
            v(:,k,P_n)=w*v(:,k,P_n)+2*rand(1,1)*(pbest(:,k,P_n)-present(:,k,P_n))+3*rand(1,1)*(gbest(:,P_n)-present(:,k,P_n));
            present(:,k,P_n)=present(:,k,P_n)+v(:,k,P_n);
            gerr(1)
        end
    end
    times=times+1;
end
%��ͼ
drawpicture(gbest);

