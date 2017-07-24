% add fractal water cut to matching.Modified CRM-EnOpt mathematic % this version can calculate connectivity of wells
% and parameters in oil production equations in mulitple wells.
% Ye=[tao,tao11,tao12,tao13,tao14,lamd1,lamd2,lamd3,lamd4,f,f1,f2,f3,f4,QL,QO],
% multi production well parameters inferring with constant BHP.
% injection; fluid prodution and oil production was arranged separately
% Data include all the fluid production volume and oil production volume and
% injection volume;data=[production,injection];
% function [output]=CRM_EnOpt_new()
% global injection;
% global fluid_pro;
% global oil_pro;
% global Num_Inj;    % number of injection well
% global Num_Pro;    % number of production well;
% global Ne;         % number of realization;
% global Num;        % total number of parameters;
% global Ye;         % realizations of parameters;

timebegin=50;               % 
timewidow=50;               % time window for computation;
alltime=size(injection,1);  % all time 
t_front=20;                 % start of calculation
t_end=500;                  % end of calculation
Num_Inj=size(injection,2);                 % initialize injection numbers
Num_Pro=size(fluid_pro,2);                 % initialize production numbers
Ne=50; %number of realization;
Num=1+Num_Inj+Num_Inj+1+Num_Inj+2;  % 1:production well tao; 3: a,b,wi; 2:fluid production and oil production; 
Param_Index=ones(Num,Num_Pro); %record the optimal resolution for each producer;
min_errors=ones(Num_Pro)*10000000000;
err_new=100000000000000000000000;


output_step=ones(2,Ne,Num_Pro);
e=ones(2,Ne,Num_Pro);
% store computation result;
Com_result=zeros(t_end,2,Ne,Num_Pro);

V_ne=ones(Ne,Ne)/Ne;

H=zeros(1,Num);
H(1,Num-1)=1;

Ho=zeros(1,Num);
Ho(1,Num)=1;

P1=zeros(Num,Num);
P1(1:2*Num_Inj+1,1:2*Num_Inj+1)=eye(2*Num_Inj+1,2*Num_Inj+1);
P2=zeros(Num,Num);
P2(2*Num_Inj+2:3*Num_Inj+2,2*Num_Inj+2:3*Num_Inj+2)=eye(Num_Inj+1,Num_Inj+1);

 err=1000000000000; %overall error;

%initializing state vector.% state vector include Ye=[tao,tao11,tao12,tao13,tao14,lamd1,lamd2,lamd3,lamd4,f,f1,f2,f3,f4,QL,QO]
for i=1:Num_Pro
    Ye(1:Num_Inj+1,1:Ne,i)=25*rand(Num_Inj+1,Ne)+1;                      %tao
    Ye(Num_Inj+1+1:Num_Inj+1+Num_Inj,1:Ne,i)=rand(Num_Inj,Ne);           %lamd
    Ye(Num_Inj+1+Num_Inj+1:3*Num_Inj+2,1:Ne,i)=ones(Num_Inj+1,Ne);       % f,f11,f12,f13
    Ye(Num_Inj+1+Num_Inj+Num_Inj+1+1:Num,1:Ne,i)=zeros(2,Ne);       % fluid production and oil production;
end

% iterative with time
% for i=1:Num_Pro
    for t=t_front+1:t_end  
       
        % calculation output
        for i=1:Num_Pro                       
            mid=Comp_Mlutiproduction(t,i);     
            Com_result(t,:,:,i)=mid;              
            Ye(Num-1:Num,:,i)=mid+random('Normal',0,1,2,Ne); % assign the Ye with output and noise.           
            for j=1:Ne      
                e(1,j,i)=fluid_pro(t,i)-Ye(Num-1,j,i);   
                e(2,j,i)=oil_pro(t,i)-Ye(Num,j,i);
                
                %calculate the error and record the optimal resolution    
                err_new=norm(e(1,:,1))*norm(e(2,:,1));
                if(err_new<min_errors(i))
                    output(:,i)=Ye(:,j,i);
                    min_errors(i)=err_new;
                end
            end                 
           
            %calculate covariance matrix;
%           Cd1=e(1,:,i)*e(1,:,i)'/(Ne-1);
            CY=(Ye(:,:,i)-Ye(:,:,i)*V_ne)*(Ye(:,:,i)-Ye(:,:,i)*V_ne)'/(Ne-1);                
            Cd1=0;
            Cd2=0; 
            % update Ye. pay attention to the noise which is impulsory.
            for j=1:Ne        
                Ye(:,j,i)=Ye(:,j,i)+P1*CY*H'/(H*CY*H'+Cd1)*e(1,j,i)+P1*0.01*random('Normal',0,1,Num,1);                 
                if oil_pro(t,i)<fluid_pro(t,i)-0.1
                    Ye(:,j,i)=Ye(:,j,i)+P2*CY*Ho'/(Ho*CY*Ho'+Cd2)*e(2,j,i)+P2*0.01*random('Normal',0,1,Num,1); 
                end    
            end        
            
        end
        
%         if(err_new<err)
%          output=Ye(:,:,i);
%         end
        err=err_new
    end
    
   for k=1:Num_Pro 
       figure(k);
    for i=1:Ne
        plot(Com_result(:,1,i,k));
        hold on;
    end    
    plot(fluid_pro(1:500,k),'r');
    
   hold on;
    for i=1:Ne
        plot(Com_result(:,2,i,k));
        hold on;
    end    
     plot(oil_pro(1:500,k),'r');
   end   
    
% end