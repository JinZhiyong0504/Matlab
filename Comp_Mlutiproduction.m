% comput ncth oil well's fluid production and oil productiong with parameters of Ye which contain Ne realization at t step;
% Ye=[tao,tao11,tao12,tao13,tao14,lamd1,lamd2,lamd3,lamd4,f,f1,f2,f3,f4,QL,QO],
function [output]=Comp_Mlutiproduction(t,nc)
global Num_Inj;    % number of injection well;
global injection;  % injection vector;
global fluid_pro;  % fluid production vector;
global Ne; % number of realization;
global Ye;  % realizations of parameters;

fluid_pro_nc=fluid_pro(:,nc);

Cal_range=10;
% 
% 
%     for i=1:Num_Inj
%         u=injrate(t,i);
%         B={[1/P(i+1)]};
%         F={[1,1/P(i+1)]};
%         m0 = idpoly(1,B,1,1,F,1,0);
%         md = c2d(m0,1);
% %         M(:,i)=sim(md,u,'InitialState',prodrate(1));
%         M(1,i)=sim(md,u,'InitialState',0);
%     end   
%     
%      Depletion_item(t)=prodrate(t-1)*exp(-1)/Ye(1));
%      
%      lamd=Ye(num_inj+2:num_inj+1+num_inj)';
%      output(1,1)=M*lamd+Depletion_item'                          % output fluid production;  
%      
%      %compute oil production
%      wl=Ye(end);
%      wl=injrate(t,:)*lamd+wi;
%      Q_o=Production(t)/(1+a*power(wi,b));
%      output(1,2)=Q_o;   


% tao which is nagative is unphsical;
for p=1:Ne
    
    for j=1:1+Num_Inj
        if Ye(j,p,nc)<1;
            Ye(j,p,nc)=1;
        end
    end
    
    for j=2+Num_Inj:2*Num_Inj+1
        if Ye(j,p,nc)<0
            Ye(j,p,nc)=0;
        end
        if Ye(j,p,nc)>1
            Ye(j,p,nc)=1;                
        end
    end  
%     
%      for j=length(Ye(:,1))-4:length(Ye(:,1))-3
%         if Ye(j,p)<0;
%             Ye(j,p)=0;
%         end
%     end
end
     
   a=0;
   b=0;
   P=0;
   lamd=0;
   M_c=0;
   M=ones(Cal_range,Num_Inj);
   output=0;
   wl=0;
   Q_o=0;
     
for n=1:Ne
%    a=Ye(2*Num_Inj+2,n,nc);
%    b=Ye(2*Num_Inj+3,n,nc);
   P=Ye(2:1+Num_Inj,n,nc);
   lamd=Ye(Num_Inj+2:2*Num_Inj+1,n,nc);
   f_wcut=Ye(2*Num_Inj+2:3*Num_Inj+2,n,nc);   
   
    for k=1:Num_Inj
        u=injection(t-Cal_range+1:t,k);
        B={[1/P(k)]};
        F={[1,1/P(k)]};
        m0 = idpoly(1,B,1,1,F,1,0);
        md = c2d(m0,1);
%         M(:,i)=sim(md,u,'InitialState',prodrate(1));
        
        M_c=sim(md,u,'InitialState',0);
        M(:,k)=M_c;
    end
    
    Depletion_item=fluid_pro_nc(t-Cal_range)*exp((-1*Cal_range)/Ye(1,n,nc)); 
          
    output(1,n)=M(Cal_range,:)*lamd+Depletion_item';                 % output fluid production;  
    
    %calculate oil production
    for i=1:length(lamd)
        Lamd_wcut(i)=lamd(i)*f_wcut(i+1);
    end
    
    output(2,n)=M(Cal_range,:)*Lamd_wcut'+Depletion_item*f_wcut(1);

%      wl=injection(1:t,:)*lamd;
%     
%      Cum_w=sum(wl);
%           
%      if Cum_w<0
%          Cum_w=1;
%      end    
     
%      Q_o=fluid_pro_nc(t)/(1+a*power(Cum_w,b));
% %      output(2,n)=Q_o;                                                % output oil production;
%      Ye(end-2,n)=Cum_w;     
     Ye(end-1:end,n)=output(:,n);
     
end   
     
     
end