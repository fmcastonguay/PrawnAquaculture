clear variables; 
close all;  

% Number of time periods
T= 100; % set max Tstar

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

%Number of collocation points for each case
Nset=60; %60; 

%Solving for steady state (DOES NOT WORK)
if false
%Solving for Steady state
SS=zeros(101,5);
i=0:0.01:1;
if true
for j=1:length(i)
   
[SS(j,:)]=SStatePrawns(i(j), i(j).*K, K, Linf, 1, beta, lambda, gamma, delta, f, K, Linf, kMax, g, aP, bP, muP, d, omega);
    
end
end

%Initial conditions
%Initial_Conditions= [0.384493173606642; 0.535437728657338; 0.999631255086669; 40; 5000]; 
Initial_Conditions= [SS(80,1); SS(80,2); SS(80,3); 40; 2500]; %for SchistoPredation.m %40 mm and 5000 (number, prawns) come from Hoover


x0ic=Initial_Conditions;
end

x0ic=[0.3896,0.5319,0.4681,40,2500]';

%CASE=1; % =1 : No feed; =2 : Supplemental Feed
%OBJ=2;  % =1 : Private Single-shot rotation; =2 : Private Infinite rotation; =3: Societal Single-Rotation; =4: Societal Infinite-Horizon
POLICY=1; % =1: No Intervention; =2 : Minimum Rotation Lenght; =3 : Limiting Rotation Lenght
GUESS=[]; %%% If first run, GUESS=[]


%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sensitivity Analyses
%%%%%%%%%%%%%%%%%%%%%%%%
% Policy Analyses

%Simple Case with no feed
if false
        OBJ=4;  %Health Infinite Horizon
        CASE=1; %No feed

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);



        OBJ=2;  %Private Infinite Horizon
        CASE=1; %No feed

[ts_Private, Topt_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Bs_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, psiNs_Private, alphaNs_Private, Ratios_Private, Ths_Private, ks_Private, Results_Private] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


%Figure
if true
fig=figure
%
subplot(211)
plot(365*ts_No,Profits_No); hold on
plot(365*ts_Private,Profits_Private); hold on
ylim([0 900])
%xlim([0 120])
title({'(A)'},'FontSize', 16)
ylabel({'Aquaculture','Profits'}, 'FontSize', 16);

%
subplot(212)
p1=plot(365*ts_No,Ws_No,'LineWidth',3); hold on
p2=plot(365*ts_Private,Ws_Private,'LineWidth',3); hold on
title({'(B)'},'FontSize', 16)
ylabel({'Infected','Snails'}, 'FontSize', 16);

    legend1=legend([p1 p2],{'Health-Maximization Objective','Profit-Maximization Objective'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.139098711038717 0.814796749318435 0.294643333980015 0.0621428569157917],...
    'Orientation','vertical',...
    'Interpreter','latex');
end

%Table
if true
H=1e6; %Per 1M people

% % % Infinite Horizon
% % No feed

load Workspace_NoFeed.mat 
H=1e6;

%Health
round(365*ts_No(end))
Ns_No(end)*50
round((1/ts_No(end))*H*(Is_No(1)-Is_No(end)))
round(H*(Is_No(1)-Is_No(end)))
round(Profits_No(end))



%Profits
round(365*ts_Private(end))
Ns_Private(end)*50
round((1/ts_Private(end))*H*(Is_Private(1)-Is_Private(end)))
round(H*(Is_Private(1)-Is_Private(end)))
round(Profits_Private(end))


load Workspace_FeedEfficiency.mat 

% % Feed

%Health
round(365*Dynamics_Health_Time(end,2))
Dynamics_Health_Ns(end,2)*50

round((1/Dynamics_Health_Time(end,2))*H*(Dynamics_Health_Is(1,2)-Dynamics_Health_Is(end,2)))
round((H*(Dynamics_Health_Is(1,2)-Dynamics_Health_Is(end,2))))
round(Dynamics_Health_Pi(end,2))


%Profits
round(365*Dynamics_NoPolicy_Time(end,2))
Dynamics_NoPolicy_Ns(end,2)*50
round((1/Dynamics_NoPolicy_Time(end,2))*H*(Dynamics_NoPolicy_Is(1,2)-Dynamics_NoPolicy_Is(end,2)))
round((H*(Dynamics_NoPolicy_Is(1,2)-Dynamics_NoPolicy_Is(end,2))))
round(Dynamics_NoPolicy_Pi(end,2))




end

end



    %New  Policy Analyses: Feed Conversion Effiency
if true
    
    [beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    
    H=1e6;
    
AlphaU=linspace(alphaU*0.5,alphaU*1.5,3); 
%AlphaU=alphaU; 

   
%%% For figure   
Convergence_Time_AlphaU_NoPolicy=zeros(1,length(AlphaU));
Difference_AvoidedInfections_AlphaU_NoPolicy=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_NoPolicy= zeros(1,length(AlphaU));
Average_Cost_AlphaU_NoPolicy=zeros(1,length(AlphaU));

TotalProfits_Convergence_NoPolicy=zeros(1,length(AlphaU));
TotalAvoidedCases_Convergence_NoPolicy=zeros(1,length(AlphaU));

Convergence_Time_AlphaU_RotLen=zeros(1,length(AlphaU));
Difference_AvoidedInfections_AlphaU_RotLen=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_RotLen= zeros(1,length(AlphaU));
Average_Cost_AlphaU_RotLen=zeros(1,length(AlphaU));


Convergence_Time_AlphaU_FeedSea=zeros(1,length(AlphaU));
Difference_AvoidedInfections_AlphaU_FeedSea=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_FeedSea= zeros(1,length(AlphaU));
Average_Cost_AlphaU_FeedSea=zeros(1,length(AlphaU));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(AlphaU)
    
    alphaU=AlphaU(i);

    %%%%%%%%%%%%%%%%%%%
    %%% No Policy
    %%%%%%%%%%%%%%%%%%%
    if true
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   %%%%Providing an Initial Guess     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        OBJ=2;  %Private Infinite Horizon
        %OBJ=1;  %Private Single-Rotation 
        CASE=1; %No feed
        GUESS=[]; %No Guess

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   %%%%Enf of Providing an Initial Guess     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        POLICY=1; %Need to specify Policy in Feed Case, here POLICY=1 means no policy
        OBJ=4; % Health Objective, Infinite Horizon
        %OBJ=3; % Health Objective, Single-Rotation
        CASE=2; %With feed
        
        %Guess
        if true
        GUESS(CASE).I=Is_No;
        GUESS(CASE).W=Ws_No;
        GUESS(CASE).X=Xs_No;
        GUESS(CASE).L=Ls_No;
        GUESS(CASE).P=Ps_No;
        GUESS(CASE).U=0;
        end
    
        [ts_Health, Us_Health, Is_Health, Ws_Health, Xs_Health, Ns_Health, Ls_Health, Ps_Health, Omegas_Health, Profits_Health, psiWs_Health, psiXs_Health, Avoided_HCosts_Health, alphaNs_Health, Ths_Health, ks_Health, Results_Health] = ...
    SchistoAquaculture_Feed(0,0,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);

%Saving Dynamics
if true
Dynamics_Health_Time(:,i)=ts_Health;
Dynamics_Health_Us(:,i)=Us_Health;
Dynamics_Health_Is(:,i)=Is_Health;
Dynamics_Health_Ws(:,i)=Ws_Health;
Dynamics_Health_Xs(:,i)=Xs_Health;
Dynamics_Health_Ns(:,i)=Ns_Health;
Dynamics_Health_Ls(:,i)=Ls_Health;
Dynamics_Health_Ps(:,i)=Ps_Health;
Dynamics_Health_OM(:,i)=Omegas_Health;
Dynamics_Health_Pi(:,i)=Profits_Health;
Dynamics_Health_PsiW(:,i)=psiWs_Health;
Dynamics_Health_Avoided(:,i)=Avoided_HCosts_Health;
Dynamics_Health_alphaNs(:,i)=alphaNs_Health;
Dynamics_Health_Th(:,i)=Ths_Health;
Dynamics_Health_ks(:,i)=ks_Health;
end



  
        OBJ=2; %Private Objective, Infinite Horizon
        %OBJ=1; %Private Objective, Single-Rotation
         
       %Guess       
        if true
        GUESS(CASE).I=Is_Health;
        GUESS(CASE).W=Ws_Health;
        GUESS(CASE).X=Xs_Health;
        GUESS(CASE).L=Ls_Health;
        GUESS(CASE).P=Ps_Health;
        GUESS(CASE).U=Us_Health;
        end

        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Health(end),Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);

%Saving Dynamics
if true
 Dynamics_NoPolicy_Time(:,i)=ts_Private;
 Dynamics_NoPolicy_Us(:,i)=Us_Private;
 Dynamics_NoPolicy_Is(:,i)=Is_Private;
 Dynamics_NoPolicy_Ws(:,i)=Ws_Private;
 Dynamics_NoPolicy_Xs(:,i)=Xs_Private;
 Dynamics_NoPolicy_Ns(:,i)=Ns_Private;
 Dynamics_NoPolicy_Ls(:,i)=Ls_Private;
 Dynamics_NoPolicy_Ps(:,i)=Ps_Private;
 Dynamics_NoPolicy_OM(:,i)=Omegas_Private;
 Dynamics_NoPolicy_Pi(:,i)=Profits_Private;
 Dynamics_NoPolicy_PsiW(:,i)=psiWs_Private;
 Dynamics_NoPolicy_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_NoPolicy_alphaNs(:,i)=alphaNs_Private;
 Dynamics_NoPolicy_Th(:,i)=Ths_Private;
 Dynamics_NoPolicy_ks(:,i)=ks_Private;
end


           %Saving for Figure
    if true
    a=ts_Health(min(find(round(Ws_Health,1)==0))); %Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_AlphaU_NoPolicy(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_NoPolicy(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_NoPolicy(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_AlphaU_NoPolicy(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_AlphaU_NoPolicy(i)=ts_Private(end)./ts_Health(end); %Ratio of rotation lengths
    
    e1=floor(Convergence_Time_AlphaU_NoPolicy(i)./ts_Health(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_AlphaU_NoPolicy(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
    e3= (Convergence_Time_AlphaU_NoPolicy(i)./ts_Health(end) - e1).*ts_Health(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_AlphaU_NoPolicy(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
    e5= Profits_Health(end).*e1 + Profits_Health(max(find(ts_Health<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
    e7= ((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_NoPolicy(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_NoPolicy(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    
    
    % %
    TotalProfits_Convergence_NoPolicy(i)=e6;
    % %
    TotalAvoidedCases_Convergence_NoPolicy(i)=e8;
    
    % %
    Average_Cost_AlphaU_NoPolicy(i) = (e6 - e5)./(e7 - e8); %Average cost of an averted case if we force the health case
    end
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%
    %%% Minimum Rotation Lenght
    %%%%%%%%%%%%%%%%%%%   
    if true
        
        POLICY=2;
        
        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Health(end),Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);


%Saving Dynamics
if true
 Dynamics_RotLen_Time(:,i)=ts_Private;
 Dynamics_RotLen_Us(:,i)=Us_Private;
 Dynamics_RotLen_Is(:,i)=Is_Private;
 Dynamics_RotLen_Ws(:,i)=Ws_Private;
 Dynamics_RotLen_Xs(:,i)=Xs_Private;
 Dynamics_RotLen_Ns(:,i)=Ns_Private;
 Dynamics_RotLen_Ls(:,i)=Ls_Private;
 Dynamics_RotLen_Ps(:,i)=Ps_Private;
 Dynamics_RotLen_OM(:,i)=Omegas_Private;
 Dynamics_RotLen_Pi(:,i)=Profits_Private;
 Dynamics_RotLen_PsiW(:,i)=psiWs_Private;
 Dynamics_RotLen_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_RotLen_alphaNs(:,i)=alphaNs_Private;
 Dynamics_RotLen_Th(:,i)=Ths_Private;
 Dynamics_RotLen_ks(:,i)=ks_Private;
end


    a=ts_Health(min(find(round(Ws_Health,1)==0))) ;%Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_AlphaU_RotLen(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_RotLen(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_RotLen(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_AlphaU_RotLen(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_AlphaU_RotLen(i)=ts_Private(end)./ts_Health(end); %Ratio of rotation lengths
    
   % e1=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Health(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
   % e3= (Convergence_Time_AlphaU_RotLen(i)./ts_Health(end) - e1).*ts_Health(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_AlphaU_RotLen(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
   % e5= Profits_Health(end).*e1 + Profits_Health(max(find(ts_Health<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
   % e7= ((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    

    % %
    Average_Cost_AlphaU_RotLen(i) = (TotalProfits_Convergence_NoPolicy(i) - e6)./(e8 - TotalAvoidedCases_Convergence_NoPolicy(i)); %Average cost of an averted case if we force the health case
    
    
    
        
    end
    

    
    %%%%%%%%%%%%%%%%%%%
    %%% Limiting Feeding Season
    %%%%%%%%%%%%%%%%%%%   
    if true
        
        POLICY=3;
        
        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Health(end),Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);


%Saving Dynamics
if true
 Dynamics_FeedSea_Time(:,i)=ts_Private;
 Dynamics_FeedSea_Us(:,i)=Us_Private;
 Dynamics_FeedSea_Is(:,i)=Is_Private;
 Dynamics_FeedSea_Ws(:,i)=Ws_Private;
 Dynamics_FeedSea_Xs(:,i)=Xs_Private;
 Dynamics_FeedSea_Ns(:,i)=Ns_Private;
 Dynamics_FeedSea_Ls(:,i)=Ls_Private;
 Dynamics_FeedSea_Ps(:,i)=Ps_Private;
 Dynamics_FeedSea_OM(:,i)=Omegas_Private;
 Dynamics_FeedSea_Pi(:,i)=Profits_Private;
 Dynamics_FeedSea_PsiW(:,i)=psiWs_Private;
 Dynamics_FeedSea_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_FeedSea_alphaNs(:,i)=alphaNs_Private;
 Dynamics_FeedSea_Th(:,i)=Ths_Private;
 Dynamics_FeedSea_ks(:,i)=ks_Private;
end


    a=ts_Health(min(find(round(Ws_Health,1)==0))) ;%Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_AlphaU_FeedSea(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_FeedSea(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_FeedSea(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_AlphaU_FeedSea(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_AlphaU_FeedSea(i)=ts_Private(end)./ts_Health(end); %Ratio of rotation lengths
    
   % e1=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Health(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_AlphaU_FeedSea(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
   % e3= (Convergence_Time_AlphaU_RotLen(i)./ts_Health(end) - e1).*ts_Health(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_AlphaU_FeedSea(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
   % e5= Profits_Health(end).*e1 + Profits_Health(max(find(ts_Health<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
   % e7= ((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_FeedSea(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    

    % %
    Average_Cost_AlphaU_FeedSea(i) = (TotalProfits_Convergence_NoPolicy(i) - e6)./(e8 - TotalAvoidedCases_Convergence_NoPolicy(i)); %Average cost of an averted case if we force the health case
    
    
    


    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure for dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
    
set(0, 'DefaultLineLineWidth', 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infect Snails by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(411)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(412)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(413)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(C)  Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(414)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)

    legend1=legend([p1 p2 p3],{'Lower Feed Efficiency','Base Case','Higher Feed Efficiency'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','vertical',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Infected Snails','(as Proportion of Steady State)'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infect Snails by degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(311)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ws(:,1)./x0ic(2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ws(:,1)./x0ic(2),'--'); hold on
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(312)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ws(:,2)./x0ic(2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ws(:,2)./x0ic(2),'--'); hold on
title({'(B) Base Case Feed Efficiency'},'FontSize', 16)
%
subplot(313)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ws(:,3)./x0ic(2),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ws(:,3)./x0ic(2),'--'); hold on
title({'(C) Higher Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Infected Snails','(as Proportion of Steady State)'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Feed by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(411)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Us(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Us(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Us(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(412)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Us(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Us(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Us(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(413)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Us(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Us(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Us(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(414)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Us(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Us(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Us(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Feed by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(311)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Us(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Us(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Us(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Us(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(312)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Us(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Us(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Us(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Us(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(313)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Us(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Us(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Us(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Us(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Aquaculture Profits by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Pi(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Pi(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Pi(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Pi(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Pi(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Pi(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Pi(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Pi(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Pi(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Pi(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Pi(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Pi(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Aquaculture Profits by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Pi(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Pi(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Pi(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Pi(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Pi(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Pi(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Pi(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Pi(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Pi(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Pi(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Pi(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Pi(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end

%%% Extra Dynamics
if false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infection by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Is(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Is(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Is(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Is(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Is(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Is(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Is(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Is(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Is(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Is(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Is(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Is(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infection by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Is(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Is(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Is(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Is(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Is(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Is(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Is(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Is(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Is(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Is(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Is(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Is(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Prawn Predation of Snails by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_PsiW(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_PsiW(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_PsiW(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_PsiW(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_PsiW(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_PsiW(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_PsiW(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_PsiW(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_PsiW(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_PsiW(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_PsiW(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_PsiW(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Prawn Predation by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_PsiW(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_PsiW(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_PsiW(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_PsiW(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_PsiW(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_PsiW(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_PsiW(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_PsiW(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_PsiW(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_PsiW(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_PsiW(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_PsiW(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Number of Prawns by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ps(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ps(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ps(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ps(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ps(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ps(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ps(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ps(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ps(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ps(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ps(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ps(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Number of Prawns by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ps(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ps(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ps(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ps(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ps(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ps(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ps(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ps(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ps(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ps(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ps(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ps(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end


end


    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure for policy evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(231)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
ylim([0 50])
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_NoPolicy,'LineWidth',3); hold on
ylim([0 2.5])
%ylim([0.3 1.05])
%plot([alphaU alphaU],[0 1.05], 'Color', 'k')

subplot(232)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([alphaU alphaU],[0 6.1], 'Color', 'k')
ylim([0 50])
%plot([alphaU alphaU],[0 40], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_RotLen,'LineWidth',3); hold on
ylim([0 2.5])
%ylim([0.3 1.05])

subplot(233)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
plot([alphaU alphaU],[0 6.1], 'Color', 'k')
ylim([0 50])
%plot([alphaU alphaU],[0 40], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_FeedSea,'LineWidth',3); hold on
ylim([0 2.5])
%ylim([0.3 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
%yyaxis left
plot(AlphaU,Difference_AvoidedInfections_AlphaU_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
%plot([alphaU alphaU],[0 10000], 'Color', 'k')
%
%yyaxis right
%plot(AlphaU,Average_Cost_AlphaU_NoPolicy,'LineWidth',3);
%ylim([0 230])
%ylim([0 375])

subplot(235)
%yyaxis left
plot(AlphaU,Difference_AvoidedInfections_AlphaU_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
%plot([alphaU alphaU],[0 10000], 'Color', 'k')
%
%yyaxis right
%plot(AlphaU,Average_Profits_RotLen,'LineWidth',3); hold on
%plot(AlphaU,Average_Cost_AlphaU_RotLen,'LineWidth',3)
%ylim([0 230])
%ylim([0 375])
%
subplot(236)
%yyaxis left
plot(AlphaU,Difference_AvoidedInfections_AlphaU_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
%plot([alphaU alphaU],[0 10000], 'Color', 'k')
%plot([alphaU alphaU],[0 25000], 'Color', 'k')
%
%yyaxis right
%plot(AlphaU,Average_Profits_FeedSea,'LineWidth',3); hold on
%plot(AlphaU,Average_Cost_AlphaU_FeedSea,'LineWidth',3)
%ylim([0 230])
%ylim([0 375])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
% saveas(gcf,'Attack_Feed_025.png'); hold off    
end

    

end


    %New  Policy Analyses: Holling Exponent for Type III Functional Response
if false
    
    [beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    
    H=1e6;
    
N_Holling=linspace(1.75,2.25,3); 

   
%%% For figure   
Convergence_Time_Holling_NoPolicy=zeros(1,length(N_Holling));
Difference_AvoidedInfections_Holling_NoPolicy=zeros(1,length(N_Holling));
Rotation_Length_Ratio_Holling_NoPolicy= zeros(1,length(N_Holling));
Average_Cost_Holling_NoPolicy=zeros(1,length(N_Holling));

TotalProfits_Convergence_NoPolicy=zeros(1,length(N_Holling));
TotalAvoidedCases_Convergence_NoPolicy=zeros(1,length(N_Holling));

Convergence_Time_Holling_RotLen=zeros(1,length(N_Holling));
Difference_AvoidedInfections_Holling_RotLen=zeros(1,length(N_Holling));
Rotation_Length_Ratio_Holling_RotLen= zeros(1,length(N_Holling));
Average_Cost_Holling_RotLen=zeros(1,length(N_Holling));


Convergence_Time_Holling_FeedSea=zeros(1,length(N_Holling));
Difference_AvoidedInfections_Holling_FeedSea=zeros(1,length(N_Holling));
Rotation_Length_Ratio_Holling_FeedSea= zeros(1,length(N_Holling));
Average_Cost_Holling_FeedSea=zeros(1,length(N_Holling));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Running the Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(N_Holling)
    
    n=N_Holling(i);

    %%%%%%%%%%%%%%%%%%%
    %%% No Policy
    %%%%%%%%%%%%%%%%%%%
    if true
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   %%%%Providing an Initial Guess     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        OBJ=2;  %Private Infinite Horizon
        CASE=1; %No feed
        GUESS=[]; %No Guess

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   %%%%Enf of Providing an Initial Guess     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        POLICY=1; %Need to specify Policy in Feed Case, here POLICY=1 means no policy
        OBJ=4; % Health Objective, Infinite Horizon
        CASE=2; %With feed
        
        %Guess
        if true
        GUESS(CASE).I=Is_No;
        GUESS(CASE).W=Ws_No;
        GUESS(CASE).X=Xs_No;
        GUESS(CASE).L=Ls_No;
        GUESS(CASE).P=Ps_No;
        GUESS(CASE).U=0;
        end
    
        [ts_Health, Us_Health, Is_Health, Ws_Health, Xs_Health, Ns_Health, Ls_Health, Ps_Health, Omegas_Health, Profits_Health, psiWs_Health, psiXs_Health, Avoided_HCosts_Health, alphaNs_Health, Ths_Health, ks_Health, Results_Health] = ...
    SchistoAquaculture_Feed(0,0,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);

%Saving Dynamics
if true
Dynamics_Health_Time(:,i)=ts_Health;
Dynamics_Health_Us(:,i)=Us_Health;
Dynamics_Health_Is(:,i)=Is_Health;
Dynamics_Health_Ws(:,i)=Ws_Health;
Dynamics_Health_Xs(:,i)=Xs_Health;
Dynamics_Health_Ns(:,i)=Ns_Health;
Dynamics_Health_Ls(:,i)=Ls_Health;
Dynamics_Health_Ps(:,i)=Ps_Health;
Dynamics_Health_OM(:,i)=Omegas_Health;
Dynamics_Health_Pi(:,i)=Profits_Health;
Dynamics_Health_PsiW(:,i)=psiWs_Health;
Dynamics_Health_Avoided(:,i)=Avoided_HCosts_Health;
Dynamics_Health_alphaNs(:,i)=alphaNs_Health;
Dynamics_Health_Th(:,i)=Ths_Health;
Dynamics_Health_ks(:,i)=ks_Health;
end



  
        OBJ=2; %Private Objective, Infinite Horizon
       %Guess       
        if true
        GUESS(CASE).I=Is_Health;
        GUESS(CASE).W=Ws_Health;
        GUESS(CASE).X=Xs_Health;
        GUESS(CASE).L=Ls_Health;
        GUESS(CASE).P=Ps_Health;
        GUESS(CASE).U=Us_Health;
        end

        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Health(end),Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);

%Saving Dynamics
if true
 Dynamics_NoPolicy_Time(:,i)=ts_Private;
 Dynamics_NoPolicy_Us(:,i)=Us_Private;
 Dynamics_NoPolicy_Is(:,i)=Is_Private;
 Dynamics_NoPolicy_Ws(:,i)=Ws_Private;
 Dynamics_NoPolicy_Xs(:,i)=Xs_Private;
 Dynamics_NoPolicy_Ns(:,i)=Ns_Private;
 Dynamics_NoPolicy_Ls(:,i)=Ls_Private;
 Dynamics_NoPolicy_Ps(:,i)=Ps_Private;
 Dynamics_NoPolicy_OM(:,i)=Omegas_Private;
 Dynamics_NoPolicy_Pi(:,i)=Profits_Private;
 Dynamics_NoPolicy_PsiW(:,i)=psiWs_Private;
 Dynamics_NoPolicy_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_NoPolicy_alphaNs(:,i)=alphaNs_Private;
 Dynamics_NoPolicy_Th(:,i)=Ths_Private;
 Dynamics_NoPolicy_ks(:,i)=ks_Private;
end


           %Saving for Figure
        if true
    a=ts_Health(min(find(round(Ws_Health,1)==0))); %Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_Holling_NoPolicy(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_Holling_NoPolicy(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_Holling_NoPolicy(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_Holling_NoPolicy(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_Holling_NoPolicy(i)=ts_Private(end)./ts_Health(end); %Ratio of rotation lengths
    
    e1=floor(Convergence_Time_Holling_NoPolicy(i)./ts_Health(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_Holling_NoPolicy(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
    e3= (Convergence_Time_Holling_NoPolicy(i)./ts_Health(end) - e1).*ts_Health(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_Holling_NoPolicy(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
    e5= Profits_Health(end).*e1 + Profits_Health(max(find(ts_Health<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
    e7= ((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_Holling_NoPolicy(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_Holling_NoPolicy(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    
    
    % %
    TotalProfits_Convergence_NoPolicy(i)=e6;
    % %
    TotalAvoidedCases_Convergence_NoPolicy(i)=e8;
    
    % %
    Average_Cost_Holling_NoPolicy(i) = (e6 - e5)./(e7 - e8); %Average cost of an averted case if we force the health case
        end
    
    end
    
    
        %%%%%%%%%%%%%%%%%%%
    %%% Minimum Rotation Lenght
    %%%%%%%%%%%%%%%%%%%   
    if true
        
        POLICY=2;
        
        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Health(end),Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);


%Saving Dynamics
if true
 Dynamics_RotLen_Time(:,i)=ts_Private;
 Dynamics_RotLen_Us(:,i)=Us_Private;
 Dynamics_RotLen_Is(:,i)=Is_Private;
 Dynamics_RotLen_Ws(:,i)=Ws_Private;
 Dynamics_RotLen_Xs(:,i)=Xs_Private;
 Dynamics_RotLen_Ns(:,i)=Ns_Private;
 Dynamics_RotLen_Ls(:,i)=Ls_Private;
 Dynamics_RotLen_Ps(:,i)=Ps_Private;
 Dynamics_RotLen_OM(:,i)=Omegas_Private;
 Dynamics_RotLen_Pi(:,i)=Profits_Private;
 Dynamics_RotLen_PsiW(:,i)=psiWs_Private;
 Dynamics_RotLen_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_RotLen_alphaNs(:,i)=alphaNs_Private;
 Dynamics_RotLen_Th(:,i)=Ths_Private;
 Dynamics_RotLen_ks(:,i)=ks_Private;
end


    a=ts_Health(min(find(round(Ws_Health,1)==0))) ;%Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_Holling_RotLen(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_Holling_RotLen(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_Holling_RotLen(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_Holling_RotLen(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_Holling_RotLen(i)=ts_Private(end)./ts_Health(end); %Ratio of rotation lengths
    
   % e1=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Health(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_Holling_RotLen(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
   % e3= (Convergence_Time_AlphaU_RotLen(i)./ts_Health(end) - e1).*ts_Health(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_Holling_RotLen(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
   % e5= Profits_Health(end).*e1 + Profits_Health(max(find(ts_Health<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
   % e7= ((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_Holling_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    

    % %
    Average_Cost_Holling_RotLen(i) = (TotalProfits_Convergence_NoPolicy(i) - e6)./(e8 - TotalAvoidedCases_Convergence_NoPolicy(i)); %Average cost of an averted case if we force the health case
    
    
    
        
    end
    

    
    %%%%%%%%%%%%%%%%%%%
    %%% Limiting Feeding Season
    %%%%%%%%%%%%%%%%%%%   
    if true
        
        
        
               %Guess       
        if true
        GUESS(CASE).I=Is_Private;
        GUESS(CASE).W=Ws_Private;
        GUESS(CASE).X=Xs_Private;
        GUESS(CASE).L=Ls_Private;
        GUESS(CASE).P=Ps_Private;
        GUESS(CASE).U=Us_Private;
        end

        
        POLICY=3;
        
        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Health(end),Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);


%Saving Dynamics
if true
 Dynamics_FeedSea_Time(:,i)=ts_Private;
 Dynamics_FeedSea_Us(:,i)=Us_Private;
 Dynamics_FeedSea_Is(:,i)=Is_Private;
 Dynamics_FeedSea_Ws(:,i)=Ws_Private;
 Dynamics_FeedSea_Xs(:,i)=Xs_Private;
 Dynamics_FeedSea_Ns(:,i)=Ns_Private;
 Dynamics_FeedSea_Ls(:,i)=Ls_Private;
 Dynamics_FeedSea_Ps(:,i)=Ps_Private;
 Dynamics_FeedSea_OM(:,i)=Omegas_Private;
 Dynamics_FeedSea_Pi(:,i)=Profits_Private;
 Dynamics_FeedSea_PsiW(:,i)=psiWs_Private;
 Dynamics_FeedSea_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_FeedSea_alphaNs(:,i)=alphaNs_Private;
 Dynamics_FeedSea_Th(:,i)=Ths_Private;
 Dynamics_FeedSea_ks(:,i)=ks_Private;
end


    a=ts_Health(min(find(round(Ws_Health,1)==0))) ;%Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_Holling_FeedSea(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_Holling_FeedSea(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_Holling_FeedSea(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_Holling_FeedSea(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_Holling_FeedSea(i)=ts_Private(end)./ts_Health(end); %Ratio of rotation lengths
    
   % e1=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Health(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_Holling_FeedSea(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
   % e3= (Convergence_Time_AlphaU_RotLen(i)./ts_Health(end) - e1).*ts_Health(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_Holling_FeedSea(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
   % e5= Profits_Health(end).*e1 + Profits_Health(max(find(ts_Health<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
   % e7= ((Is_Health(1)-Is_Health(end))./ts_Health(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_Holling_FeedSea(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    

    % %
    Average_Cost_Holling_FeedSea(i) = (TotalProfits_Convergence_NoPolicy(i) - e6)./(e8 - TotalAvoidedCases_Convergence_NoPolicy(i)); %Average cost of an averted case if we force the health case
    
    
    


    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure for dynamics  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infect Snails by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(411)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(412)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(413)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(C)  Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(414)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)

    legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Infected Snails','(as Proportion of Steady State)'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infect Snails by degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(311)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ws(:,1)./x0ic(2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ws(:,1)./x0ic(2),'--'); hold on
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(312)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ws(:,2)./x0ic(2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ws(:,2)./x0ic(2),'--'); hold on
title({'(B) Base Case Feed Efficiency'},'FontSize', 16)
%
subplot(313)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ws(:,3)./x0ic(2),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ws(:,3)./x0ic(2),'--'); hold on
title({'(C) Higher Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Infected Snails','(as Proportion of Steady State)'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Feed by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(411)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Us(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Us(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Us(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(412)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Us(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Us(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Us(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(413)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Us(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Us(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Us(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(414)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Us(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Us(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Us(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Feed by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(311)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Us(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Us(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Us(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Us(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(312)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Us(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Us(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Us(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Us(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(313)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Us(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Us(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Us(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Us(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Aquaculture Profits by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Pi(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Pi(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Pi(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Pi(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Pi(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Pi(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Pi(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Pi(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Pi(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Pi(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Pi(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Pi(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Aquaculture Profits by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Pi(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Pi(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Pi(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Pi(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Pi(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Pi(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Pi(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Pi(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Pi(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Pi(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Pi(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Pi(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end

    

if false
    
    
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infection by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Is(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Is(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Is(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Is(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Is(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Is(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Is(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Is(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Is(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Is(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Is(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Is(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infection by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Is(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Is(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Is(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Is(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Is(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Is(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Is(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Is(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Is(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Is(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Is(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Is(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Prawn Predation of Snails by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_PsiW(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_PsiW(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_PsiW(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_PsiW(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_PsiW(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_PsiW(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_PsiW(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_PsiW(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_PsiW(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_PsiW(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_PsiW(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_PsiW(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Prawn Predation by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_PsiW(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_PsiW(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_PsiW(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_PsiW(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_PsiW(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_PsiW(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_PsiW(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_PsiW(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_PsiW(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_PsiW(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_PsiW(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_PsiW(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Number of Prawns by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ps(:,1)); hold on
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ps(:,2)); hold on
plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ps(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ps(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ps(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ps(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ps(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ps(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ps(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ps(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ps(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ps(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Number of Prawns by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Health_Time(:,1),Dynamics_Health_Ps(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ps(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ps(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ps(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Health_Time(:,2),Dynamics_Health_Ps(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ps(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ps(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ps(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Health_Time(:,3),Dynamics_Health_Ps(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ps(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ps(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ps(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end



end

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure for policy evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Policy Evaluation Figure
if true
fig=figure
subplot(231)
yyaxis left
plot(N_Holling,Convergence_Time_Holling_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
%
yyaxis right
plot(N_Holling,Rotation_Length_Ratio_Holling_NoPolicy,'LineWidth',3); hold on
%ylim([0.45 1.05])
%ylim([0.3 1.05])
%plot([N_Holling(2) N_Holling(2)],[0 1.05], 'Color', 'k')

subplot(232)
yyaxis left
plot(N_Holling,Convergence_Time_Holling_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([N_Holling(2) N_Holling(2)],[0 6.1], 'Color', 'k')
%plot([N_Holling(2) N_Holling(2)],[0 40], 'Color', 'k')
%
yyaxis right
plot(N_Holling,Rotation_Length_Ratio_Holling_RotLen,'LineWidth',3); hold on
%ylim([0.45 1.05])
%ylim([0.3 1.05])

subplot(233)
yyaxis left
plot(N_Holling,Convergence_Time_Holling_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
%plot([N_Holling(2) N_Holling(2)],[0 6.1], 'Color', 'k')
%plot([N_Holling(2) N_Holling(2)],[0 40], 'Color', 'k')
%
yyaxis right
plot(N_Holling,Rotation_Length_Ratio_Holling_FeedSea,'LineWidth',3); hold on
%ylim([0.45 1.05])
%ylim([0.3 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
%yyaxis left
plot(N_Holling,Difference_AvoidedInfections_Holling_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
%plot([N_Holling(2) N_Holling(2)],[0 3000], 'Color', 'k')
%plot([N_Holling(2) N_Holling(2)],[0 10000], 'Color', 'k')
%
%yyaxis right
%plot(N_Holling,Average_Cost_Holling_NoPolicy,'LineWidth',3);
%ylim([0 230])
%ylim([0 375])

subplot(235)
%yyaxis left
plot(N_Holling,Difference_AvoidedInfections_Holling_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
%plot([N_Holling(2) N_Holling(2)],[0 3000], 'Color', 'k')
%plot([N_Holling(2) N_Holling(2)],[0 10000], 'Color', 'k')
%
%yyaxis right
%plot(N_Holling,Average_Profits_RotLen,'LineWidth',3); hold on
%plot(N_Holling,Average_Cost_Holling_RotLen,'LineWidth',3)
%ylim([0 230])
%ylim([0 375])
%
subplot(236)
%yyaxis left
plot(N_Holling,Difference_AvoidedInfections_Holling_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
%plot([N_Holling(2) N_Holling(2)],[0 10000], 'Color', 'k')
%plot([N_Holling(2) N_Holling(2)],[0 25000], 'Color', 'k')
%
%yyaxis right
%plot(N_Holling,Average_Profits_FeedSea,'LineWidth',3); hold on
%plot(N_Holling,Average_Cost_Holling_FeedSea,'LineWidth',3)
%ylim([0 230])
%ylim([0 375])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
  %  title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
% saveas(gcf,'Attack_Feed_025.png'); hold off    
end

    


end





















%Dynamics of the policies
if false

    
    
for Nset=50
%Attack rate on supplemental feed
if true

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    AlphaU=linspace(alphaU*0.5,alphaU*1.5,3);
    
    Dynamics_L_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Opt_AlphaU=zeros(Nset+1,length(AlphaU));
    
    
    Dynamics_L_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_P_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_W_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_NoFeed_AlphaU=zeros(Nset+1,length(AlphaU));

    
    Dynamics_L_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_AlphaU=zeros(Nset+1,length(AlphaU));
    
    
    Dynamics_L_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_X_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_I_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));

    

    Dynamics_L_Feed_FeedSea_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed_FeedSea_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed_FeedSea_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_FeedSea_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_FeedSea_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_X_Feed_FeedSea_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_I_Feed_FeedSea_AlphaU=zeros(Nset+1,length(AlphaU));



    
        for i=1:length(AlphaU)
        
        alphaU=AlphaU(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        
        OBJ=3;
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_Opt(:,i)=Ls_No;
    Dynamics_P_Opt(:,i)=Ps_No;
    Dynamics_W_Opt(:,i)=Ws_No;
    Dynamics_Time_Opt_AlphaU(:,i)=ts_No;
    
    
        OBJ=1;
        CASE=1;

[ts_NoFeed, Topt_NoFeed, Is_NoFeed, Ws_NoFeed, Xs_NoFeed, Ns_NoFeed, Ls_NoFeed, Ps_NoFeed, Bs_NoFeed, Omegas_NoFeed, Profits_NoFeed, psiWs_NoFeed, psiXs_NoFeed, psiNs_NoFeed, alphaNs_NoFeed, Ratios_NoFeed, Ths_NoFeed, ks_NoFeed, Results_NoFeed] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_NoFeed(:,i)=Ls_NoFeed;
    Dynamics_P_NoFeed(:,i)=Ps_NoFeed;
    Dynamics_W_NoFeed(:,i)=Ws_NoFeed;
    Dynamics_Time_NoFeed_AlphaU(:,i)=ts_NoFeed;

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed(:,i)=Ls;
    Dynamics_P_Feed(:,i)=Ps;
    Dynamics_W_Feed(:,i)=Ws;
    Dynamics_U_Feed_AlphaU(:,i)=Us;
    Dynamics_Time_Feed_AlphaU(:,i)=ts;
    
        end   
        
        
        
        %%% Standardized Rotation Length
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        
       
        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_RotLen(:,i)=Ls;
    Dynamics_P_Feed_RotLen(:,i)=Ps;
    Dynamics_W_Feed_RotLen(:,i)=Ws;
    Dynamics_U_Feed_RotLen_AlphaU(:,i)=Us;
    Dynamics_Time_Feed_RotLen_AlphaU(:,i)=ts;
    Dynamics_X_Feed_RotLen_AlphaU(:,i)=Xs;
    Dynamics_I_Feed_RotLen_AlphaU(:,i)=Is;
    
        end   
      
        
        
        %%% Limited Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=4;
        
     
        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_FeedSea_AlphaU(:,i)=Ls;
    Dynamics_P_Feed_FeedSea_AlphaU(:,i)=Ps;
    Dynamics_W_Feed_FeedSea_AlphaU(:,i)=Ws;
    Dynamics_U_Feed_FeedSea_AlphaU(:,i)=Us;
    Dynamics_Time_Feed_FeedSea_AlphaU(:,i)=ts;
    Dynamics_X_Feed_FeedSea_AlphaU(:,i)=Xs;
    Dynamics_I_Feed_FeedSea_AlphaU(:,i)=Is;
        end   

        
        
        end
    
       %%%figure 
if true
    
    
    %Colors
    bC=[0, 0.4470, 0.7410];
    rC=[0.8500, 0.3250, 0.0980];
    yC=[0.9290, 0.6940, 0.1250];
        
       fig= figure
     
        subplot(331)      
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,1)*365,Dynamics_W_Feed(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1]) 
        xlim([0 300])
        title({'No Policy','(A)'},'FontSize', 16)
        ylabel({'Infected Snails','(Low Feed Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'
        
        subplot(334)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,2)*365,Dynamics_W_Feed(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        ylabel({'Infected Snails','(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(D)'},'FontSize', 16)
        
        subplot(337)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,3)*365,Dynamics_W_Feed(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        ylabel({'Infected Snails','(High Feed Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(G)'},'FontSize', 16)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(332)      
        p1=plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        p2=plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        p3=plot(Dynamics_Time_Feed_RotLen_AlphaU(:,1)*365,Dynamics_W_Feed_RotLen(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1]) 
        xlim([0 300])
        title({'Minimum Length','(B)'},'FontSize', 16)
        legend1=legend([p1 p2 p3],{'Societal Optimum ~~~','Aquaculture Optimum w/o Feed~~~','Aquaculture Optimum w/ Feed'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
        set(legend1,...
    'Position',[0.0938584639410563 0.00827803823612821 0.837143993377685 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

        
        subplot(335)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen_AlphaU(:,2)*365,Dynamics_W_Feed_RotLen(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(E)'},'FontSize', 16)
         
         
        subplot(338)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen_AlphaU(:,3)*365,Dynamics_W_Feed_RotLen(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(H)'},'FontSize', 16)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(333)      
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU(:,1)*365,Dynamics_W_Feed_FeedSea_AlphaU(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'Limited Feeding','(C)'},'FontSize', 16)
               
        subplot(336)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU(:,2)*365,Dynamics_W_Feed_FeedSea_AlphaU(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(F)'},'FontSize', 16)
        
        
        subplot(339)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU(:,3)*365,Dynamics_W_Feed_FeedSea_AlphaU(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(I)'},'FontSize', 16)
        

            han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Time (days)',''}, 'FontSize', 16);
    
    
   % saveas(gcf,'SA_DynamicsW_AttackFeed.png'); hold off 
        
        
        
        
        
        
end
    
end  
    
%Holling's Type III exponent
if true

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    N=linspace(1.75,2.25,3);
    
    Dynamics_L_Opt=zeros(Nset+1,length(N));
    Dynamics_P_Opt=zeros(Nset+1,length(N));
    Dynamics_W_Opt=zeros(Nset+1,length(N));
    Dynamics_Time_Opt=zeros(Nset+1,length(N));

    Dynamics_L_NoFeed=zeros(Nset+1,length(N));
    Dynamics_P_NoFeed=zeros(Nset+1,length(N));
    Dynamics_W_NoFeed=zeros(Nset+1,length(N));
    Dynamics_Time_NoFeed=zeros(Nset+1,length(N));
    
    
    Dynamics_L_Feed=zeros(Nset+1,length(N));
    Dynamics_P_Feed=zeros(Nset+1,length(N));
    Dynamics_W_Feed=zeros(Nset+1,length(N));
    Dynamics_U_Feed=zeros(Nset+1,length(N));
    Dynamics_Time_Feed=zeros(Nset+1,length(N));
        

    Dynamics_L_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_P_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_W_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_U_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_Time_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_X_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_I_Feed_RotLen=zeros(Nset+1,length(N));
    
    Dynamics_L_Feed_FeedSea=zeros(Nset+1,length(N));
    Dynamics_P_Feed_FeedSea=zeros(Nset+1,length(N));
    Dynamics_W_Feed_FeedSea=zeros(Nset+1,length(N));
    Dynamics_U_Feed_FeedSea=zeros(Nset+1,length(N));
    Dynamics_Time_Feed_FeedSea=zeros(Nset+1,length(N));
    Dynamics_X_Feed_FeedSea=zeros(Nset+1,length(N));
    Dynamics_I_Feed_FeedSea=zeros(Nset+1,length(N));


    
        for i=1:length(N)
        
        n=N(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        
        OBJ=3;
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_Opt(:,i)=Ls_No;
    Dynamics_P_Opt(:,i)=Ps_No;
    Dynamics_W_Opt(:,i)=Ws_No;
    Dynamics_Time_Opt(:,i)=ts_No;
    
    
        OBJ=1;
        CASE=1;

[ts_NoFeed, Topt_NoFeed, Is_NoFeed, Ws_NoFeed, Xs_NoFeed, Ns_NoFeed, Ls_NoFeed, Ps_NoFeed, Bs_NoFeed, Omegas_NoFeed, Profits_NoFeed, psiWs_NoFeed, psiXs_NoFeed, psiNs_NoFeed, alphaNs_NoFeed, Ratios_NoFeed, Ths_NoFeed, ks_NoFeed, Results_NoFeed] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_NoFeed(:,i)=Ls_NoFeed;
    Dynamics_P_NoFeed(:,i)=Ps_NoFeed;
    Dynamics_W_NoFeed(:,i)=Ws_NoFeed;
    Dynamics_Time_NoFeed(:,i)=ts_NoFeed;

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed(:,i)=Ls;
    Dynamics_P_Feed(:,i)=Ps;
    Dynamics_W_Feed(:,i)=Ws;
    Dynamics_U_Feed(:,i)=Us;
    Dynamics_Time_Feed(:,i)=ts;
    
        end   
        
        
        
        %%% Standardized Rotation Length
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_RotLen(:,i)=Ls;
    Dynamics_P_Feed_RotLen(:,i)=Ps;
    Dynamics_W_Feed_RotLen(:,i)=Ws;
    Dynamics_U_Feed_RotLen(:,i)=Us;
    Dynamics_Time_Feed_RotLen(:,i)=ts;
    Dynamics_X_Feed_RotLen(:,i)=Xs;
    Dynamics_I_Feed_RotLen(:,i)=Is;
    
        end   
      
        
        
        %%% Limited Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=4;
        
        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_FeedSea(:,i)=Ls;
    Dynamics_P_Feed_FeedSea(:,i)=Ps;
    Dynamics_W_Feed_FeedSea(:,i)=Ws;
    Dynamics_U_Feed_FeedSea(:,i)=Us;
    Dynamics_Time_Feed_FeedSea(:,i)=ts;
    Dynamics_X_Feed_FeedSea(:,i)=Xs;
    Dynamics_I_Feed_FeedSea(:,i)=Is;
        end   

        
        
        end
    

    
        
        
       %%%figure 
if true
    
    
    %Colors
    bC=[0, 0.4470, 0.7410];
    rC=[0.8500, 0.3250, 0.0980];
    yC=[0.9290, 0.6940, 0.1250];
        
       fig= figure
     
        subplot(331)      
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,1)*365,Dynamics_W_Feed(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1]) 
        xlim([0 300])
        title({'No Policy','(A)'},'FontSize', 16)
        ylabel({'Infected Snails','(Specialized Prawns)'}, 'FontSize', 14) %,'interpreter','latex'
        
        subplot(334)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,2)*365,Dynamics_W_Feed(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        ylabel({'Infected Snails','(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(D)'},'FontSize', 16)
        
        subplot(337)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,3)*365,Dynamics_W_Feed(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        ylabel({'Infected Snails','(Generalist Prawns)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(G)'},'FontSize', 16)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(332)      
        p1=plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        p2=plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        p3=plot(Dynamics_Time_Feed_RotLen(:,1)*365,Dynamics_W_Feed_RotLen(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1]) 
        xlim([0 300])
        title({'Minimum Length','(B)'},'FontSize', 16)
        legend1=legend([p1 p2 p3],{'Societal Optimum ~~~','Aquaculture Optimum w/o Feed~~~','Aquaculture Optimum w/ Feed'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
        set(legend1,...
    'Position',[0.0938584639410563 0.00827803823612821 0.837143993377685 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

        
        subplot(335)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen(:,2)*365,Dynamics_W_Feed_RotLen(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(E)'},'FontSize', 16)
         
         
        subplot(338)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen(:,3)*365,Dynamics_W_Feed_RotLen(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(H)'},'FontSize', 16)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(333)      
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea(:,1)*365,Dynamics_W_Feed_FeedSea(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1]) 
        xlim([0 300])
        title({'Limited Feeding','(C)'},'FontSize', 16)
               
        subplot(336)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea(:,2)*365,Dynamics_W_Feed_FeedSea(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(F)'},'FontSize', 16)
        
        
        subplot(339)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea(:,3)*365,Dynamics_W_Feed_FeedSea(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(I)'},'FontSize', 16)
        

            han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    %xlabel(han,{'','','','Time (days)'}, 'FontSize', 16);
    xlabel(han,{'Time (days)',''}, 'FontSize', 16);
        
        
   %saveas(gcf,'SA_DynamicsW_Holling.png'); hold off 
        
        
        %figure for feed
        
    
        
end
    
        
        
end


end




%For Nset=60
if true
for Nset=60
    
    %Attack rate on supplemental feed
if true

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    AlphaU=linspace(alphaU*0.5,alphaU*1.5,3);
    
    Dynamics_L_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Opt_AlphaU=zeros(Nset+1,length(AlphaU));
    
    
    Dynamics_L_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_P_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_W_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_NoFeed_AlphaU=zeros(Nset+1,length(AlphaU));

    
    Dynamics_L_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_AlphaU=zeros(Nset+1,length(AlphaU));
    
    
    Dynamics_L_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));
    

    Dynamics_L_Feed_FeedSea_AlphaU2=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed_FeedSea_AlphaU2=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed_FeedSea_AlphaU2=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_FeedSea_AlphaU2=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_FeedSea_AlphaU2=zeros(Nset+1,length(AlphaU));
    Dynamics_X_Feed_FeedSea_AlphaU2=zeros(Nset+1,length(AlphaU));
    Dynamics_I_Feed_FeedSea_AlphaU2=zeros(Nset+1,length(AlphaU));

    
        for i=1:length(AlphaU)
        
        alphaU=AlphaU(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        GUESS=[];
        OBJ=3;
        CASE=1;
        
[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_Opt(:,i)=Ls_No;
    Dynamics_P_Opt(:,i)=Ps_No;
    Dynamics_W_Opt(:,i)=Ws_No;
    Dynamics_Time_Opt_AlphaU(:,i)=ts_No;
    
    
        OBJ=1;
        CASE=1;

[ts_NoFeed, Topt_NoFeed, Is_NoFeed, Ws_NoFeed, Xs_NoFeed, Ns_NoFeed, Ls_NoFeed, Ps_NoFeed, Bs_NoFeed, Omegas_NoFeed, Profits_NoFeed, psiWs_NoFeed, psiXs_NoFeed, psiNs_NoFeed, alphaNs_NoFeed, Ratios_NoFeed, Ths_NoFeed, ks_NoFeed, Results_NoFeed] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_NoFeed(:,i)=Ls_NoFeed;
    Dynamics_P_NoFeed(:,i)=Ps_NoFeed;
    Dynamics_W_NoFeed(:,i)=Ws_NoFeed;
    Dynamics_Time_NoFeed_AlphaU(:,i)=ts_NoFeed;

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed(:,i)=Ls;
    Dynamics_P_Feed(:,i)=Ps;
    Dynamics_W_Feed(:,i)=Ws;
    Dynamics_U_Feed_AlphaU(:,i)=Us;
    Dynamics_Time_Feed_AlphaU(:,i)=ts;
    
        end   
        
        
        
        %%% Standardized Rotation Length
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        
       
        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_RotLen(:,i)=Ls;
    Dynamics_P_Feed_RotLen(:,i)=Ps;
    Dynamics_W_Feed_RotLen(:,i)=Ws;
    Dynamics_U_Feed_RotLen_AlphaU(:,i)=Us;
    Dynamics_Time_Feed_RotLen_AlphaU(:,i)=ts;
    
        end   
      
        
        
        %%% Limited Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        
      
            GUESS(CASE).I=Dynamics_I_Feed_FeedSea_AlphaU(:,1);
            GUESS(CASE).W=Dynamics_W_Feed_FeedSea_AlphaU(:,1);
            GUESS(CASE).X=Dynamics_X_Feed_FeedSea_AlphaU(:,1);
            GUESS(CASE).L=Dynamics_L_Feed_FeedSea_AlphaU(:,1);
            GUESS(CASE).P=Dynamics_P_Feed_FeedSea_AlphaU(:,1);
            GUESS(CASE).U=Dynamics_U_Feed_FeedSea_AlphaU(:,1);
    

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_FeedSea_AlphaU2(:,i)=Ls;
    Dynamics_P_Feed_FeedSea_AlphaU2(:,i)=Ps;
    Dynamics_W_Feed_FeedSea_AlphaU2(:,i)=Ws;
    Dynamics_U_Feed_FeedSea_AlphaU2(:,i)=Us;
    Dynamics_Time_Feed_FeedSea_AlphaU2(:,i)=ts;
    Dynamics_X_Feed_FeedSea_AlphaU2(:,i)=Xs;
    Dynamics_I_Feed_FeedSea_AlphaU2(:,i)=Is;

        end   

        
        
        end
    
       %%%figure 
if true
    
    
    %Colors
    bC=[0, 0.4470, 0.7410];
    rC=[0.8500, 0.3250, 0.0980];
    yC=[0.9290, 0.6940, 0.1250];
        
       fig= figure
     
        subplot(331)      
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,1)*365,Dynamics_W_Feed(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1]) 
        xlim([0 300]) 
        title({'No Policy','(A)'},'FontSize', 16)
        ylabel({'Infected Snails','(Low Feed Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'
        
        subplot(334)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,2)*365,Dynamics_W_Feed(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        ylabel({'Infected Snails','(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(D)'},'FontSize', 16)
        
        subplot(337)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,3)*365,Dynamics_W_Feed(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        ylabel({'Infected Snails','(High Feed Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(G)'},'FontSize', 16)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(332)      
        p1=plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        p2=plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        p3=plot(Dynamics_Time_Feed_RotLen_AlphaU(:,1)*365,Dynamics_W_Feed_RotLen(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'Minimum Length','(B)'},'FontSize', 16)
        legend1=legend([p1 p2 p3],{'Societal Optimum ~~~','Aquaculture Optimum w/o Feed~~~','Aquaculture Optimum w/ Feed'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
        set(legend1,...
    'Position',[0.0938584639410563 0.00827803823612821 0.837143993377685 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

        
        subplot(335)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen_AlphaU(:,2)*365,Dynamics_W_Feed_RotLen(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(E)'},'FontSize', 16)
         
         
        subplot(338)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen_AlphaU(:,3)*365,Dynamics_W_Feed_RotLen(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(H)'},'FontSize', 16)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(333)      
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU2(:,1)*365,Dynamics_W_Feed_FeedSea_AlphaU2(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])  
        xlim([0 300]) 
        title({'Limited Feeding','(C)'},'FontSize', 16)
               
        subplot(336)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU2(:,2)*365,Dynamics_W_Feed_FeedSea_AlphaU2(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(F)'},'FontSize', 16)
        
        
        subplot(339)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU2(:,3)*365,Dynamics_W_Feed_FeedSea_AlphaU2(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(I)'},'FontSize', 16)
        

            han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Time (days)',''}, 'FontSize', 16);
    
    
    saveas(gcf,'SA_DynamicsW_AttackFeed.png'); hold off 
        
        
        
        
        
        
end
    
end  
    


%Holling's Type III exponent
if true

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    N=linspace(1.75,2.25,3);
    
    Dynamics_L_Opt=zeros(Nset+1,length(N));
    Dynamics_P_Opt=zeros(Nset+1,length(N));
    Dynamics_W_Opt=zeros(Nset+1,length(N));
    Dynamics_Time_Opt=zeros(Nset+1,length(N));

    Dynamics_L_NoFeed=zeros(Nset+1,length(N));
    Dynamics_P_NoFeed=zeros(Nset+1,length(N));
    Dynamics_W_NoFeed=zeros(Nset+1,length(N));
    Dynamics_Time_NoFeed=zeros(Nset+1,length(N));
    
    
    Dynamics_L_Feed=zeros(Nset+1,length(N));
    Dynamics_P_Feed=zeros(Nset+1,length(N));
    Dynamics_W_Feed=zeros(Nset+1,length(N));
    Dynamics_U_Feed=zeros(Nset+1,length(N));
    Dynamics_Time_Feed=zeros(Nset+1,length(N));
        

    Dynamics_L_Feed_RotLen2=zeros(Nset+1,length(N));
    Dynamics_P_Feed_RotLen2=zeros(Nset+1,length(N));
    Dynamics_W_Feed_RotLen2=zeros(Nset+1,length(N));
    Dynamics_U_Feed_RotLen2=zeros(Nset+1,length(N));
    Dynamics_Time_Feed_RotLen2=zeros(Nset+1,length(N));
    Dynamics_X_Feed_RotLen2=zeros(Nset+1,length(N));
    Dynamics_I_Feed_RotLen2=zeros(Nset+1,length(N));

 
 
    Dynamics_L_Feed_FeedSea2=zeros(Nset+1,length(N));
    Dynamics_P_Feed_FeedSea2=zeros(Nset+1,length(N));
    Dynamics_W_Feed_FeedSea2=zeros(Nset+1,length(N));
    Dynamics_U_Feed_FeedSea2=zeros(Nset+1,length(N));
    Dynamics_Time_Feed_FeedSea2=zeros(Nset+1,length(N));
    Dynamics_X_Feed_FeedSea2=zeros(Nset+1,length(N));
    Dynamics_I_Feed_FeedSea2=zeros(Nset+1,length(N));

    
        for i=1:length(N)
        
        n=N(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        
        OBJ=3;
        CASE=1;
        GUESS=[];

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_Opt(:,i)=Ls_No;
    Dynamics_P_Opt(:,i)=Ps_No;
    Dynamics_W_Opt(:,i)=Ws_No;
    Dynamics_Time_Opt(:,i)=ts_No;
    
    
        OBJ=1;
        CASE=1;

[ts_NoFeed, Topt_NoFeed, Is_NoFeed, Ws_NoFeed, Xs_NoFeed, Ns_NoFeed, Ls_NoFeed, Ps_NoFeed, Bs_NoFeed, Omegas_NoFeed, Profits_NoFeed, psiWs_NoFeed, psiXs_NoFeed, psiNs_NoFeed, alphaNs_NoFeed, Ratios_NoFeed, Ths_NoFeed, ks_NoFeed, Results_NoFeed] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_NoFeed(:,i)=Ls_NoFeed;
    Dynamics_P_NoFeed(:,i)=Ps_NoFeed;
    Dynamics_W_NoFeed(:,i)=Ws_NoFeed;
    Dynamics_Time_NoFeed(:,i)=ts_NoFeed;

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed(:,i)=Ls;
    Dynamics_P_Feed(:,i)=Ps;
    Dynamics_W_Feed(:,i)=Ws;
    Dynamics_U_Feed(:,i)=Us;
    Dynamics_Time_Feed(:,i)=ts;
    
        end   
        
        
        
        %%% Standardized Rotation Length
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        

        OBJ=1;
        CASE=2;
        
        
        GUESS(CASE).I=Dynamics_I_Feed_RotLen(:,1);
        GUESS(CASE).W=Dynamics_W_Feed_RotLen(:,1);
        GUESS(CASE).X=Dynamics_X_Feed_RotLen(:,1);
        GUESS(CASE).L=Dynamics_L_Feed_RotLen(:,1);
        GUESS(CASE).P=Dynamics_P_Feed_RotLen(:,1);
        GUESS(CASE).U=Dynamics_U_Feed_RotLen(:,1);

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_RotLen2(:,i)=Ls;
    Dynamics_P_Feed_RotLen2(:,i)=Ps;
    Dynamics_W_Feed_RotLen2(:,i)=Ws;
    Dynamics_U_Feed_RotLen2(:,i)=Us;
    Dynamics_Time_Feed_RotLen2(:,i)=ts;
    
        end   
      
        
        
        %%% Limited Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        
        OBJ=1;
        CASE=2;
        
        
            GUESS(CASE).I=Dynamics_I_Feed_FeedSea(:,i);
            GUESS(CASE).W=Dynamics_W_Feed_FeedSea(:,i);
            GUESS(CASE).X=Dynamics_X_Feed_FeedSea(:,i);
            GUESS(CASE).L=Dynamics_L_Feed_FeedSea(:,i);
            GUESS(CASE).P=Dynamics_P_Feed_FeedSea(:,i);
            GUESS(CASE).U=Dynamics_U_Feed_FeedSea(:,i);
        
        

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_FeedSea2(:,i)=Ls;
    Dynamics_P_Feed_FeedSea2(:,i)=Ps;
    Dynamics_W_Feed_FeedSea2(:,i)=Ws;
    Dynamics_U_Feed_FeedSea2(:,i)=Us;
    Dynamics_Time_Feed_FeedSea2(:,i)=ts;
    Dynamics_X_Feed_FeedSea2(:,i)=Xs;
    Dynamics_I_Feed_FeedSea2(:,i)=Is;
        end   

        
        
        end
    

    
        
        
       %%%figure 
if true
    
    
    %Colors
    bC=[0, 0.4470, 0.7410];
    rC=[0.8500, 0.3250, 0.0980];
    yC=[0.9290, 0.6940, 0.1250];
        
       fig= figure
     
        subplot(331)      
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,1)*365,Dynamics_W_Feed(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'No Policy','(A)'},'FontSize', 16)
        ylabel({'Infected Snails','(Specialized Prawns)'}, 'FontSize', 14) %,'interpreter','latex'
        
        subplot(334)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,2)*365,Dynamics_W_Feed(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        ylabel({'Infected Snails','(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(D)'},'FontSize', 16)
        
        subplot(337)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,3)*365,Dynamics_W_Feed(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        ylabel({'Infected Snails','(Generalist Prawns)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(G)'},'FontSize', 16)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(332)      
        p1=plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        p2=plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        p3=plot(Dynamics_Time_Feed_RotLen2(:,1)*365,Dynamics_W_Feed_RotLen(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])  
        xlim([0 300])
        title({'Minimum Length','(B)'},'FontSize', 16)
        legend1=legend([p1 p2 p3],{'Societal Optimum ~~~','Aquaculture Optimum w/o Feed~~~','Aquaculture Optimum w/ Feed'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
        set(legend1,...
    'Position',[0.0938584639410563 0.00827803823612821 0.837143993377685 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

        
        subplot(335)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen2(:,2)*365,Dynamics_W_Feed_RotLen(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(E)'},'FontSize', 16)
         
         
        subplot(338)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen2(:,3)*365,Dynamics_W_Feed_RotLen(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(H)'},'FontSize', 16)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(333)      
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea2(:,1)*365,Dynamics_W_Feed_FeedSea2(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])   
        xlim([0 300])
        title({'Limited Feeding','(C)'},'FontSize', 16)
               
        subplot(336)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea2(:,2)*365,Dynamics_W_Feed_FeedSea2(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(F)'},'FontSize', 16)
        
        
        subplot(339)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea2(:,3)*365,Dynamics_W_Feed_FeedSea2(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(I)'},'FontSize', 16)
        

            han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    %xlabel(han,{'','','','Time (days)'}, 'FontSize', 16);
    xlabel(han,{'Time (days)',''}, 'FontSize', 16);
        
        
   saveas(gcf,'SA_DynamicsW_Holling.png'); hold off 
        
        
        %figure for feed
        
    
        
end
    
        
        
end


    
    
end



end

%For Nset=70
if false
for Nset=70
    
    %Attack rate on supplemental feed
if true

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    AlphaU=linspace(alphaU*0.5,alphaU*1.5,3);
    
    Dynamics_L_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Opt=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Opt_AlphaU=zeros(Nset+1,length(AlphaU));
    
    
    Dynamics_L_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_P_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_W_NoFeed=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_NoFeed_AlphaU=zeros(Nset+1,length(AlphaU));

    
    Dynamics_L_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_AlphaU=zeros(Nset+1,length(AlphaU));
    
    
    Dynamics_L_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed_RotLen=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_RotLen_AlphaU=zeros(Nset+1,length(AlphaU));
    

    Dynamics_L_Feed_FeedSea_AlphaU3=zeros(Nset+1,length(AlphaU));
    Dynamics_P_Feed_FeedSea_AlphaU3=zeros(Nset+1,length(AlphaU));
    Dynamics_W_Feed_FeedSea_AlphaU3=zeros(Nset+1,length(AlphaU));
    Dynamics_U_Feed_FeedSea_AlphaU3=zeros(Nset+1,length(AlphaU));
    Dynamics_Time_Feed_FeedSea_AlphaU3=zeros(Nset+1,length(AlphaU));
    Dynamics_X_Feed_FeedSea_AlphaU3=zeros(Nset+1,length(AlphaU));
    Dynamics_I_Feed_FeedSea_AlphaU3=zeros(Nset+1,length(AlphaU));

    
        for i=1:length(AlphaU)
        
        alphaU=AlphaU(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        GUESS=[];
        OBJ=3;
        CASE=1;
        
[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_Opt(:,i)=Ls_No;
    Dynamics_P_Opt(:,i)=Ps_No;
    Dynamics_W_Opt(:,i)=Ws_No;
    Dynamics_Time_Opt_AlphaU(:,i)=ts_No;
    
    
        OBJ=1;
        CASE=1;

[ts_NoFeed, Topt_NoFeed, Is_NoFeed, Ws_NoFeed, Xs_NoFeed, Ns_NoFeed, Ls_NoFeed, Ps_NoFeed, Bs_NoFeed, Omegas_NoFeed, Profits_NoFeed, psiWs_NoFeed, psiXs_NoFeed, psiNs_NoFeed, alphaNs_NoFeed, Ratios_NoFeed, Ths_NoFeed, ks_NoFeed, Results_NoFeed] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_NoFeed(:,i)=Ls_NoFeed;
    Dynamics_P_NoFeed(:,i)=Ps_NoFeed;
    Dynamics_W_NoFeed(:,i)=Ws_NoFeed;
    Dynamics_Time_NoFeed_AlphaU(:,i)=ts_NoFeed;

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed(:,i)=Ls;
    Dynamics_P_Feed(:,i)=Ps;
    Dynamics_W_Feed(:,i)=Ws;
    Dynamics_U_Feed_AlphaU(:,i)=Us;
    Dynamics_Time_Feed_AlphaU(:,i)=ts;
    
        end   
        
        
        
        %%% Standardized Rotation Length
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        
       
        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_RotLen(:,i)=Ls;
    Dynamics_P_Feed_RotLen(:,i)=Ps;
    Dynamics_W_Feed_RotLen(:,i)=Ws;
    Dynamics_U_Feed_RotLen_AlphaU(:,i)=Us;
    Dynamics_Time_Feed_RotLen_AlphaU(:,i)=ts;
    
        end   
      
        
        
        %%% Limited Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        
      
            GUESS(CASE).I=Dynamics_I_Feed_FeedSea_AlphaU2(:,i);
            GUESS(CASE).W=Dynamics_W_Feed_FeedSea_AlphaU2(:,i);
            GUESS(CASE).X=Dynamics_X_Feed_FeedSea_AlphaU2(:,i);
            GUESS(CASE).L=Dynamics_L_Feed_FeedSea_AlphaU2(:,i);
            GUESS(CASE).P=Dynamics_P_Feed_FeedSea_AlphaU2(:,i);
            GUESS(CASE).U=Dynamics_U_Feed_FeedSea_AlphaU2(:,i);
    

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_FeedSea_AlphaU3(:,i)=Ls;
    Dynamics_P_Feed_FeedSea_AlphaU3(:,i)=Ps;
    Dynamics_W_Feed_FeedSea_AlphaU3(:,i)=Ws;
    Dynamics_U_Feed_FeedSea_AlphaU3(:,i)=Us;
    Dynamics_Time_Feed_FeedSea_AlphaU3(:,i)=ts;
    Dynamics_X_Feed_FeedSea_AlphaU3(:,i)=Xs;
    Dynamics_I_Feed_FeedSea_AlphaU3(:,i)=Is;

        end   

        
        
        end
    
       %%%figure 
if true
    
    
    %Colors
    bC=[0, 0.4470, 0.7410];
    rC=[0.8500, 0.3250, 0.0980];
    yC=[0.9290, 0.6940, 0.1250];
        
       fig= figure
     
        subplot(331)      
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,1)*365,Dynamics_W_Feed(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1]) 
        xlim([0 300]) 
        title({'No Policy','(A)'},'FontSize', 16)
        ylabel({'Infected Snails','(Low Feed Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'
        
        subplot(334)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,2)*365,Dynamics_W_Feed(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        ylabel({'Infected Snails','(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(D)'},'FontSize', 16)
        
        subplot(337)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_AlphaU(:,3)*365,Dynamics_W_Feed(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        ylabel({'Infected Snails','(High Feed Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(G)'},'FontSize', 16)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(332)      
        p1=plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        p2=plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        p3=plot(Dynamics_Time_Feed_RotLen_AlphaU(:,1)*365,Dynamics_W_Feed_RotLen(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'Minimum Length','(B)'},'FontSize', 16)
        legend1=legend([p1 p2 p3],{'Societal Optimum ~~~','Aquaculture Optimum w/o Feed~~~','Aquaculture Optimum w/ Feed'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
        set(legend1,...
    'Position',[0.0938584639410563 0.00827803823612821 0.837143993377685 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

        
        subplot(335)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen_AlphaU(:,2)*365,Dynamics_W_Feed_RotLen(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(E)'},'FontSize', 16)
         
         
        subplot(338)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen_AlphaU(:,3)*365,Dynamics_W_Feed_RotLen(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(H)'},'FontSize', 16)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(333)      
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU3(:,1)*365,Dynamics_W_Feed_FeedSea_AlphaU3(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])  
        xlim([0 300]) 
        title({'Limited Feeding','(C)'},'FontSize', 16)
               
        subplot(336)
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU3(:,2)*365,Dynamics_W_Feed_FeedSea_AlphaU3(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(F)'},'FontSize', 16)
        
        
        subplot(339)   
        plot(Dynamics_Time_Opt_AlphaU(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed_AlphaU(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea_AlphaU3(:,3)*365,Dynamics_W_Feed_FeedSea_AlphaU3(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'(I)'},'FontSize', 16)
        

            han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Time (days)',''}, 'FontSize', 16);
    
    
   % saveas(gcf,'SA_DynamicsW_AttackFeed.png'); hold off 
        
        
        
        
        
        
end
    
end  
    


%Holling's Type III exponent
if true

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    N=linspace(1.75,2.25,3);
    
    Dynamics_L_Opt=zeros(Nset+1,length(N));
    Dynamics_P_Opt=zeros(Nset+1,length(N));
    Dynamics_W_Opt=zeros(Nset+1,length(N));
    Dynamics_Time_Opt=zeros(Nset+1,length(N));

    Dynamics_L_NoFeed=zeros(Nset+1,length(N));
    Dynamics_P_NoFeed=zeros(Nset+1,length(N));
    Dynamics_W_NoFeed=zeros(Nset+1,length(N));
    Dynamics_Time_NoFeed=zeros(Nset+1,length(N));
    
    
    Dynamics_L_Feed=zeros(Nset+1,length(N));
    Dynamics_P_Feed=zeros(Nset+1,length(N));
    Dynamics_W_Feed=zeros(Nset+1,length(N));
    Dynamics_U_Feed=zeros(Nset+1,length(N));
    Dynamics_Time_Feed=zeros(Nset+1,length(N));
        

    Dynamics_L_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_P_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_W_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_U_Feed_RotLen=zeros(Nset+1,length(N));
    Dynamics_Time_Feed_RotLen=zeros(Nset+1,length(N));
    
    Dynamics_L_Feed_FeedSea3=zeros(Nset+1,length(N));
    Dynamics_P_Feed_FeedSea3=zeros(Nset+1,length(N));
    Dynamics_W_Feed_FeedSea3=zeros(Nset+1,length(N));
    Dynamics_U_Feed_FeedSea3=zeros(Nset+1,length(N));
    Dynamics_Time_Feed_FeedSea3=zeros(Nset+1,length(N));
    Dynamics_X_Feed_FeedSea3=zeros(Nset+1,length(N));
    Dynamics_I_Feed_FeedSea3=zeros(Nset+1,length(N));

    
        for i=1:length(N)
        
        n=N(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        
        OBJ=3;
        CASE=1;
        GUESS=[];

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_Opt(:,i)=Ls_No;
    Dynamics_P_Opt(:,i)=Ps_No;
    Dynamics_W_Opt(:,i)=Ws_No;
    Dynamics_Time_Opt(:,i)=ts_No;
    
    
        OBJ=1;
        CASE=1;

[ts_NoFeed, Topt_NoFeed, Is_NoFeed, Ws_NoFeed, Xs_NoFeed, Ns_NoFeed, Ls_NoFeed, Ps_NoFeed, Bs_NoFeed, Omegas_NoFeed, Profits_NoFeed, psiWs_NoFeed, psiXs_NoFeed, psiNs_NoFeed, alphaNs_NoFeed, Ratios_NoFeed, Ths_NoFeed, ks_NoFeed, Results_NoFeed] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    Dynamics_L_NoFeed(:,i)=Ls_NoFeed;
    Dynamics_P_NoFeed(:,i)=Ps_NoFeed;
    Dynamics_W_NoFeed(:,i)=Ws_NoFeed;
    Dynamics_Time_NoFeed(:,i)=ts_NoFeed;

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed(:,i)=Ls;
    Dynamics_P_Feed(:,i)=Ps;
    Dynamics_W_Feed(:,i)=Ws;
    Dynamics_U_Feed(:,i)=Us;
    Dynamics_Time_Feed(:,i)=ts;
    
        end   
        
        
        
        %%% Standardized Rotation Length
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        

        OBJ=1;
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_RotLen(:,i)=Ls;
    Dynamics_P_Feed_RotLen(:,i)=Ps;
    Dynamics_W_Feed_RotLen(:,i)=Ws;
    Dynamics_U_Feed_RotLen(:,i)=Us;
    Dynamics_Time_Feed_RotLen(:,i)=ts;
    
        end   
      
        
        
        %%% Limited Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        
        OBJ=1;
        CASE=2;
        
        
            GUESS(CASE).I=Dynamics_I_Feed_FeedSea2(:,i);
            GUESS(CASE).W=Dynamics_W_Feed_FeedSea2(:,i);
            GUESS(CASE).X=Dynamics_X_Feed_FeedSea2(:,i);
            GUESS(CASE).L=Dynamics_L_Feed_FeedSea2(:,i);
            GUESS(CASE).P=Dynamics_P_Feed_FeedSea2(:,i);
            GUESS(CASE).U=Dynamics_U_Feed_FeedSea2(:,i);
        
        

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_NoFeed,Ws_NoFeed(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    Dynamics_L_Feed_FeedSea3(:,i)=Ls;
    Dynamics_P_Feed_FeedSea3(:,i)=Ps;
    Dynamics_W_Feed_FeedSea3(:,i)=Ws;
    Dynamics_U_Feed_FeedSea3(:,i)=Us;
    Dynamics_Time_Feed_FeedSea3(:,i)=ts;
    Dynamics_X_Feed_FeedSea3(:,i)=Xs;
    Dynamics_I_Feed_FeedSea3(:,i)=Is;
        end   

        
        
        end
    

    
        
        
       %%%figure 
if true
    
    
    %Colors
    bC=[0, 0.4470, 0.7410];
    rC=[0.8500, 0.3250, 0.0980];
    yC=[0.9290, 0.6940, 0.1250];
        
       fig= figure
     
        subplot(331)      
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,1)*365,Dynamics_W_Feed(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        title({'No Policy','(A)'},'FontSize', 16)
        ylabel({'Infected Snails','(Specialized Prawns)'}, 'FontSize', 14) %,'interpreter','latex'
        
        subplot(334)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,2)*365,Dynamics_W_Feed(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300]) 
        ylabel({'Infected Snails','(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(D)'},'FontSize', 16)
        
        subplot(337)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed(:,3)*365,Dynamics_W_Feed(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        ylabel({'Infected Snails','(Generalist Prawns)'}, 'FontSize', 14) %,'interpreter','latex'
        title({'(G)'},'FontSize', 16)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(332)      
        p1=plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        p2=plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        p3=plot(Dynamics_Time_Feed_RotLen(:,1)*365,Dynamics_W_Feed_RotLen(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])  
        xlim([0 300])
        title({'Minimum Length','(B)'},'FontSize', 16)
        legend1=legend([p1 p2 p3],{'Societal Optimum ~~~','Aquaculture Optimum w/o Feed~~~','Aquaculture Optimum w/ Feed'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
        set(legend1,...
    'Position',[0.0938584639410563 0.00827803823612821 0.837143993377685 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

        
        subplot(335)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen(:,2)*365,Dynamics_W_Feed_RotLen(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(E)'},'FontSize', 16)
         
         
        subplot(338)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_RotLen(:,3)*365,Dynamics_W_Feed_RotLen(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(H)'},'FontSize', 16)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(333)      
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea3(:,1)*365,Dynamics_W_Feed_FeedSea3(:,1)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])   
        xlim([0 300])
        title({'Limited Feeding','(C)'},'FontSize', 16)
               
        subplot(336)
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea3(:,2)*365,Dynamics_W_Feed_FeedSea3(:,2)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(F)'},'FontSize', 16)
        
        
        subplot(339)   
        plot(Dynamics_Time_Opt(:,2)*365,Dynamics_W_Opt(:,2)./x0ic(2),'LineWidth',3,'color',bC,'LineStyle','-'); hold on
        plot(Dynamics_Time_NoFeed(:,2)*365,Dynamics_W_NoFeed(:,2)./x0ic(2),'LineWidth',3,'color',rC,'LineStyle','-'); hold on
        plot(Dynamics_Time_Feed_FeedSea3(:,3)*365,Dynamics_W_Feed_FeedSea3(:,3)./x0ic(2),'LineWidth',3,'color',yC,'LineStyle','-'); hold on
        ylim([0 1])
        xlim([0 300])
        title({'(I)'},'FontSize', 16)
        

            han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    %xlabel(han,{'','','','Time (days)'}, 'FontSize', 16);
    xlabel(han,{'Time (days)',''}, 'FontSize', 16);
        
        
   %saveas(gcf,'SA_DynamicsW_Holling.png'); hold off 
        
        
        %figure for feed
        
    
        
end
    
        
        
end


    
    
end



end



% Figure for Dynamics of Feeding
if true
                 
            bC=[0, 0.4470, 0.7410];
            rC=[0.8500, 0.3250, 0.0980];
            yC=[0.9290, 0.6940, 0.1250];
            yC2=[0.0290, 0.6940, 0.1250];
            yC3=[0.9290, 0.3940, 0.8250];
            
         fig=figure
        subplot(3,2,1)
        plot(365.*Dynamics_Time_Feed_AlphaU(:,1),Dynamics_U_Feed_AlphaU(:,1),'LineWidth',3,'color',bC); hold on %,'LineStyle','--'
        plot(365.*Dynamics_Time_Feed_RotLen_AlphaU(:,1),Dynamics_U_Feed_RotLen_AlphaU(:,1),'LineWidth',3,'color',rC); hold on
        plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU(:,1),Dynamics_U_Feed_FeedSea_AlphaU(:,1),'LineWidth',3,'color',yC); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU2(:,1),Dynamics_U_Feed_FeedSea_AlphaU2(:,1),'LineWidth',3,'color',yC2); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU3(:,1),Dynamics_U_Feed_FeedSea_AlphaU3(:,1),'LineWidth',3,'color',yC3); hold on
        title({'Feed Conversion Efficiency','(A)'},'FontSize', 16)
        ylabel({'(Low Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'
        ylim([0 3])  
        %xlim([0 1])
        
        
        subplot(3,2,3)
        plot(365.*Dynamics_Time_Feed_AlphaU(:,2),Dynamics_U_Feed_AlphaU(:,2),'LineWidth',3,'color',bC); hold on
        plot(365.*Dynamics_Time_Feed_RotLen_AlphaU(:,2),Dynamics_U_Feed_RotLen_AlphaU(:,2),'LineWidth',3,'color',rC); hold on
        plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU(:,2),Dynamics_U_Feed_FeedSea_AlphaU(:,2),'LineWidth',3,'color',yC); hold on
       % plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU2(:,2),Dynamics_U_Feed_FeedSea_AlphaU2(:,2),'LineWidth',3,'color',yC2); hold on
       % plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU3(:,2),Dynamics_U_Feed_FeedSea_AlphaU3(:,2),'LineWidth',3,'color',yC3); hold on
        ylim([0 3])  
        %xlim([0 1])
        title({'(C)'},'FontSize', 16)
        ylabel({'(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        
            
        subplot(3,2,5)
        plot(365.*Dynamics_Time_Feed_AlphaU(:,3),Dynamics_U_Feed_AlphaU(:,3),'LineWidth',3,'color',bC); hold on
        plot(365.*Dynamics_Time_Feed_RotLen_AlphaU(:,3),Dynamics_U_Feed_RotLen_AlphaU(:,3),'LineWidth',3,'color',rC); hold on
        plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU(:,3),Dynamics_U_Feed_FeedSea_AlphaU(:,3),'LineWidth',3,'color',yC); hold on
       % plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU2(:,3),Dynamics_U_Feed_FeedSea_AlphaU2(:,3),'LineWidth',3,'color',yC2); hold on
       % plot(365.*Dynamics_Time_Feed_FeedSea_AlphaU3(:,3),Dynamics_U_Feed_FeedSea_AlphaU3(:,3),'LineWidth',3,'color',yC3); hold on
        ylabel({'(High Efficiency)'}, 'FontSize', 14) %,'interpreter','latex'      
        ylim([0 3])  
        %xlim([0 1])
        title({'(E)'},'FontSize', 16)
        
        subplot(3,2,2)
        plot(365.*Dynamics_Time_Feed(:,1),Dynamics_U_Feed(:,1),'LineWidth',3,'color',bC); hold on
        plot(365.*Dynamics_Time_Feed_RotLen(:,1),Dynamics_U_Feed_RotLen(:,1),'LineWidth',3,'color',rC); hold on
        plot(365.*Dynamics_Time_Feed_FeedSea(:,1),Dynamics_U_Feed_FeedSea(:,1),'LineWidth',3,'color',yC); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea2(:,1),Dynamics_U_Feed_FeedSea2(:,1),'LineWidth',3,'color',yC2); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea3(:,1),Dynamics_U_Feed_FeedSea2(:,1),'LineWidth',3,'color',yC3); hold on
        title({'Holling''s Type III Exponent','(B)'},'FontSize', 16)
        ylabel({'(Specialized)'}, 'FontSize', 14) %,'interpreter','latex'
        ylim([0 3])  
        %xlim([0 1])
        
        
        subplot(3,2,4)
        plot(365.*Dynamics_Time_Feed(:,2),Dynamics_U_Feed(:,2),'LineWidth',3,'color',bC); hold on
        plot(365.*Dynamics_Time_Feed_RotLen(:,2),Dynamics_U_Feed_RotLen(:,2),'LineWidth',3,'color',rC); hold on
        plot(365.*Dynamics_Time_Feed_FeedSea(:,2),Dynamics_U_Feed_FeedSea(:,2),'LineWidth',3,'color',yC); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea2(:,2),Dynamics_U_Feed_FeedSea2(:,2),'LineWidth',3,'color',yC2); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea3(:,2),Dynamics_U_Feed_FeedSea3(:,2),'LineWidth',3,'color',yC3); hold on
        ylabel({'(Base Case)'}, 'FontSize', 14) %,'interpreter','latex'
        ylim([0 3])  
        %xlim([0 1])
        title({'(D)'},'FontSize', 16)
            
        subplot(3,2,6)
        p1=plot(365.*Dynamics_Time_Feed(:,3),Dynamics_U_Feed(:,3),'LineWidth',3,'color',bC); hold on
        p2=plot(365.*Dynamics_Time_Feed_RotLen(:,3),Dynamics_U_Feed_RotLen(:,3),'LineWidth',3,'color',rC); hold on
        p3=plot(365.*Dynamics_Time_Feed_FeedSea(:,3),Dynamics_U_Feed_FeedSea(:,3),'LineWidth',3,'color',yC); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea2(:,3),Dynamics_U_Feed_FeedSea2(:,3),'LineWidth',3,'color',yC2); hold on
        %plot(365.*Dynamics_Time_Feed_FeedSea3(:,3),Dynamics_U_Feed_FeedSea3(:,3),'LineWidth',3,'color',yC3); hold on
        ylabel({'(Generalist)'}, 'FontSize', 14) %,'interpreter','latex'
        ylim([0 3])  
        title({'(F)'},'FontSize', 16)
        legend1=legend([p1 p2 p3],{'No Policy ~~~','Minimum Length~~~','Limited Feeding'},'Interpreter','latex','Orientation','horizontal','Location','northeast');      
            set(legend1,...
    'Position',[0.244573208720414 0.00583103707097239 0.562054184504917 0.0352380951245626],...
    'Orientation','horizontal',...
    'Interpreter','latex');




            han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
   ylabel(han,{'Quantity of Feed',''}, 'FontSize', 16);
    xlabel(han,{'Time (days)',''}, 'FontSize', 16);
    
  saveas(gcf,'SA_DynamicsU.png'); hold off 
    
    
    
end


end













%%%%%%%%%%%%%%%%%%%%%%%%
%%% One set of results
%%%%%%%%%%%%%%%%%%%%%%%%
if false
Nset=60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Single-Rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
%Profit

OBJ=1;  
CASE=1;
GUESS=[];

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=1;
CASE=2;
        if true
            GUESS(CASE).I=Is_No;
            GUESS(CASE).W=Ws_No;
            GUESS(CASE).X=Xs_No;
            GUESS(CASE).L=Ls_No;
            GUESS(CASE).P=Ps_No;
            GUESS(CASE).U=0;
        end


[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


%Health

OBJ=5;  
CASE=1;
GUESS=[];

[ts_Health, Topt_Health, Is_Health, Ws_Health, Xs_Health, Ns_Health, Ls_Health, Ps_Health, Bs_Health, Omegas_Health, Profits_Health, psiWs_Health, psiXs_Health, psiNs_Health, alphaNs_Health, Ratios_Health, Ths_Health, ks_Health, Results_Health] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=5;
CASE=2;

        if true
            GUESS(CASE).I=Is_Health;
            GUESS(CASE).W=Ws_Health;
            GUESS(CASE).X=Xs_Health;
            GUESS(CASE).L=Ls_Health;
            GUESS(CASE).P=Ps_Health;
            GUESS(CASE).U=Us;
        end
        

[ts_HealthFeed, Topt_HealthFeed, Us_HealthFeed, Costs_HealthFeed, Is_HealthFeed, Ws_HealthFeed, Xs_HealthFeed, Ns_HealthFeed, Ls_HealthFeed, Ps_HealthFeed, Bs_HealthFeed, Omegas_HealthFeed, Profits_HealthFeed, psiWs_HealthFeed, psiXs_HealthFeed, psiNs_HealthFeed, alphaNs_HealthFeed, Ratios_HealthFeed, Ths_HealthFeed, ks_HealthFeed, ks_I_HealthFeed, ks_E_HealthFeed, Results_HealthFeed] = ...
    SchistoAquaculture_Feed(Topt_Health,Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);



    
%Optimal

OBJ=3;  
CASE=1;
GUESS=[];

[ts_Opt, Topt_Opt, Is_Opt, Ws_Opt, Xs_Opt, Ns_Opt, Ls_Opt, Ps_Opt, Bs_Opt, Omegas_Opt, Profits_Opt, psiWs_Opt, psiXs_Opt, psiNs_Opt, alphaNs_Opt, Ratios_Opt, Ths_Opt, ks_Opt, Results_Opt] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=3;
CASE=2;

        if true
            GUESS(CASE).I=Is_HealthFeed;
            GUESS(CASE).W=Ws_HealthFeed;
            GUESS(CASE).X=Xs_HealthFeed;
            GUESS(CASE).L=Ls_HealthFeed;
            GUESS(CASE).P=Ps_HealthFeed;
            GUESS(CASE).U=Us_HealthFeed;
        end
        

[ts_OptFeed, Topt_OptFeed, Us_OptFeed, Costs_OptFeed, Is_OptFeed, Ws_OptFeed, Xs_OptFeed, Ns_OptFeed, Ls_OptFeed, Ps_OptFeed, Bs_OptFeed, Omegas_OptFeed, Profits_OptFeed, psiWs_OptFeed, psiXs_OptFeed, psiNs_OptFeed, alphaNs_OptFeed, Ratios_OptFeed, Ths_OptFeed, ks_OptFeed, ks_I_OptFeed, ks_E_OptFeed, Results_OptFeed] = ...
    SchistoAquaculture_Feed(Topt_Opt,Ws_Opt(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Infinite Horizon
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true


%Profit

    
OBJ=2;  
CASE=1;
GUESS=[];

[ts_No_Inf, Topt_No_Inf, Is_No_Inf, Ws_No_Inf, Xs_No_Inf, Ns_No_Inf, Ls_No_Inf, Ps_No_Inf, Bs_No_Inf, Omegas_No_Inf, Profits_No_Inf, psiWs_No_Inf, psiXs_No_Inf, psiNs_No_Inf, alphaNs_No_Inf, Ratios_No_Inf, Ths_No_Inf, ks_No_Inf, Results_No_Inf] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

OBJ=2;
CASE=2;
        if true
            GUESS(CASE).I=Is_No_Inf;
            GUESS(CASE).W=Ws_No_Inf;
            GUESS(CASE).X=Xs_No_Inf;
            GUESS(CASE).L=Ls_No_Inf;
            GUESS(CASE).P=Ps_No_Inf;
            GUESS(CASE).U=0;
        end


[ts_Inf, Topt_Inf, Us_Inf, Costs_Inf, Is_Inf, Ws_Inf, Xs_Inf, Ns_Inf, Ls_Inf, Ps_Inf, Bs_Inf, Omegas_Inf, Profits_Inf, psiWs_Inf, psiXs_Inf, psiNs_Inf, alphaNs_Inf, Ratios_Inf, Ths_Inf, ks_Inf, ks_I_Inf, ks_E_Inf, Results_Inf] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);



%Health

OBJ=6;  
CASE=1;
GUESS=[];

[ts_Health_Inf, Topt_Health_Inf, Is_Health_Inf, Ws_Health_Inf, Xs_Health_Inf, Ns_Health_Inf, Ls_Health_Inf, Ps_Health_Inf, Bs_Health_Inf, Omegas_Health_Inf, Profits_Health_Inf, psiWs_Health_Inf, psiXs_Health_Inf, psiNs_Health_Inf, alphaNs_Health_Inf, Ratios_Health_Inf, Ths_Health_Inf, ks_Health_Inf, Results_Health_Inf] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=6;
CASE=2;

        if true
            GUESS(CASE).I=Is_Health_Inf;
            GUESS(CASE).W=Ws_Health_Inf;
            GUESS(CASE).X=Xs_Health_Inf;
            GUESS(CASE).L=Ls_Health_Inf;
            GUESS(CASE).P=Ps_Health_Inf;
            GUESS(CASE).U=Us_Inf;
        end
        

[ts_HealthFeed_Inf, Topt_HealthFeed_Inf, Us_HealthFeed_Inf, Costs_HealthFeed_Inf, Is_HealthFeed_Inf, Ws_HealthFeed_Inf, Xs_HealthFeed_Inf, Ns_HealthFeed_Inf, Ls_HealthFeed_Inf, Ps_HealthFeed_Inf, Bs_HealthFeed_Inf, Omegas_HealthFeed_Inf, Profits_HealthFeed_Inf, psiWs_HealthFeed_Inf, psiXs_HealthFeed_Inf, psiNs_HealthFeed_Inf, alphaNs_HealthFeed_Inf, Ratios_HealthFeed_Inf, Ths_HealthFeed_Inf, ks_HealthFeed_Inf, ks_I_HealthFeed_Inf, ks_E_HealthFeed_Inf, Results_HealthFeed_Inf] = ...
    SchistoAquaculture_Feed(Topt_Health_Inf,Ws_Health_Inf(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);



%Optimal
OBJ=4;  
CASE=1;
GUESS=[];

[ts_Opt_Inf, Topt_Opt_Inf, Is_Opt_Inf, Ws_Opt_Inf, Xs_Opt_Inf, Ns_Opt_Inf, Ls_Opt_Inf, Ps_Opt_Inf, Bs_Opt_Inf, Omegas_Opt_Inf, Profits_Opt_Inf, psiWs_Opt_Inf, psiXs_Opt_Inf, psiNs_Opt_Inf, alphaNs_Opt_Inf, Ratios_Opt_Inf, Ths_Opt_Inf, ks_Opt_Inf, Results_Opt_Inf] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=4;
CASE=2;

        if true
            GUESS(CASE).I=Is_Inf;
            GUESS(CASE).W=Ws_Inf;
            GUESS(CASE).X=Xs_Inf;
            GUESS(CASE).L=Ls_Inf;
            GUESS(CASE).P=Ps_Inf;
            GUESS(CASE).U=Us_Inf;
        end

[ts_OptFeed_Inf, Topt_OptFeed_Inf, Us_OptFeed_Inf, Costs_OptFeed_Inf, Is_OptFeed_Inf, Ws_OptFeed_Inf, Xs_OptFeed_Inf, Ns_OptFeed_Inf, Ls_OptFeed_Inf, Ps_OptFeed_Inf, Bs_OptFeed_Inf, Omegas_OptFeed_Inf, Profits_OptFeed_Inf, psiWs_OptFeed_Inf, psiXs_OptFeed_Inf, psiNs_OptFeed_Inf, alphaNs_OptFeed_Inf, Ratios_OptFeed_Inf, Ths_OptFeed_Inf, ks_OptFeed_Inf, ks_I_OptFeed_Inf, ks_E_OptFeed_Inf, Results_OptFeed_Inf] = ...
    SchistoAquaculture_Feed(Topt_Opt_Inf,Ws_Opt_Inf(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Second Run 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Single-Rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
%Profit

OBJ=1;  
CASE=1;
        if true
            GUESS(CASE).I=Is_No;
            GUESS(CASE).W=Ws_No;
            GUESS(CASE).X=Xs_No;
            GUESS(CASE).L=Ls_No;
            GUESS(CASE).P=Ps_No;
        end


[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=1;
CASE=2;
        if true
            GUESS(CASE).I=Is;
            GUESS(CASE).W=Ws;
            GUESS(CASE).X=Xs;
            GUESS(CASE).L=Ls;
            GUESS(CASE).P=Ps;
            GUESS(CASE).U=Us;
        end


[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


%Health

OBJ=5;  
CASE=1;
        if true
            GUESS(CASE).I=Is_Health;
            GUESS(CASE).W=Ws_Health;
            GUESS(CASE).X=Xs_Health;
            GUESS(CASE).L=Ls_Health;
            GUESS(CASE).P=Ps_Health;
        end
        

[ts_Health, Topt_Health, Is_Health, Ws_Health, Xs_Health, Ns_Health, Ls_Health, Ps_Health, Bs_Health, Omegas_Health, Profits_Health, psiWs_Health, psiXs_Health, psiNs_Health, alphaNs_Health, Ratios_Health, Ths_Health, ks_Health, Results_Health] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=5;
CASE=2;

        if true
            GUESS(CASE).I=Is_HealthFeed;
            GUESS(CASE).W=Ws_HealthFeed;
            GUESS(CASE).X=Xs_HealthFeed;
            GUESS(CASE).L=Ls_HealthFeed;
            GUESS(CASE).P=Ps_HealthFeed;
            GUESS(CASE).U=Us_HealthFeed;
        end
        

[ts_HealthFeed, Topt_HealthFeed, Us_HealthFeed, Costs_HealthFeed, Is_HealthFeed, Ws_HealthFeed, Xs_HealthFeed, Ns_HealthFeed, Ls_HealthFeed, Ps_HealthFeed, Bs_HealthFeed, Omegas_HealthFeed, Profits_HealthFeed, psiWs_HealthFeed, psiXs_HealthFeed, psiNsHealthFeed, alphaNs_HealthFeed, Ratios_HealthFeed, Ths_HealthFeed, ks_HealthFeed, ks_I_HealthFeed, ks_E_HealthFeed, Results_HealthFeed] = ...
    SchistoAquaculture_Feed(Topt_Health,Ws_Health(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);



    
%Optimal

OBJ=3;  
CASE=1;
        if true
            GUESS(CASE).I=Is_Opt;
            GUESS(CASE).W=Ws_Opt;
            GUESS(CASE).X=Xs_Opt;
            GUESS(CASE).L=Ls_Opt;
            GUESS(CASE).P=Ps_Opt;
        end
        
[ts_Opt, Topt_Opt, Is_Opt, Ws_Opt, Xs_Opt, Ns_Opt, Ls_Opt, Ps_Opt, Bs_Opt, Omegas_Opt, Profits_Opt, psiWs_Opt, psiXs_Opt, psiNs_Opt, alphaNs_Opt, Ratios_Opt, Ths_Opt, ks_Opt, Results_Opt] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=3;
CASE=2;

        if true
            GUESS(CASE).I=Is_OptFeed;
            GUESS(CASE).W=Ws_OptFeed;
            GUESS(CASE).X=Xs_OptFeed;
            GUESS(CASE).L=Ls_OptFeed;
            GUESS(CASE).P=Ps_OptFeed;
            GUESS(CASE).U=Us_OptFeed;
        end
        

[ts_OptFeed, Topt_OptFeed, Us_OptFeed, Costs_OptFeed, Is_OptFeed, Ws_OptFeed, Xs_OptFeed, Ns_OptFeed, Ls_OptFeed, Ps_OptFeed, Bs_OptFeed, Omegas_OptFeed, Profits_OptFeed, psiWs_OptFeed, psiXs_OptFeed, psiNs_OptFeed, alphaNs_OptFeed, Ratios_OptFeed, Ths_OptFeed, ks_OptFeed, ks_I_OptFeed, ks_E_OptFeed, Results_OptFeed] = ...
    SchistoAquaculture_Feed(Topt_Opt,Ws_Opt(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Infinite Horizon
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true


%Profit

    
OBJ=2;  
CASE=1;
        if true
            GUESS(CASE).I=Is_No_Inf;
            GUESS(CASE).W=Ws_No_Inf;
            GUESS(CASE).X=Xs_No_Inf;
            GUESS(CASE).L=Ls_No_Inf;
            GUESS(CASE).P=Ps_No_Inf;
        end

[ts_No_Inf, Topt_No_Inf, Is_No_Inf, Ws_No_Inf, Xs_No_Inf, Ns_No_Inf, Ls_No_Inf, Ps_No_Inf, Bs_No_Inf, Omegas_No_Inf, Profits_No_Inf, psiWs_No_Inf, psiXs_No_Inf, psiNs_No_Inf, alphaNs_No_Inf, Ratios_No_Inf, Ths_No_Inf, ks_No_Inf, Results_No_Inf] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

OBJ=2;
CASE=2;
        if true
            GUESS(CASE).I=Is_Inf;
            GUESS(CASE).W=Ws_Inf;
            GUESS(CASE).X=Xs_Inf;
            GUESS(CASE).L=Ls_Inf;
            GUESS(CASE).P=Ps_Inf;
            GUESS(CASE).U=Us_Inf;
        end


[ts_Inf, Topt_Inf, Us_Inf, Costs_Inf, Is_Inf, Ws_Inf, Xs_Inf, Ns_Inf, Ls_Inf, Ps_Inf, Bs_Inf, Omegas_Inf, Profits_Inf, psiWs_Inf, psiXs_Inf, psiNs_Inf, alphaNs_Inf, Ratios_Inf, Ths_Inf, ks_Inf, ks_I_Inf, ks_E_Inf, Results_Inf] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);



%Health

OBJ=6;  
CASE=1;
        if true
            GUESS(CASE).I=Is_Health_Inf;
            GUESS(CASE).W=Ws_Health_Inf;
            GUESS(CASE).X=Xs_Health_Inf;
            GUESS(CASE).L=Ls_Health_Inf;
            GUESS(CASE).P=Ps_Health_Inf;
        end

[ts_Health_Inf, Topt_Health_Inf, Is_Health_Inf, Ws_Health_Inf, Xs_Health_Inf, Ns_Health_Inf, Ls_Health_Inf, Ps_Health_Inf, Bs_Health_Inf, Omegas_Health_Inf, Profits_Health_Inf, psiWs_Health_Inf, psiXs_Health_Inf, psiNs_Health_Inf, alphaNs_Health_Inf, Ratios_Health_Inf, Ths_Health_Inf, ks_Health_Inf, Results_Health_Inf] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=6;
CASE=2;

        if true
            GUESS(CASE).I=Is_HealthFeed_Inf;
            GUESS(CASE).W=Ws_HealthFeed_Inf;
            GUESS(CASE).X=Xs_HealthFeed_Inf;
            GUESS(CASE).L=Ls_HealthFeed_Inf;
            GUESS(CASE).P=Ps_HealthFeed_Inf;
            GUESS(CASE).U=Us_HealthFeed_Inf;
        end
        

[ts_HealthFeed_Inf, Topt_HealthFeed_Inf, Us_HealthFeed_Inf, Costs_HealthFeed_Inf, Is_HealthFeed_Inf, Ws_HealthFeed_Inf, Xs_HealthFeed_Inf, Ns_HealthFeed_Inf, Ls_HealthFeed_Inf, Ps_HealthFeed_Inf, Bs_HealthFeed_Inf, Omegas_HealthFeed_Inf, Profits_HealthFeed_Inf, psiWs_HealthFeed_Inf, psiXs_HealthFeed_Inf, psiNsHealthFeed_Inf, alphaNs_HealthFeed_Inf, Ratios_HealthFeed_Inf, Ths_HealthFeed_Inf, ks_HealthFeed_Inf, ks_I_HealthFeed_Inf, ks_E_HealthFeed_Inf, Results_HealthFeed_Inf] = ...
    SchistoAquaculture_Feed(Topt_Health_Inf,Ws_Health_Inf(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);



%Optimal
OBJ=4;  
CASE=1;
        if true
            GUESS(CASE).I=Is_Opt_Inf;
            GUESS(CASE).W=Ws_Opt_Inf;
            GUESS(CASE).X=Xs_Opt_Inf;
            GUESS(CASE).L=Ls_Opt_Inf;
            GUESS(CASE).P=Ps_Opt_Inf;
        end


[ts_Opt_Inf, Topt_Opt_Inf, Is_Opt_Inf, Ws_Opt_Inf, Xs_Opt_Inf, Ns_Opt_Inf, Ls_Opt_Inf, Ps_Opt_Inf, Bs_Opt_Inf, Omegas_Opt_Inf, Profits_Opt_Inf, psiWs_Opt_Inf, psiXs_Opt_Inf, psiNs_Opt_Inf, alphaNs_Opt_Inf, Ratios_Opt_Inf, Ths_Opt_Inf, ks_Opt_Inf, Results_Opt_Inf] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


OBJ=4;
CASE=2;

        if true
            GUESS(CASE).I=Is_OptFeed_Inf;
            GUESS(CASE).W=Ws_OptFeed_Inf;
            GUESS(CASE).X=Xs_OptFeed_Inf;
            GUESS(CASE).L=Ls_OptFeed_Inf;
            GUESS(CASE).P=Ps_OptFeed_Inf;
            GUESS(CASE).U=Us_OptFeed_Inf;
        end

[ts_OptFeed_Inf, Topt_OptFeed_Inf, Us_OptFeed_Inf, Costs_OptFeed_Inf, Is_OptFeed_Inf, Ws_OptFeed_Inf, Xs_OptFeed_Inf, Ns_OptFeed_Inf, Ls_OptFeed_Inf, Ps_OptFeed_Inf, Bs_OptFeed_Inf, Omegas_OptFeed_Inf, Profits_OptFeed_Inf, psiWs_OptFeed_Inf, psiXs_OptFeed_Inf, psiNs_OptFeed_Inf, alphaNs_OptFeed_Inf, Ratios_OptFeed_Inf, Ths_OptFeed_Inf, ks_OptFeed_Inf, ks_I_OptFeed_Inf, ks_E_OptFeed_Inf, Results_OptFeed_Inf] = ...
    SchistoAquaculture_Feed(Topt_Opt_Inf,Ws_Opt_Inf(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

end
end



%%% Figure Single-Rotation Harvest & Feed
if true
    
fig=figure
subplot(231)
plot(ts_Opt*365,Profits_Opt,'LineWidth',3); hold on
plot(ts_OptFeed*365,Profits_OptFeed,'LineWidth',3); hold on
%xlim([0 400])
%plot([0 400],[0 0],'LineWidth',1, 'color','k')
%
subplot(232)
plot(ts_Health.*365,Profits_Health,'LineWidth',3);hold on
plot(ts_HealthFeed.*365,Profits_HealthFeed,'LineWidth',3);
%xlim([0 400])
%plot([0 400],[0 0],'LineWidth',1, 'color','k')
%
subplot(233)
plot(ts_No*365,Profits_No,'LineWidth',3); hold on
plot(ts*365,Profits,'LineWidth',3);
%xlim([0 400])
%plot([0 400],[0 0],'LineWidth',1, 'color','k')
%%%%%%%%%%%%%%%%%%%%%
%%% Infected Snails
%%%%%%%%%%%%%%%%%%%%%
subplot(234)
plot(ts_Opt*365,Ws_Opt./Ws_Opt(1),'LineWidth',3);hold on
plot(ts_OptFeed*365,Ws_OptFeed./Ws_OptFeed(1),'LineWidth',3);hold on
%xlim([0 200])
%
subplot(235)
plot(ts_Health.*365,Ws_Health./Ws_Health(1),'LineWidth',3);hold on
plot(ts_HealthFeed.*365,Ws_HealthFeed./Ws_HealthFeed(1),'LineWidth',3);
%xlim([0 200])
%
subplot(236)
plot(ts_No*365,Ws_No./Ws_No(1),'LineWidth',3); hold on
plot(ts*365,Ws./Ws(1),'LineWidth',3);
%xlim([0 200])

 saveas(gcf,'SingleRotation.png'); hold off 

%%%%%%%%%%%%%%%%%%%%%
%%% Feeding Paths
%%%%%%%%%%%%%%%%%%%%%
figure
subplot(311)
plot(ts_Opt*365,Us_OptFeed,'LineWidth',3);hold on
%xlim([0 400])
%
subplot(312)
plot(ts_HealthFeed.*365,Us_HealthFeed,'LineWidth',3);
%xlim([0 400])

%
subplot(313)
plot(ts*365,Us,'LineWidth',3);
%xlim([0 400])

 saveas(gcf,'SingleRotation_Feed.png'); hold off 
end


%%%%% Figure Infinite-horizon Harvest & Feed
if true
    
fig=figure
subplot(231)
plot(ts_Opt_Inf*365,Profits_Opt_Inf,'LineWidth',3); hold on
plot(ts_OptFeed_Inf*365,Profits_OptFeed_Inf,'LineWidth',3); hold on
%xlim([0 200])
%plot([0 200],[0 0],'LineWidth',1, 'color','k')
%
subplot(232)
plot(ts_Health_Inf.*365,Profits_Health_Inf,'LineWidth',3);hold on
plot(ts_HealthFeed_Inf.*365,Profits_HealthFeed_Inf,'LineWidth',3);
%xlim([0 1200])
%plot([0 1200],[0 0],'LineWidth',1, 'color','k')
%
subplot(233)
plot(ts_No_Inf*365,Profits_No_Inf,'LineWidth',3); hold on
plot(ts_Inf*365,Profits_Inf,'LineWidth',3);
%xlim([0 200])
%plot([0 200],[0 0],'LineWidth',1, 'color','k')
%%%%%%%%%%%%%%%%%%%%%
%%% Infected Snails
%%%%%%%%%%%%%%%%%%%%%
subplot(234)
plot(ts_Opt_Inf*365,Ws_Opt_Inf./Ws_Opt_Inf(1),'LineWidth',3);hold on
plot(ts_OptFeed_Inf*365,Ws_OptFeed_Inf./Ws_OptFeed_Inf(1),'LineWidth',3);hold on
%xlim([0 200])
%
subplot(235)
plot(ts_Health_Inf.*365,Ws_Health_Inf./Ws_Health_Inf(1),'LineWidth',3);hold on
plot(ts_HealthFeed_Inf.*365,Ws_HealthFeed_Inf./Ws_HealthFeed_Inf(1),'LineWidth',3);
%xlim([0 1200])
%
subplot(236)
plot(ts_No_Inf*365,Ws_No_Inf./Ws_No_Inf(1),'LineWidth',3); hold on
plot(ts_Inf*365,Ws_Inf./Ws_Inf(1),'LineWidth',3);
%xlim([0 200])

 saveas(gcf,'InfiniteHorizon.png'); hold off 
 
%%%%%%%%%%%%%%%%%%%%%
%%% Feeding Paths
%%%%%%%%%%%%%%%%%%%%%
figure
subplot(311)
plot(ts_OptFeed_Inf*365,Us_OptFeed_Inf,'LineWidth',3);hold on
%xlim([0 400])
%
subplot(312)
plot(ts_HealthFeed_Inf.*365,Us_HealthFeed_Inf,'LineWidth',3);
%xlim([0 400])

%
subplot(313)
plot(ts_Inf*365,Us_Inf,'LineWidth',3);
%xlim([0 400])
 saveas(gcf,'InfiniteHorizon_Feed.png'); hold off 
end


%Information for new Table
if true
H=1e6;

% % % Single Rotation
% % No feed

%Health
round(ts_Health(end)*365)
round(H*(Is_Health(1)-Is_Health(end)))
round(Profits_Health(end))
round(psiNs_Health(end))

%Profits
round(ts_No(end)*365)
round(H*(Is_No(1)-Is_No(end)))
round(Profits_No(end))
round(psiNs_No(end))


%Social
round(ts_Opt(end)*365)
round(H*(Is_Opt(1)-Is_Opt(end)))
round(Profits_Opt(end))
round(psiNs_Opt(end))

% % Feed

%Health
round(ts_HealthFeed(end)*365)
round(H*(Is_HealthFeed(1)-Is_HealthFeed(end)))
round(Profits_HealthFeed(end))
round(psiNs_HealthFeed(end))


%Profits
round(ts(end)*365)
round(H*(Is(1)-Is(end)))
round(Profits(end))
round(psiNs(end))


%Social
round(ts_OptFeed(end)*365)
round(H*(Is_OptFeed(1)-Is_OptFeed(end)))
round(Profits_OptFeed(end))
round(psiNs_OptFeed(end))




% % % % % % % % % % % % % % % % % 
% % % Infinite Horizon % % % % % %
% % % % % % % % % % % % % % % % %

% % No feed

%Health
round(ts_Health_Inf(end)*365)
round(H*(Is_Health_Inf(1)-Is_Health_Inf(end)))
round(Profits_Health_Inf(end))
round(psiNs_Health_Inf(end))

%Profits
round(ts_No_Inf(end)*365)
round(H*(Is_No_Inf(1)-Is_No_Inf(end)))
round(Profits_No_Inf(end))
round(psiNs_No_Inf(end))

%Social
round(ts_Opt_Inf(end)*365)
round(H*(Is_Opt_Inf(1)-Is_Opt_Inf(end)))
round(Profits_Opt_Inf(end))
round(psiNs_Opt_Inf(end))


% % Feed

%Health
round(ts_HealthFeed_Inf(end)*365)
round(H*(Is_HealthFeed_Inf(1)-Is_HealthFeed_Inf(end)))
round(Profits_HealthFeed_Inf(end))
round(psiNs_HealthFeed_Inf(end))

%Profits
round(ts_Inf(end)*365)
round(H*(Is_Inf(1)-Is_Inf(end)))
round(Profits_Inf(end))
round(psiNs_Inf(end))

%Social
round(ts_OptFeed_Inf(end)*365)
round(H*(Is_OptFeed_Inf(1)-Is_OptFeed_Inf(end)))
round(Profits_OptFeed_Inf(end))
round(psiNs_OptFeed_Inf(end))


end




%Information for Table

if false
    
    
    
H=1e6;

%Number of years until convergence
NYUC1=((1-Ws_Opt(end)/Ws_Opt(1))/(1-Ws(end)/Ws(1)))*ts(end) ;%One-Rotation w/ Feed
NYUC2=((1-Ws_Opt_Inf(end)/Ws_Opt_Inf(1))/(1-Ws_Inf(end)/Ws_Inf(1)))*ts_Inf(end); %Infinite Horizon w/ Feed

NYUC1_No=((1-Ws_Opt(end)/Ws_Opt(1))/(1-Ws_No(end)/Ws_No(1)))*ts(end) ;%One-Rotation w/o Feed
NYUC2_No=((1-Ws_Opt_Inf(end)/Ws_Opt_Inf(1))/(1-Ws_No_Inf(end)/Ws_No_Inf(1)))*ts_Inf(end); %Infinite Horizon w/o Feed

%Number of rotations during convergence
NRDC1=NYUC1./ts_Opt(end); %One-Rotation w/ Feed
NRDC2=NYUC2./ts_Opt_Inf(end); %Infinite Rotation w/ Feed

NRDC1_No=NYUC1_No./ts_Opt(end); %One-Rotation w/ Feed
NRDC2_No=NYUC2_No./ts_Opt_Inf(end); %Infinite Rotation w/ Feed

%Ratio of rotation lengths
%%% Column 5 of Table
RRL1_No=ts_Opt(end)./ts_No(end);
RRL1_Feed=ts_Opt(end)./ts(end);
RRL2_No=ts_Opt_Inf(end)./ts_No_Inf(end);
RRL2_Feed=ts_Opt_Inf(end)./ts_Inf(end);

%Human health impact of divergence

HHI1=((Is(end)-Is(1))*H).*NRDC1; %One rotation w/ Feed
HHI2=((Is_Inf(end)-Is_Inf(1))*H).*NRDC2; %Infinite Horizon w/ Feed

HHI1_No=((Is_No(end)-Is_No(1))*H).*NRDC1_No; %One rotation w/o Feed
HHI2_No=((Is_No_Inf(end)-Is_No_Inf(1))*H).*NRDC2_No; %Infinite Horizon w/o Feed


HHI1_Opt=((Is_Opt(end)-Is_Opt(1))*H).*NRDC1; %One rotation w/ Feed
HHI2_Opt=((Is_Opt_Inf(end)-Is_Opt_Inf(1))*H).*NRDC2; %Infinite Horizon w/ Feed 

HHI1_Opt_No=((Is_Opt(end)-Is_Opt(1))*H).*NRDC1_No; %One rotation w/o Feed
HHI2_Opt_No=((Is_Opt_Inf(end)-Is_Opt_Inf(1))*H).*NRDC2_No; %Infinite Horizon w/o Feed 



%%% Column 4
HHI_One=round(HHI1-HHI1_Opt);
HHI_Inf=round(HHI2-HHI2_Opt); 

HHI_One_No=round(HHI1_No-HHI1_Opt_No);
HHI_Inf_No=round(HHI2_No-HHI2_Opt_No);
    
    
end





        %Paper Figure for dynamics   
if false 
        
fig=figure
subplot(221)
plot(ts_Opt*365,Profits_Opt,'LineWidth',3); hold on
plot(ts_OptFeed*365,Profits_OptFeed,'LineWidth',3); hold on
plot(ts_Health.*365,Profits_Health,'LineWidth',3);hold on
plot(ts_HealthFeed.*365,Profits_HealthFeed,'LineWidth',3);
plot(ts_No*365,Profits_No,'LineWidth',3); hold on
plot(ts*365,Profits,'LineWidth',3);

title({'One-Rotation Horizon','(A)'},'FontSize', 16)
xlim([0 365*ts_Opt(end)+10])
ylim([-500 max(Profits)+100])    
plot([0 365*ts_Opt(end)+10],[0 0],'LineWidth',1, 'color','k')
ylabel({'Discounted Aquaculture','Profits (USD)'}, 'FontSize', 14)
        
subplot(222)
plot(ts_Opt_Inf*365,Profits_Opt_Inf,'LineWidth',3); hold on
plot(ts_OptFeed_Inf*365,Profits_OptFeed_Inf,'LineWidth',3); hold on
plot(ts_No_Inf*365,Profits_No_Inf,'LineWidth',3);
plot(ts_Inf*365,Profits_Inf,'LineWidth',3);


title({'Infinite Horizon','(B)'},'FontSize', 16)
%xlim([0 365*ts(end)+10])
ylim([-500 max(Profits)+100])    
plot([0 365*ts(end)+10],[0 0],'LineWidth',1, 'color','k')
%ylabel({'Discounted Profits','(USD)'}, 'FontSize', 14)
        
 
subplot(223)
p1=plot(ts_Opt*365,Ws_Opt./Ws_Opt(1),'LineWidth',3);hold on
p2=plot(ts_OptFeed*365,Ws_OptFeed./Ws_OptFeed(1),'LineWidth',3);hold on
p3=plot(ts_No*365,Ws_No./Ws_No(1),'LineWidth',3);
p4=plot(ts*365,Ws./Ws(1),'LineWidth',3);
title('(C)','FontSize', 16)
xlim([0 365*ts_Opt(end)+10])
ylim([0 1])
ylabel({'Infected Snails','(as prop. of Steady State)'}, 'FontSize', 14)
    legend1=legend([p1 p2 p3 p4],{'Societal Optimum w/o Feed~~~','Societal Optimum w/ Feed~~~','Aquaculture Optimum w/o Feed~~~','Aquaculture Optimum w/ Feed'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');



subplot(224)
plot(ts_Opt_Inf*365,Ws_Opt_Inf./Ws_Opt_Inf(1),'LineWidth',3); hold on
plot(ts_OptFeed_Inf*365,Ws_OptFeed_Inf./Ws_OptFeed_Inf(1),'LineWidth',3); hold on
plot(ts_No_Inf*365,Ws_No_Inf./Ws_No_Inf(1),'LineWidth',3,'LineStyle','--');
plot(ts_Inf*365,Ws_Inf./Ws_Inf(1),'LineWidth',3);
title('(D)','FontSize', 16)
xlim([0 365*ts(end)+10])
ylim([0 1])
%ylabel({'Infected Snails','(as prop. of Steady State)'}, 'FontSize', 14)


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    
    
    %saveas(gcf,'Aquaculture.png'); hold off 

end


   
    
end





%%% Things to delete:
if false

    
    
% Attack rate of Prawns on Non-Snail Prey
%Solved with Nset=40, Gauss collocation points in 6 or 11 steps
if false
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

AlphaY=linspace(alphaY*0.5,alphaY*1.5,6); %21 worked well May 13, 2021 %31 worked wwell, May 13 2021

Convergence_Time_AlphaY_NoPolicy=zeros(1,length(AlphaY));
Difference_Infections_AlphaY_NoPolicy=zeros(1,length(AlphaY));
Rotation_Length_Ratio_AlphaY_NoPolicy= zeros(1,length(AlphaY));
Profits_AlphaY_NoPolicy=zeros(1,length(AlphaY));
RotLen_AlphaY_NoPolicy=zeros(1,length(AlphaY));
Average_Cost_AlphaY_NoPolicy=zeros(1,length(AlphaY));

Convergence_Time_AlphaY_RotLen=zeros(1,length(AlphaY));
Difference_Infections_AlphaY_RotLen=zeros(1,length(AlphaY));
Rotation_Length_Ratio_AlphaY_RotLen= zeros(1,length(AlphaY));
Profits_AlphaY_RotLen=zeros(1,length(AlphaY));
RotLen_AlphaY_RotLen=zeros(1,length(AlphaY));
Average_Cost_AlphaY_RotLen=zeros(1,length(AlphaY));


Convergence_Time_AlphaY_FeedSea=zeros(1,length(AlphaY));
Difference_Infections_AlphaY_FeedSea=zeros(1,length(AlphaY));
Rotation_Length_Ratio_AlphaY_FeedSea= zeros(1,length(AlphaY));
Profits_AlphaY_FeedSea=zeros(1,length(AlphaY));
RotLen_AlphaY_FeedSea=zeros(1,length(AlphaY));
Average_Cost_AlphaY_FeedSea=zeros(1,length(AlphaY));


    for i=1:length(AlphaY)
        
        alphaY=AlphaY(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        GUESS=[];
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

GUESS=[];

        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_AlphaY_NoPolicy(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaY_NoPolicy(i)./ts_No(end);
        Rotation_Length_Ratio_AlphaY_NoPolicy(i)= ts(end)./ts_No(end);
               
               
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaY_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_AlphaY_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_AlphaY_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_NoPolicy_5000=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);       
        
        Profits_AlphaY_NoPolicy(i)=Profits(end);
        RotLen_AlphaY_NoPolicy(i)=ts(end);
        Average_Cost_AlphaY_NoPolicy(i)=(Profits(end))./ts(end);
        
        
        end
        
        %%% Minimum Rotation Lenght
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        GUESS=[]; 
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end


[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

  
GUESS=[];

        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end
        
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_AlphaY_RotLen(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaY_RotLen(i)./ts_No(end);      
        Rotation_Length_Ratio_AlphaY_RotLen(i)= ts(end)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaY_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_AlphaY_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;

        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_AlphaY_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_RotLen_5000=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);
        
        Profits_AlphaY_RotLen(i)=Profits(end);
        RotLen_AlphaY_RotLen(i)=ts(end);
        Average_Cost_AlphaY_RotLen(i)=round((Average_Cost_AlphaY_NoPolicy(i)-Profits(end)./ts(end)),0)./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_RotLen_5000));
        %Average_Cost_AlphaY_RotLen(i)=(round(Average_Cost_AlphaY_NoPolicy(i)-Profits(end),1))./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_RotLen_5000));
        
        
        
        end
        
        %%% Limiting Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        GUESS=[];
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end


[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

GUESS=[];

        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end        
        
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_AlphaY_FeedSea(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaY_FeedSea(i)./ts_No(end);
        Rotation_Length_Ratio_AlphaY_FeedSea(i)= ts(end)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaY_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_AlphaY_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_AlphaY_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_FeedSea_5000=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);
        
        Profits_AlphaY_FeedSea(i)=Profits(end);
        RotLen_AlphaY_FeedSea(i)=ts(end);
        Average_Cost_AlphaY_FeedSea(i)=round((Average_Cost_AlphaY_NoPolicy(i)-Profits(end)./ts(end)),0)./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_FeedSea_5000));
        if isnan(Average_Cost_AlphaY_FeedSea(i))==1
           Average_Cost_AlphaY_FeedSea(i)=0   ;
        end
        
        %Average_Cost_AlphaY_FeedSea(i)=(round(Average_Cost_AlphaY_NoPolicy(i)-Profits(end),1%))./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_FeedSea_5000));
         
        
        end
        
    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()


%%%% Figure
if true
fig=figure
subplot(231)
yyaxis left
plot(AlphaY,Convergence_Time_AlphaY_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
plot([alphaY alphaY],[0 2.5], 'Color', 'k')
%plot([alphaY alphaY],[0 1.5], 'Color', 'k')
%
yyaxis right
plot(AlphaY,Rotation_Length_Ratio_AlphaY_NoPolicy,'LineWidth',3); hold on
ylim([0.2 1.05])
%ylim([0.4 1.05])

subplot(232)
yyaxis left
plot(AlphaY,Convergence_Time_AlphaY_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([alphaY alphaY],[0 2.5], 'Color', 'k')
%plot([alphaY alphaY],[0 1.5], 'Color', 'k')

%
yyaxis right
plot(AlphaY,Rotation_Length_Ratio_AlphaY_RotLen,'LineWidth',3); hold on
ylim([0.2 1.05])
%ylim([0.4 1.05])
%
subplot(233)
yyaxis left
plot(AlphaY,Convergence_Time_AlphaY_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
plot([alphaY alphaY],[0 2.5], 'Color', 'k')
%plot([alphaY alphaY],[0 1.5], 'Color', 'k')
%
yyaxis right
plot(AlphaY,Rotation_Length_Ratio_AlphaY_FeedSea,'LineWidth',3); hold on
ylim([0.2 1.05])
%ylim([0.4 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
yyaxis left
plot(AlphaY,-Difference_Infections_AlphaY_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
plot([alphaY alphaY],[0 4000], 'Color', 'k')
%plot([alphaY alphaY],[0 15000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_NoPolicy,'LineWidth',3);
ylim([0 125])

subplot(235)
yyaxis left
plot(AlphaY,-Difference_Infections_AlphaY_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
plot([alphaY alphaY],[0 4000], 'Color', 'k')
%plot([alphaY alphaY],[0 15000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_RotLen,'LineWidth',3); hold on
plot(AlphaY,Average_Cost_AlphaY_RotLen,'LineWidth',3)
%ylim([0 250])
ylim([0 125])

subplot(236)
yyaxis left
plot(AlphaY,-Difference_Infections_AlphaY_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
plot([alphaY alphaY],[0 4000], 'Color', 'k')
%plot([alphaY alphaY],[0 15000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_FeedSea,'LineWidth',3); hold on
plot(AlphaY,Average_Cost_AlphaY_FeedSea,'LineWidth',3)
ylim([0 125])
%ylim([0 250])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Attack Rate of Prawns on Non-Snail Prey'}, 'FontSize', 14);
%saveas(gcf,'Attack_NonSnail_TypeII_Infinite.png'); hold off    
end

  
end

    
    
% Handling Time of Prawns on Supplemental Feed
if false

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

THU=linspace(ThU*0.25,ThU,21); %21 worked well May 13, 2021 %31 worked wwell, May 13 2021

Convergence_Time_THU_NoPolicy=zeros(1,length(THU));
Difference_Infections_THU_NoPolicy=zeros(1,length(THU));
Rotation_Length_Ratio_THU_NoPolicy= zeros(1,length(THU));
Average_Profits_THU_NoPolicy=zeros(1,length(THU));
Average_Cost_THU_NoPolicy=zeros(1,length(THU));

Convergence_Time_THU_RotLen=zeros(1,length(THU));
Difference_Infections_THU_RotLen=zeros(1,length(THU));
Rotation_Length_Ratio_THU_RotLen= zeros(1,length(THU));
Average_Profits_THU_RotLen=zeros(1,length(THU));
Average_Cost_THU_RotLen=zeros(1,length(THU));


Convergence_Time_THU_FeedSea=zeros(1,length(THU));
Difference_Infections_THU_FeedSea=zeros(1,length(THU));
Rotation_Length_Ratio_THU_FeedSea= zeros(1,length(THU));
Average_Profits_THU_FeedSea=zeros(1,length(THU));
Average_Cost_THU_FeedSea=zeros(1,length(THU));


    for i=1:length(THU)
        
        ThU=THU(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        GUESS=[];
        CASE=1;
        
if i>=2
    
    ThU=THU(i-1);

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);



if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end
end



[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

GUESS=[];

        CASE=2;
if i>=2
    
    ThU=THU(i-1);
        
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end
end


[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_THU_NoPolicy(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_THU_NoPolicy(i)./ts_No(end);
        Rotation_Length_Ratio_THU_NoPolicy(i)= ts(end)./ts_No(end);
               
               
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_THU_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_THU_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_THU_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_NoPolicy_5000=Convergence_Infections_WO-Convergence_Infections_W;       
        
        Average_Profits_THU_NoPolicy(i)=(Profits(end)-Profits_No(end))./(-Difference_Infections_NoPolicy_5000);
        Average_Cost_THU_NoPolicy(i)=Profits(end);
        
        
        end
        
        %%% Minimum Rotation Lenght
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        GUESS=[]; 
        CASE=1;

 if i>=2
    
    ThU=THU(i-1);
[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end
 end

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

  
GUESS=[];

        CASE=2;

 if i>=2
    
    ThU=THU(i-1);
    
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end
 end
 
 
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_THU_RotLen(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_THU_RotLen(i)./ts_No(end);      
        Rotation_Length_Ratio_THU_RotLen(i)= ts(end)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_THU_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_THU_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;

        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_THU_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_RotLen_5000=Convergence_Infections_WO-Convergence_Infections_W;
        
        Average_Profits_THU_RotLen(i)=(Profits(end)-Profits_No(end))./(-Difference_Infections_RotLen_5000);
        Average_Cost_THU_RotLen(i)=(Average_Cost_THU_NoPolicy(i)-Profits(end))./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_RotLen_5000));
        
        
        
        end
        
        %%% Limiting Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        GUESS=[];
        CASE=1;
        
if i>=2
    
    ThU=THU(i-1);

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);



if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end
end

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

GUESS=[];

        CASE=2;

 if i>=2
    
    ThU=THU(i-1);
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end        
   
 end

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_THU_FeedSea(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_THU_FeedSea(i)./ts_No(end);
        Rotation_Length_Ratio_THU_FeedSea(i)= ts(end)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_THU_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_THU_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_THU_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_FeedSea_5000=Convergence_Infections_WO-Convergence_Infections_W;
        
        Average_Profits_THU_FeedSea(i)=(Profits(end)-Profits_No(end))./(-Difference_Infections_FeedSea_5000);
        Average_Cost_THU_FeedSea(i)=(Average_Cost_THU_NoPolicy(i)-Profits(end))./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_FeedSea_5000));
        
         
        
        end
        
    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()


%%%% Figure
if true
fig=figure
subplot(231)
yyaxis left
plot(THU,Convergence_Time_THU_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
plot([ThU ThU],[0 3], 'Color', 'k')
%
yyaxis right
plot(THU,Rotation_Length_Ratio_THU_NoPolicy,'LineWidth',3); hold on
ylim([0.45 1.05])

subplot(232)
yyaxis left
plot(THU,Convergence_Time_THU_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([ThU ThU],[0 3], 'Color', 'k')
%
yyaxis right
plot(THU,Rotation_Length_Ratio_THU_RotLen,'LineWidth',3); hold on
ylim([0.45 1.05])

subplot(233)
yyaxis left
plot(THU,Convergence_Time_THU_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
plot([ThU ThU],[0 3], 'Color', 'k')
%
yyaxis right
plot(THU,Rotation_Length_Ratio_THU_FeedSea,'LineWidth',3); hold on
ylim([0.45 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
yyaxis left
plot(THU,-Difference_Infections_THU_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
plot([ThU ThU],[0 3000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_NoPolicy,'LineWidth',3);
ylim([0 3000])

subplot(235)
yyaxis left
plot(THU,-Difference_Infections_THU_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
plot([ThU ThU],[0 3000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_RotLen,'LineWidth',3); hold on
plot(THU,Average_Cost_THU_RotLen,'LineWidth',3)
ylim([0 3000])

subplot(236)
yyaxis left
plot(THU,-Difference_Infections_THU_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
plot([ThU ThU],[0 3000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_FeedSea,'LineWidth',3); hold on
plot(THU,Average_Cost_THU_FeedSea,'LineWidth',3)
%ylim([0 100])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Prawns Handling Time of Feed'}, 'FontSize', 14);
% saveas(gcf,'Attack_NonSnail.png'); hold off    
end

  
end


% Handling Time of Prawns on Supplemental Feed
if false

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

THY=linspace(ThY*0.001,ThY,6); %21 worked well May 13, 2021 %31 worked wwell, May 13 2021

Convergence_Time_THY_NoPolicy=zeros(1,length(THY));
Difference_Infections_THY_NoPolicy=zeros(1,length(THY));
Rotation_Length_Ratio_THY_NoPolicy= zeros(1,length(THY));
Average_Profits_THY_NoPolicy=zeros(1,length(THY));
Average_Cost_THY_NoPolicy=zeros(1,length(THY));

Convergence_Time_THY_RotLen=zeros(1,length(THY));
Difference_Infections_THY_RotLen=zeros(1,length(THY));
Rotation_Length_Ratio_THY_RotLen= zeros(1,length(THY));
Average_Profits_THY_RotLen=zeros(1,length(THY));
Average_Cost_THY_RotLen=zeros(1,length(THY));


Convergence_Time_THY_FeedSea=zeros(1,length(THY));
Difference_Infections_THY_FeedSea=zeros(1,length(THY));
Rotation_Length_Ratio_THY_FeedSea= zeros(1,length(THY));
Average_Profits_THY_FeedSea=zeros(1,length(THY));
Average_Cost_THY_FeedSea=zeros(1,length(THY));


    for i=1:length(THY)
        
        ThY=THY(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        GUESS=[];
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

GUESS=[];

        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_THY_NoPolicy(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_THY_NoPolicy(i)./ts_No(end);
        Rotation_Length_Ratio_THY_NoPolicy(i)= ts(end)./ts_No(end);
               
               
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_THY_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_THY_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_THY_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_NoPolicy_5000=Convergence_Infections_WO-Convergence_Infections_W;       
        
        Average_Profits_THY_NoPolicy(i)=(Profits(end)-Profits_No(end))./(-Difference_Infections_NoPolicy_5000);
        Average_Cost_THY_NoPolicy(i)=Profits(end);
        
        
        end
        
        %%% Minimum Rotation Lenght
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        GUESS=[]; 
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end


[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

  
GUESS=[];

        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end
        
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_THY_RotLen(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_THY_RotLen(i)./ts_No(end);      
        Rotation_Length_Ratio_THY_RotLen(i)= ts(end)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_THY_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_THY_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;

        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_THY_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_RotLen_5000=Convergence_Infections_WO-Convergence_Infections_W;
        
        Average_Profits_THY_RotLen(i)=(Profits(end)-Profits_No(end))./(-Difference_Infections_RotLen_5000);
        Average_Cost_THY_RotLen(i)=(Average_Cost_THY_NoPolicy(i)-Profits(end))./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_RotLen_5000));
        
        
        
        end
        
        %%% Limiting Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        GUESS=[];
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

if true
GUESS(CASE).I=Is_No;
GUESS(CASE).W=Ws_No;
GUESS(CASE).X=Xs_No;
GUESS(CASE).L=Ls_No;
GUESS(CASE).P=Ps_No;
end


[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

GUESS=[];

        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

if true
GUESS(CASE).I=Is;
GUESS(CASE).W=Ws;
GUESS(CASE).X=Xs;
GUESS(CASE).L=Ls;
GUESS(CASE).P=Ps;
GUESS(CASE).U=Us;
end        
        
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_THY_FeedSea(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_THY_FeedSea(i)./ts_No(end);
        Rotation_Length_Ratio_THY_FeedSea(i)= ts(end)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_THY_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_THY_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_THY_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_FeedSea_5000=Convergence_Infections_WO-Convergence_Infections_W;
        
        Average_Profits_THY_FeedSea(i)=(Profits(end)-Profits_No(end))./(-Difference_Infections_FeedSea_5000);
        Average_Cost_THY_FeedSea(i)=(Average_Cost_THY_NoPolicy(i)-Profits(end))./(-(Difference_Infections_NoPolicy_5000-Difference_Infections_FeedSea_5000));
        
         
        
        end
        
    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()


%%%% Figure
if true
fig=figure
subplot(231)
yyaxis left
plot(THY,Convergence_Time_THY_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
plot([ThY ThY],[0 3], 'Color', 'k')
%
yyaxis right
plot(THY,Rotation_Length_Ratio_THY_NoPolicy,'LineWidth',3); hold on
ylim([0.45 1.05])

subplot(232)
yyaxis left
plot(THY,Convergence_Time_THY_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([ThY ThY],[0 3], 'Color', 'k')
%
yyaxis right
plot(THY,Rotation_Length_Ratio_THY_RotLen,'LineWidth',3); hold on
ylim([0.45 1.05])

subplot(233)
yyaxis left
plot(THY,Convergence_Time_THY_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
plot([ThY ThY],[0 3], 'Color', 'k')
%
yyaxis right
plot(THY,Rotation_Length_Ratio_THY_FeedSea,'LineWidth',3); hold on
ylim([0.45 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
yyaxis left
plot(THY,-Difference_Infections_THY_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
plot([ThY ThY],[0 3000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_NoPolicy,'LineWidth',3);
ylim([0 3000])

subplot(235)
yyaxis left
plot(THY,-Difference_Infections_THY_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
plot([ThY ThY],[0 3000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_RotLen,'LineWidth',3); hold on
plot(THY,Average_Cost_THY_RotLen,'LineWidth',3)
%ylim([0 100])

subplot(236)
yyaxis left
plot(THY,-Difference_Infections_THY_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
plot([ThY ThY],[0 3000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_FeedSea,'LineWidth',3); hold on
plot(THY,Average_Cost_THY_FeedSea,'LineWidth',3)
%ylim([0 100])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Prawns Handling Time of Non Snail Prey'}, 'FontSize', 14);
% saveas(gcf,'Attack_NonSnail.png'); hold off    
end

  
end









% Attack rate of Prawns on Non-Snail Prey
if false

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

AlphaY=linspace(alphaY*0.5,alphaY*1.5,31);

Convergence_Time_AlphaY_NoPolicy=zeros(1,length(AlphaY));
Difference_Infections_AlphaY_NoPolicy=zeros(1,length(AlphaY));

Convergence_Time_AlphaY_RotLen=zeros(1,length(AlphaY));
Difference_Infections_AlphaY_RotLen=zeros(1,length(AlphaY));

Convergence_Time_AlphaY_FeedSea=zeros(1,length(AlphaY));
Difference_Infections_AlphaY_FeedSea=zeros(1,length(AlphaY));


    for i=1:length(AlphaY)
        
        alphaY=AlphaY(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
   SchistoAquaculture_Feed(Topt_No,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

        
        Convergence_Time_AlphaY_NoPolicy(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaY_NoPolicy(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaY_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_AlphaY_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        end
        
        %%% Minimum Rotation Lenght
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
   SchistoAquaculture_Feed(Topt_No,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

        
        Convergence_Time_AlphaY_RotLen(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaY_RotLen(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaY_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_AlphaY_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;
        end
        
        %%% Limiting Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
   SchistoAquaculture_Feed(Topt_No,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

        
        Convergence_Time_AlphaY_FeedSea(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaY_FeedSea(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaY_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_AlphaY_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        end
        
    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()

fig=figure
subplot(231)
plot(AlphaY,Convergence_Time_AlphaY_NoPolicy,'LineWidth',3); hold on
title({'(A)','No Policy'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
plot([alphaY alphaY],[0 3], 'Color', 'k')

subplot(232)
plot(AlphaY,Convergence_Time_AlphaY_RotLen,'LineWidth',3); hold on
title({'(B)',' Minimum Lenght'},'FontSize', 16)
plot([alphaY alphaY],[0 3], 'Color', 'k')

subplot(233)
plot(AlphaY,Convergence_Time_AlphaY_FeedSea,'LineWidth',3); hold on
title({'(C)', 'Limiting Feeding'},'FontSize', 16)
plot([alphaY alphaY],[0 3], 'Color', 'k')


subplot(234)
plot(AlphaY,-Difference_Infections_AlphaY_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
plot([alphaY alphaY],[0 4000], 'Color', 'k')

subplot(235)
plot(AlphaY,-Difference_Infections_AlphaY_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
plot([alphaY alphaY],[0 4000], 'Color', 'k')

subplot(236)
plot(AlphaY,-Difference_Infections_AlphaY_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
plot([alphaY alphaY],[0 4000], 'Color', 'k')

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Attack Rate of Prawns on Non-Snail Prey'}, 'FontSize', 14);
  saveas(gcf,'Attack_NonSnail.png'); hold off    
  
  
end

   
% Attack rate of Prawns on Feed
if false

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

Eta=linspace(eta*0.25,eta*1.75,31);

Convergence_Time_Eta_NoPolicy=zeros(1,length(Eta));
Difference_Infections_Eta_NoPolicy=zeros(1,length(Eta));

Convergence_Time_Eta_RotLen=zeros(1,length(Eta));
Difference_Infections_Eta_RotLen=zeros(1,length(Eta));

Convergence_Time_Eta_FeedSea=zeros(1,length(Eta));
Difference_Infections_Eta_FeedSea=zeros(1,length(Eta));


    for i=1:length(Eta)
        
        eta=Eta(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
   SchistoAquaculture_Feed(Topt_No,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

        
        Convergence_Time_Eta_NoPolicy(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_Eta_NoPolicy(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_Eta_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_Eta_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        end
        
        %%% Minimum Rotation Lenght
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
   SchistoAquaculture_Feed(Topt_No,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

        
        Convergence_Time_Eta_RotLen(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_Eta_RotLen(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_Eta_RotLen(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_Eta_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;
        end
        
        %%% Limiting Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        
        CASE=1;

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
   SchistoAquaculture_Feed(Topt_No,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);

        
        Convergence_Time_Eta_FeedSea(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_Eta_FeedSea(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_Eta_FeedSea(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_Eta_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        end
        
    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()

fig=figure
subplot(231)
plot(Eta,Convergence_Time_Eta_NoPolicy,'LineWidth',3); hold on
title({'(A)','No Policy'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
plot([eta eta],[0 4], 'Color', 'k')

subplot(232)
plot(Eta,Convergence_Time_Eta_RotLen,'LineWidth',3); hold on
title({'(B)',' Minimum Lenght'},'FontSize', 16)
plot([eta eta],[0 4], 'Color', 'k')

subplot(233)
plot(Eta,Convergence_Time_Eta_FeedSea,'LineWidth',3); hold on
title({'(C)', 'Limiting Feeding'},'FontSize', 16)
plot([eta eta],[0 4], 'Color', 'k')


subplot(234)
plot(Eta,-Difference_Infections_Eta_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
plot([eta eta],[0 3000], 'Color', 'k')

subplot(235)
plot(Eta,-Difference_Infections_Eta_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
plot([eta eta],[0 3000], 'Color', 'k')

subplot(236)
plot(Eta,-Difference_Infections_Eta_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
plot([eta eta],[0 3000], 'Color', 'k')

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Parameter Mimicking Variations in Search Costs','Associated with Inanimate Food'}, 'FontSize', 14);
    saveas(gcf,'Search_Costs.png'); hold off    
  
  
end





   

% Density-dependent reduction in growth
if false
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()


G=linspace(0,g*2,11);

Convergence_Time_G=zeros(1,length(G));
Difference_Infections_G=zeros(1,length(G));


    for i=1:length(G)
        
        g=G(i);
        
        CASE=1;

        [ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
            SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

        [ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
            SchistoAquaculture_Feed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    
        
        Convergence_Time_G(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_G(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_G(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_G(i)=Convergence_Infections_WO-Convergence_Infections_W;

    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()
  
if true
fig=figure
subplot(211)
plot(G,Convergence_Time_G,'LineWidth',3); hold on
title('(A) Convergence Time','FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
plot([g g],[1.9 2.4], 'Color', 'k')


subplot(212)
plot(G,-Difference_Infections_G,'LineWidth',3); hold
title('(B) Impact on Human Health of Divergence', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
plot([g g],[1500 2000], 'Color', 'k')

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Parameter Mimicking Density-Dependent','Reduction in Growth with Biomass'}, 'FontSize', 14);
     saveas(gcf,'Density_Dependent_Growth.png'); hold off    
  

end


end

% Discount Rate
if false
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()
   

R=linspace(0,r*2,6);

Convergence_Time_R=zeros(1,length(R));
Difference_Infections_R=zeros(1,length(R));


    for i=1:length(R)
        
        r=R(i);
        
        CASE=1;

        [ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
            SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

    
        CASE=2;

        [ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
            SchistoAquaculture_Feed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


    
        
        Convergence_Time_R(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_R(i)./ts_No(end);
        
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_R(i);
        Convergence_Infections_WO=((Is_No(end)-Is_No(1))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        
        Difference_Infections_R(i)=Convergence_Infections_WO-Convergence_Infections_W;

    end
    

if true
fig=figure
subplot(211)
plot(R,Convergence_Time_R,'LineWidth',3);
title('(A) Convergence Time','FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)

subplot(212)
plot(R,-Difference_Infections_R,'LineWidth',3)
title('(B) Impact on Human Health of Divergence', 'FontSize', 16)
ylabel({'Additional','Human Cases'}, 'FontSize', 14)

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Discount Rate'}, 'FontSize', 14);
   

end


end


%Figure for last rotation simulated
if false
    
    if true
fig=figure
subplot(223)
plot(ts_No*365,Ws_No./Ws_No(1),'LineWidth',3);
title('(A) Infected Snails w/o Feed')
xlim([0 365*ts_No(end)+10])
ylim([0 1])
subplot(221)
plot(ts_No*365,Profits_No,'LineWidth',3);
title('(C) Discounted Profits (USD) w/o Feed')
xlim([0 365*ts_No(end)+10])
ylim([0 max(Profits)+100])
subplot(224)
plot(ts*365,Ws./Ws(1),'LineWidth',3);
title('(B) Infected Snails w/ Feed')
xlim([0 365*ts_No(end)+10])
ylim([0 1])
subplot(222)
plot(ts*365,Profits,'LineWidth',3);
title('(D) Discounted Profits (USD) w/ Feed')
xlim([0 365*ts_No(end)+10])
ylim([0 max(Profits)+100])

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
   
    
%saveas(gcf,'Aquaculture.png'); hold off    
    
%saveas(gcf,'Aquaculture_Infinite.png'); hold off

    end

    
    if false
%%%%%%%%%%%%%%
fig=figure
subplot(221)
plot(ts_No*365,Ws_No./Ws_No(1)); hold on
plot(ts*365,Ws./Ws(1));
title({'Infected Snails','(as Proportion of Initial Level)'})
xlim([0 365*ts_No(end)+10])
ylim([0 1])
%
subplot(223)
plot(ts_No*365,Profits_No); hold on
plot(ts*365,Profits);
title('Discounted Profits (USD)')
xlim([0 365*ts_No(end)+10])
ylim([0 max(Profits)+100])
%
H=1e6;
subplot(222)
p1=plot(ts_No*365,(Is_No-Is_No(1)).*H); hold on
p2=plot(ts*365,(Is-Is(1)).*H); hold on
xlim([0 365*ts_No(end)+10])
%ylim([0 1])
title({'Reduction in Infected Humans','per Million People'})
    legend1=legend([p1 p2],{'Without Feeding~~~','With Feeding'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.320840017407273 0.492857142857143 0.384517125449869 0.0445661601554252],...
    'Orientation','horizontal',...
    'Interpreter','latex');

%
subplot(224)
plot(ts_No*365,psiWs_No); hold on
plot(ts*365,psiWs); hold on
title({'Per-Capita Prawn','Predation of Snails'})
xlim([0 365*ts_No(end)+10])

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
   
saveas(gcf,'Aquaculture2.png'); hold off

    figure
    plot(ts,Us);
    
%ylim([0 max(Profits)+100])

    end

end














%Figures
if false
  

if true %Prawns
figure
subplot(221)
plot(ts*365,Ls);
title('Average Lenght of Prawns (mm)')
xlim([0 365*ts(end)])
subplot(222)
plot(ts*365,Ps);
title('Number of Prawns')
xlim([0 365*ts(end)])
subplot(223)
plot(ts*365,Bs)
title('Average Body Size of Prawns (g)')
xlim([0 365*ts(end)])
subplot(224)
plot(ts*365,Omegas);
title('Biomass of Prawns (kg)')
xlim([0 365*ts(end)])
end

if true %Human and Snails
figure
subplot(221)
plot(ts*365,Is);
title('Infected Humans')
subplot(222)
plot(ts*365,Ws);
title('Infected Snails')
subplot(223)
plot(ts*365,Xs);
title('Healthy Snails')
subplot(224)
plot(ts*365,Ns);
title('Total Snails')
end

if OBJ==1 %Paper Figure

figure
subplot(211)
plot(ts*365,Ws./Ws(1));
title('Infected Snails')
xlim([0 365*ts(end)])
ylim([0 1])
subplot(212)
plot(ts*365,Profits);
title('One Shot Discounted Profits (USD)')
xlim([0 365*ts(end)])
ylim([0 max(Profits)+100])
    
elseif OBJ==2
    
figure
subplot(211)
plot(ts*365,Ws./Ws(1));
title('Infected Snails')
xlim([0 365*ts(end)])
subplot(212)
plot(ts*365,Profits);
title('Discounted Profits of One Rotation for Infinite Horizon  (USD)')
xlim([0 365*ts(end)])
ylim([0 max(Profits)+100])
    
end



if true %Attack Rate and Prawn/Snail Ratio

figure
subplot(311)
plot(ts*365, alphaNs)
title('Attack Rate')
subplot(312)
plot(ts*365,Ratios);
title('Prawn/Snail Ratio')
subplot(313)
plot(ts*365,Ths);
title('Handling Time')

end

if CASE==2
    
    figure
    subplot(221)
    plot(ts*365,Us);
    title('Supplemental Feed')
    subplot(223)
    plot(ts*365,Us./Omegas);
    title('Supplemental Feed/Prawn Biomass')
    subplot(222)
    plot(ts*365,Costs);
    title('Cumulative Costs')
    subplot(224)
    plot(ts*365,exp(-r.*ts).*(cU.*Us));
    title('Current Costs (qualitative only)')

end

    
if true %Per-Capita Attack rates

figure
plot(ts*365,psiWs); hold on
plot(ts*365,psiXs); hold on
plot(ts*365,psiNs); hold on
xlim([0 365*ts(end)])
%ylim([0 1.1])

end

end

end

    %New  Policy Analyses: Feed Conversion Effiency with Social objective
if false
    
    [beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    
    H=1e6;
    
AlphaU=linspace(alphaU*0.5,alphaU*1.5,3); 


   
%%% For figure   
Convergence_Time_AlphaU_NoPolicy=zeros(1,length(AlphaU));
Difference_AvoidedInfections_AlphaU_NoPolicy=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_NoPolicy= zeros(1,length(AlphaU));
Average_Cost_AlphaU_NoPolicy=zeros(1,length(AlphaU));

TotalProfits_Convergence_NoPolicy=zeros(1,length(AlphaU));
TotalAvoidedCases_Convergence_NoPolicy=zeros(1,length(AlphaU));

Convergence_Time_AlphaU_RotLen=zeros(1,length(AlphaU));
Difference_AvoidedInfections_AlphaU_RotLen=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_RotLen= zeros(1,length(AlphaU));
Average_Cost_AlphaU_RotLen=zeros(1,length(AlphaU));


Convergence_Time_AlphaU_FeedSea=zeros(1,length(AlphaU));
Difference_AvoidedInfections_AlphaU_FeedSea=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_FeedSea= zeros(1,length(AlphaU));
Average_Cost_AlphaU_FeedSea=zeros(1,length(AlphaU));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(AlphaU)
    
    alphaU=AlphaU(i);

    %%%%%%%%%%%%%%%%%%%
    %%% No Policy
    %%%%%%%%%%%%%%%%%%%
    if true
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   %%%%Providing an Initial Guess     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        OBJ=2;  %Private Infinite Horizon
        CASE=1; %No feed
        GUESS=[]; %No Guess

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   %%%%Enf of Providing an Initial Guess     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        POLICY=1; %Need to specify Policy in Feed Case, here POLICY=1 means no policy
        OBJ=6; % Societal Objective, Infinite Horizon
        CASE=2; %With feed
        
        %Guess
        if true
        GUESS(CASE).I=Is_No;
        GUESS(CASE).W=Ws_No;
        GUESS(CASE).X=Xs_No;
        GUESS(CASE).L=Ls_No;
        GUESS(CASE).P=Ps_No;
        GUESS(CASE).U=0;
        end
    
        [ts_Social, Us_Social, Is_Social, Ws_Social, Xs_Social, Ns_Social, Ls_Social, Ps_Social, Omegas_Social, Profits_Social, psiWs_Social, psiXs_Social, Avoided_HCosts_Social, alphaNs_Social, Ths_Social, ks_Social, Results_Social] = ...
    SchistoAquaculture_Feed(0,0,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);

%Saving Dynamics
if true
Dynamics_Social_Time(:,i)=ts_Social;
Dynamics_Social_Us(:,i)=Us_Social;
Dynamics_Social_Is(:,i)=Is_Social;
Dynamics_Social_Ws(:,i)=Ws_Social;
Dynamics_Social_Xs(:,i)=Xs_Social;
Dynamics_Social_Ns(:,i)=Ns_Social;
Dynamics_Social_Ls(:,i)=Ls_Social;
Dynamics_Social_Ps(:,i)=Ps_Social;
Dynamics_Social_OM(:,i)=Omegas_Social;
Dynamics_Social_Pi(:,i)=Profits_Social;
Dynamics_Social_PsiW(:,i)=psiWs_Social;
Dynamics_Social_Avoided(:,i)=Avoided_HCosts_Social;
Dynamics_Social_alphaNs(:,i)=alphaNs_Social;
Dynamics_Social_Th(:,i)=Ths_Social;
Dynamics_Social_ks(:,i)=ks_Social;
end



  
        OBJ=2; %Private Objective, Infinite Horizon
       %Guess       
        if true
        GUESS(CASE).I=Is_Social;
        GUESS(CASE).W=Ws_Social;
        GUESS(CASE).X=Xs_Social;
        GUESS(CASE).L=Ls_Social;
        GUESS(CASE).P=Ps_Social;
        GUESS(CASE).U=Us_Social;
        end

        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Social(end),Ws_Social(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);

%Saving Dynamics
if true
 Dynamics_NoPolicy_Time(:,i)=ts_Private;
 Dynamics_NoPolicy_Us(:,i)=Us_Private;
 Dynamics_NoPolicy_Is(:,i)=Is_Private;
 Dynamics_NoPolicy_Ws(:,i)=Ws_Private;
 Dynamics_NoPolicy_Xs(:,i)=Xs_Private;
 Dynamics_NoPolicy_Ns(:,i)=Ns_Private;
 Dynamics_NoPolicy_Ls(:,i)=Ls_Private;
 Dynamics_NoPolicy_Ps(:,i)=Ps_Private;
 Dynamics_NoPolicy_OM(:,i)=Omegas_Private;
 Dynamics_NoPolicy_Pi(:,i)=Profits_Private;
 Dynamics_NoPolicy_PsiW(:,i)=psiWs_Private;
 Dynamics_NoPolicy_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_NoPolicy_alphaNs(:,i)=alphaNs_Private;
 Dynamics_NoPolicy_Th(:,i)=Ths_Private;
 Dynamics_NoPolicy_ks(:,i)=ks_Private;
end


           %Saving for Policy Analysis Figure
        if false
    a=ts_Social(min(find(round(Ws_Social,1)==0))); %Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_AlphaU_NoPolicy(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Social(1)-Is_Social(end))./ts_Social(end)).*Convergence_Time_AlphaU_NoPolicy(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_NoPolicy(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_AlphaU_NoPolicy(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_AlphaU_NoPolicy(i)=ts_Private(end)./ts_Social(end); %Ratio of rotation lengths
    
    e1=floor(Convergence_Time_AlphaU_NoPolicy(i)./ts_Social(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_AlphaU_NoPolicy(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
    e3= (Convergence_Time_AlphaU_NoPolicy(i)./ts_Social(end) - e1).*ts_Social(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_AlphaU_NoPolicy(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
    e5= Profits_Social(end).*e1 + Profits_Social(max(find(ts_Social<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
    e7= ((Is_Social(1)-Is_Social(end))./ts_Social(end)).*Convergence_Time_AlphaU_NoPolicy(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_NoPolicy(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    
    
    % %
    TotalProfits_Convergence_NoPolicy(i)=e6;
    % %
    TotalAvoidedCases_Convergence_NoPolicy(i)=e8;
    
    % %
    Average_Cost_AlphaU_NoPolicy(i) = (e6 - e5)./(e7 - e8); %Average cost of an averted case if we force the health case
        end
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%
    %%% Minimum Rotation Lenght
    %%%%%%%%%%%%%%%%%%%   
    if true
        
        POLICY=2;
        
        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Social(end),Ws_Social(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);


%Saving Dynamics
if true
 Dynamics_RotLen_Time(:,i)=ts_Private;
 Dynamics_RotLen_Us(:,i)=Us_Private;
 Dynamics_RotLen_Is(:,i)=Is_Private;
 Dynamics_RotLen_Ws(:,i)=Ws_Private;
 Dynamics_RotLen_Xs(:,i)=Xs_Private;
 Dynamics_RotLen_Ns(:,i)=Ns_Private;
 Dynamics_RotLen_Ls(:,i)=Ls_Private;
 Dynamics_RotLen_Ps(:,i)=Ps_Private;
 Dynamics_RotLen_OM(:,i)=Omegas_Private;
 Dynamics_RotLen_Pi(:,i)=Profits_Private;
 Dynamics_RotLen_PsiW(:,i)=psiWs_Private;
 Dynamics_RotLen_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_RotLen_alphaNs(:,i)=alphaNs_Private;
 Dynamics_RotLen_Th(:,i)=Ths_Private;
 Dynamics_RotLen_ks(:,i)=ks_Private;
end
        
        %Saving for Policy Analysis Figure
        if false
    a=ts_Social(min(find(round(Ws_Social,1)==0))) ;%Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_AlphaU_RotLen(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Social(1)-Is_Social(end))./ts_Social(end)).*Convergence_Time_AlphaU_RotLen(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_RotLen(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_AlphaU_RotLen(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_AlphaU_RotLen(i)=ts_Private(end)./ts_Social(end); %Ratio of rotation lengths
    
   % e1=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Social(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
   % e3= (Convergence_Time_AlphaU_RotLen(i)./ts_Social(end) - e1).*ts_Social(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_AlphaU_RotLen(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
   % e5= Profits_Social(end).*e1 + Profits_Social(max(find(ts_Social<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
   % e7= ((Is_Social(1)-Is_Social(end))./ts_Social(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    

    % %
    Average_Cost_AlphaU_RotLen(i) = (TotalProfits_Convergence_NoPolicy(i) - e6)./(e8 - TotalAvoidedCases_Convergence_NoPolicy(i)); %Average cost of an averted case if we force the health case
        end
    
    
        
    end
    

    
    %%%%%%%%%%%%%%%%%%%
    %%% Limiting Feeding Season
    %%%%%%%%%%%%%%%%%%%   
    if true
        
        POLICY=3;
        
        [ts_Private, Us_Private, Is_Private, Ws_Private, Xs_Private, Ns_Private, Ls_Private, Ps_Private, Omegas_Private, Profits_Private, psiWs_Private, psiXs_Private, Avoided_HCosts_Private, alphaNs_Private, Ths_Private, ks_Private, Results_Private] = ...
    SchistoAquaculture_Feed(ts_Social(end),Ws_Social(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS);


%Saving Dynamics
if true
 Dynamics_FeedSea_Time(:,i)=ts_Private;
 Dynamics_FeedSea_Us(:,i)=Us_Private;
 Dynamics_FeedSea_Is(:,i)=Is_Private;
 Dynamics_FeedSea_Ws(:,i)=Ws_Private;
 Dynamics_FeedSea_Xs(:,i)=Xs_Private;
 Dynamics_FeedSea_Ns(:,i)=Ns_Private;
 Dynamics_FeedSea_Ls(:,i)=Ls_Private;
 Dynamics_FeedSea_Ps(:,i)=Ps_Private;
 Dynamics_FeedSea_OM(:,i)=Omegas_Private;
 Dynamics_FeedSea_Pi(:,i)=Profits_Private;
 Dynamics_FeedSea_PsiW(:,i)=psiWs_Private;
 Dynamics_FeedSea_Avoided(:,i)=Avoided_HCosts_Private;
 Dynamics_FeedSea_alphaNs(:,i)=alphaNs_Private;
 Dynamics_FeedSea_Th(:,i)=Ths_Private;
 Dynamics_FeedSea_ks(:,i)=ks_Private;
end

        %Saving for Policy Analysis Figure
        if false
    a=ts_Social(min(find(round(Ws_Social,1)==0))) ;%Time period W(t) reaches the vicinity of zero in health case
    b=max(find(ts_Private<=a)); %Collocation point equivalent to time period for which W(t) reaches the vicinity of zero in health case
    c=1./(1-(Ws_Private(b)./Ws_Private(1))); %Factor multiplying the time it takes for which W(t) to reache the vicinity of zero in health case
    
    % %
    Convergence_Time_AlphaU_FeedSea(i)=ts_Private(b).*c; % Convergence time (in years)
    
    d1=((Is_Social(1)-Is_Social(end))./ts_Social(end)).*Convergence_Time_AlphaU_FeedSea(i).*H; %Avoided infections (per 1M) during convergence in Health Case
    d2=((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_FeedSea(i).*H; %Avoided infections (per 1M) during convergence in Private Case
    
    % %
    Difference_AvoidedInfections_AlphaU_FeedSea(i)=d1-d2;  %Difference in avoided infections
    % %
    Rotation_Length_Ratio_AlphaU_FeedSea(i)=ts_Private(end)./ts_Social(end); %Ratio of rotation lengths
    
   % e1=floor(Convergence_Time_AlphaU_RotLen(i)./ts_Social(end));  %Number of full Health Case rotations during convergence
    e2=floor(Convergence_Time_AlphaU_FeedSea(i)./ts_Private(end)); %Number of full Private Case rotations during convergence
   % e3= (Convergence_Time_AlphaU_RotLen(i)./ts_Social(end) - e1).*ts_Social(end) ;%Amount of time (in years) spent in the last rotation before convergence in Health case
    e4= (Convergence_Time_AlphaU_FeedSea(i)./ts_Private(end) - e2).*ts_Private(end); %Amount of time (in years) spent in the last rotation before convergence in Private case  
   % e5= Profits_Social(end).*e1 + Profits_Social(max(find(ts_Social<=e3))); %Profits earned during convergence in the health case
    e6= Profits_Private(end).*e2 + Profits_Private(max(find(ts_Private<=e4))); %Profits earned during convergence in the private case 
   % e7= ((Is_Social(1)-Is_Social(end))./ts_Social(end)).*Convergence_Time_AlphaU_RotLen(i).*5000; %Avoided infections (per 5K) during convergence in Health Case
    e8= ((Is_Private(1)-Is_Private(end))./ts_Private(end)).*Convergence_Time_AlphaU_FeedSea(i).*5000; %Avoided infections (per 5K) during convergence in Private Case
    

    % %
    Average_Cost_AlphaU_FeedSea(i) = (TotalProfits_Convergence_NoPolicy(i) - e6)./(e8 - TotalAvoidedCases_Convergence_NoPolicy(i)); %Average cost of an averted case if we force the health case
    
        end
    


    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure for dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
    
set(0, 'DefaultLineLineWidth', 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infect Snails by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(411)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(412)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(413)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(C)  Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(414)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)

    legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Infected Snails','(as Proportion of Steady State)'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infect Snails by degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(311)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ws(:,1)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ws(:,1)./x0ic(2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ws(:,1)./x0ic(2),'--'); hold on
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(312)
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ws(:,2)./x0ic(2),'LineWidth',3); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ws(:,2)./x0ic(2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ws(:,2)./x0ic(2),'--'); hold on
title({'(B) Base Case Feed Efficiency'},'FontSize', 16)
%
subplot(313)
p1=plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ws(:,3)./x0ic(2),'LineWidth',3); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ws(:,3)./x0ic(2),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ws(:,3)./x0ic(2),'--'); hold on
title({'(C) Higher Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Infected Snails','(as Proportion of Steady State)'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Feed by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(411)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Us(:,1)); hold on
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Us(:,2)); hold on
plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Us(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(412)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Us(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Us(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Us(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(413)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Us(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Us(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Us(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(414)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Us(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Us(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Us(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Feed by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(311)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Us(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Us(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Us(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Us(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(312)
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Us(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Us(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Us(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Us(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(313)
p1=plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Us(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Us(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Us(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Us(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Aquaculture Profits by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Pi(:,1)); hold on
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Pi(:,2)); hold on
plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Pi(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Pi(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Pi(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Pi(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Pi(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Pi(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Pi(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Pi(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Pi(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Pi(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Aquaculture Profits by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Pi(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Pi(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Pi(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Pi(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Pi(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Pi(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Pi(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Pi(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Pi(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Pi(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Pi(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Pi(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end


if false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infection by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Is(:,1)); hold on
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Is(:,2)); hold on
plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Is(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Is(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Is(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Is(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Is(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Is(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Is(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Is(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Is(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Is(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Infection by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Is(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Is(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Is(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Is(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Is(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Is(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Is(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Is(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Is(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Is(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Is(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Is(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Prawn Predation of Snails by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_PsiW(:,1)); hold on
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_PsiW(:,2)); hold on
plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_PsiW(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_PsiW(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_PsiW(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_PsiW(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_PsiW(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_PsiW(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_PsiW(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_PsiW(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_PsiW(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_PsiW(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Prawn Predation by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_PsiW(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_PsiW(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_PsiW(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_PsiW(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_PsiW(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_PsiW(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_PsiW(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_PsiW(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_PsiW(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_PsiW(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_PsiW(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_PsiW(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Number of Prawns by Objective/Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(141)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Ps(:,1)); hold on
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Ps(:,2)); hold on
plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Ps(:,3)); hold on
title({'(A) Health Objective'},'FontSize', 16)
%
subplot(142)
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ps(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ps(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ps(:,3)); hold on
title({'(B) Profit-Maximization'},'FontSize', 16)
%
subplot(143)
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ps(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ps(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ps(:,3)); hold on
title({'(C) Profit-Maximization & Minimum Rotation Length'},'FontSize', 16)
%
subplot(144)
p1=plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ps(:,1)); hold on
p2=plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ps(:,2)); hold on
p3=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ps(:,3)); hold on
title({'(D) Profit-Maximization & Limited Feeding Season'},'FontSize', 16)


   legend1=legend([p1 p2 p3],{'Lower Feed Efficiency~~~','Base Case~~~','Higher Feed Efficiency~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Amount of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Number of Prawns by Degree of Feed Conversion Efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
fig=figure
subplot(131)
plot(365*Dynamics_Social_Time(:,1),Dynamics_Social_Ps(:,1)); hold on
plot(365*Dynamics_NoPolicy_Time(:,1),Dynamics_NoPolicy_Ps(:,1)); hold on
plot(365*Dynamics_RotLen_Time(:,1),Dynamics_RotLen_Ps(:,1),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,1),Dynamics_FeedSea_Ps(:,1),'--'); hold on
%ylim([0 2.5])
title({'(A) Lower Feed Efficiency'},'FontSize', 16)
%
subplot(132)
plot(365*Dynamics_Social_Time(:,2),Dynamics_Social_Ps(:,2)); hold on
plot(365*Dynamics_NoPolicy_Time(:,2),Dynamics_NoPolicy_Ps(:,2)); hold on
plot(365*Dynamics_RotLen_Time(:,2),Dynamics_RotLen_Ps(:,2),'--'); hold on
plot(365*Dynamics_FeedSea_Time(:,2),Dynamics_FeedSea_Ps(:,2),'--'); hold on
title({'(B) Base Case Efficiency'},'FontSize', 16)
%
subplot(133)
p1=plot(365*Dynamics_Social_Time(:,3),Dynamics_Social_Ps(:,3)); hold on
p2=plot(365*Dynamics_NoPolicy_Time(:,3),Dynamics_NoPolicy_Ps(:,3)); hold on
p3=plot(365*Dynamics_RotLen_Time(:,3),Dynamics_RotLen_Ps(:,3),'--'); hold on
p4=plot(365*Dynamics_FeedSea_Time(:,3),Dynamics_FeedSea_Ps(:,3),'--'); hold on
title({'(C) Higher Feed Efficiency'},'FontSize', 16)

    legend1=legend([p1 p2 p3 p4],{'Health Objective~~~','Profit-Maximization~~~','PMAX \& Min. Rot. Len.~~~','PMAX \& Lim. Feed. Sea.~~~'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,{'Quantity of Feed'}, 'FontSize', 16);
    xlabel(han,'Time (days)', 'FontSize', 16);
    %title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    
end


end


    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure for policy evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
fig=figure
subplot(231)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_NoPolicy,'LineWidth',3); hold on
%ylim([0.45 1.05])
%ylim([0.3 1.05])
%plot([alphaU alphaU],[0 1.05], 'Color', 'k')

subplot(232)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([alphaU alphaU],[0 6.1], 'Color', 'k')
%plot([alphaU alphaU],[0 40], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_RotLen,'LineWidth',3); hold on
%ylim([0.45 1.05])
%ylim([0.3 1.05])

subplot(233)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
plot([alphaU alphaU],[0 6.1], 'Color', 'k')
%plot([alphaU alphaU],[0 40], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_FeedSea,'LineWidth',3); hold on
%ylim([0.45 1.05])
%ylim([0.3 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
yyaxis left
plot(AlphaU,Difference_AvoidedInfections_AlphaU_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
%plot([alphaU alphaU],[0 10000], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Average_Cost_AlphaU_NoPolicy,'LineWidth',3);
%ylim([0 230])
%ylim([0 375])

subplot(235)
yyaxis left
plot(AlphaU,Difference_AvoidedInfections_AlphaU_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
%plot([alphaU alphaU],[0 10000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_RotLen,'LineWidth',3); hold on
plot(AlphaU,Average_Cost_AlphaU_RotLen,'LineWidth',3)
%ylim([0 230])
%ylim([0 375])
%
subplot(236)
yyaxis left
plot(AlphaU,Difference_AvoidedInfections_AlphaU_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
%plot([alphaU alphaU],[0 10000], 'Color', 'k')
%plot([alphaU alphaU],[0 25000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_FeedSea,'LineWidth',3); hold on
plot(AlphaU,Average_Cost_AlphaU_FeedSea,'LineWidth',3)
%ylim([0 230])
%ylim([0 375])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
    title(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
% saveas(gcf,'Attack_Feed_025.png'); hold off    
end

    

end


    % Policy Analyses
if false

    H=1e6; %Measure of infection, i.e. Nb of infection per H population


    
 
% Attack rate of Prawns on Feed
%Solved with Nset=5, Gauss collocation points in 5 steps
if true

[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

AlphaU=linspace(alphaU*0.5,alphaU*1.5,3); %21 worked well May 13, 2021 %31 worked wwell, May 13 2021

Convergence_Time_AlphaU_NoPolicy=zeros(1,length(AlphaU));
Difference_Infections_AlphaU_NoPolicy=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_NoPolicy= zeros(1,length(AlphaU));
Profits_AlphaU_NoPolicy=zeros(1,length(AlphaU));
Average_Cost_AlphaU_NoPolicy=zeros(1,length(AlphaU));

Convergence_Time_AlphaU_RotLen=zeros(1,length(AlphaU));
Difference_Infections_AlphaU_RotLen=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_RotLen= zeros(1,length(AlphaU));
Profits_AlphaU_RotLen=zeros(1,length(AlphaU));
Average_Cost_AlphaU_RotLen=zeros(1,length(AlphaU));


Convergence_Time_AlphaU_FeedSea=zeros(1,length(AlphaU));
Difference_Infections_AlphaU_FeedSea=zeros(1,length(AlphaU));
Rotation_Length_Ratio_AlphaU_FeedSea= zeros(1,length(AlphaU));
Profits_AlphaU_FeedSea=zeros(1,length(AlphaU));
Average_Cost_AlphaU_FeedSea=zeros(1,length(AlphaU));


    for i=1:length(AlphaU)
        
        alphaU=AlphaU(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;

        OBJ=5;

        CASE=1;
        GUESS=[];
        
[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


        CASE=2;
        
        if true
            GUESS(CASE).I=Is_No;
            GUESS(CASE).W=Ws_No;
            GUESS(CASE).X=Xs_No;
            GUESS(CASE).L=Ls_No;
            GUESS(CASE).P=Ps_No;
            GUESS(CASE).U=0;
        end
        


[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_AlphaU_NoPolicy(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaU_NoPolicy(i)./ts_No(end);
        Rotation_Length_Ratio_AlphaU_NoPolicy(i)= ts(end)./ts_No(end);
               
        if OBJ==1
        
        Convergence_Infections_W=((Is(end)-Is(1))*H);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_AlphaU_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000);
            
        elseif OBJ==2
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaU_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_AlphaU_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_AlphaU_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        end
        
        
        if OBJ==1
        Difference_Infections_AlphaU_NoPolicy_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W);       
        elseif OBJ==2
        Difference_Infections_AlphaU_NoPolicy_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);       
        end
        
        
        Profits_AlphaU_NoPolicy(i)=Profits(end);
        
        if OBJ==1
        Average_Cost_AlphaU_NoPolicy(i)=Profits(end);  
        elseif OBJ==2
        Average_Cost_AlphaU_NoPolicy(i)=Profits(end)./ts(end);
        end
        
        end
        
        %%% Minimum Rotation Lenght
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        

        CASE=1;
        GUESS=[];
[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


        CASE=2;
        
        if true
            GUESS(CASE).I=Is;
            GUESS(CASE).W=Ws;
            GUESS(CASE).X=Xs;
            GUESS(CASE).L=Ls;
            GUESS(CASE).P=Ps;
            GUESS(CASE).U=Us;
        end
        
        
[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_AlphaU_RotLen(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaU_RotLen(i)./ts_No(end);      
        Rotation_Length_Ratio_AlphaU_RotLen(i)= ts(end)./ts_No(end);
        
        if OBJ==1
        Convergence_Infections_W=((Is(end)-Is(1))*H);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H);     
        Difference_Infections_AlphaU_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;

        
        Convergence_Infections_W=((Is(end)-Is(1))*5000);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000);       
            
        elseif OBJ==2
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaU_RotLen(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_AlphaU_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;

        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_AlphaU_RotLen(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000)*Nb_Rotations_NoFeed_DuringConvergence;   
        end
        
        if OBJ==1
        Difference_Infections_AlphaU_RotLen_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W);            
        elseif OBJ==2
        Difference_Infections_AlphaU_RotLen_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);
        end
        
        
        Profits_AlphaU_RotLen(i)=Profits(end);
        
        if OBJ==1
        Average_Cost_AlphaU_RotLen(i)=round((Average_Cost_AlphaU_NoPolicy(i)-Profits(end)),0)./((round(Difference_Infections_AlphaU_NoPolicy_5000(i),0)-round(Difference_Infections_AlphaU_RotLen_5000(i),0)));    
        elseif OBJ==2
        Average_Cost_AlphaU_RotLen(i)=round((Average_Cost_AlphaU_NoPolicy(i)-Profits(end)./ts(end)),0)./((round(Difference_Infections_AlphaU_NoPolicy_5000(i),0)-round(Difference_Infections_AlphaU_RotLen_5000(i),0)));
        end
        
        
        
        
        end
        
        %%% Limiting Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        

        CASE=1;
        GUESS=[];
        
[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


        CASE=2;

        if true
            GUESS(CASE).I=Is;
            GUESS(CASE).W=Ws;
            GUESS(CASE).X=Xs;
            GUESS(CASE).L=Ls;
            GUESS(CASE).P=Ps;
            GUESS(CASE).U=Us;
        end

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_AlphaU_FeedSea(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_AlphaU_FeedSea(i)./ts_No(end);
        Rotation_Length_Ratio_AlphaU_FeedSea(i)= ts(end)./ts_No(end);

        if OBJ==1
        Convergence_Infections_W=((Is(end)-Is(1))*H);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H);
        Difference_Infections_AlphaU_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000);            
            
        elseif OBJ==2
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_AlphaU_FeedSea(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_AlphaU_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_AlphaU_FeedSea(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        end
        
        
        
        if OBJ==1
        Difference_Infections_AlphaU_FeedSea_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W);    
        elseif OBJ==2
        Difference_Infections_AlphaU_FeedSea_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);
        end
        
        
        Profits_AlphaU_FeedSea(i)=Profits(end);
        
        if OBJ==1
        Average_Cost_AlphaU_FeedSea(i)=round((Average_Cost_AlphaU_NoPolicy(i)-Profits(end)),0)./((round(Difference_Infections_AlphaU_NoPolicy_5000(i),0)-round(Difference_Infections_AlphaU_FeedSea_5000(i),0)));    
        elseif OBJ==2
        Average_Cost_AlphaU_FeedSea(i)=round((Average_Cost_AlphaU_NoPolicy(i)-Profits(end)./ts(end)),0)./((round(Difference_Infections_AlphaU_NoPolicy_5000(i),0)-round(Difference_Infections_AlphaU_FeedSea_5000(i),0)));
        end
        
        if isnan(Average_Cost_AlphaU_FeedSea(i))==1
           Average_Cost_AlphaU_FeedSea(i)=0   ;
        end
              
        
        end
        
    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()


%%%% Figure
if true
fig=figure
subplot(231)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
%plot([alphaU alphaU],[0 40], 'Color', 'k')
plot([alphaU alphaU],[0 6.1], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_NoPolicy,'LineWidth',3); hold on
%ylim([0.45 1.05])
ylim([0.3 1.05])

subplot(232)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([alphaU alphaU],[0 6.1], 'Color', 'k')
%plot([alphaU alphaU],[0 40], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_RotLen,'LineWidth',3); hold on
%ylim([0.45 1.05])
ylim([0.3 1.05])

subplot(233)
yyaxis left
plot(AlphaU,Convergence_Time_AlphaU_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
plot([alphaU alphaU],[0 6.1], 'Color', 'k')
%plot([alphaU alphaU],[0 40], 'Color', 'k')
%
yyaxis right
plot(AlphaU,Rotation_Length_Ratio_AlphaU_FeedSea,'LineWidth',3); hold on
%ylim([0.45 1.05])
ylim([0.3 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
yyaxis left
plot(AlphaU,Difference_Infections_AlphaU_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
plot([alphaU alphaU],[0 10000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_NoPolicy,'LineWidth',3);
ylim([0 230])
%ylim([0 375])

subplot(235)
yyaxis left
plot(AlphaU,Difference_Infections_AlphaU_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
plot([alphaU alphaU],[0 10000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_RotLen,'LineWidth',3); hold on
plot(AlphaU,Average_Cost_AlphaU_RotLen,'LineWidth',3)
ylim([0 230])
%ylim([0 375])
%
subplot(236)
yyaxis left
plot(AlphaU,Difference_Infections_AlphaU_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
plot([alphaU alphaU],[0 10000], 'Color', 'k')
%plot([alphaU alphaU],[0 25000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_FeedSea,'LineWidth',3); hold on
plot(AlphaU,Average_Cost_AlphaU_FeedSea,'LineWidth',3)
ylim([0 230])
%ylim([0 375])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Attack Rate of Prawns on Supplemental Feed'}, 'FontSize', 14);
% saveas(gcf,'Attack_Feed_025.png'); hold off    
end

  
end



%Exponent of the Holling type functional response
if true
    
    [beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters();

    
    N=linspace(1.75,2.25,3);

Convergence_Time_N_NoPolicy=zeros(1,length(N));
Difference_Infections_N_NoPolicy=zeros(1,length(N));
Rotation_Length_Ratio_N_NoPolicy= zeros(1,length(N));
Profits_N_NoPolicy=zeros(1,length(N));
Average_Cost_N_NoPolicy=zeros(1,length(N));

Convergence_Time_N_RotLen=zeros(1,length(N));
Difference_Infections_N_RotLen=zeros(1,length(N));
Rotation_Length_Ratio_N_RotLen= zeros(1,length(N));
Profits_N_RotLen=zeros(1,length(N));
Average_Cost_N_RotLen=zeros(1,length(N));


Convergence_Time_N_FeedSea=zeros(1,length(N));
Difference_Infections_N_FeedSea=zeros(1,length(N));
Rotation_Length_Ratio_N_FeedSea= zeros(1,length(N));
Profits_N_FeedSea=zeros(1,length(N));
Average_Cost_N_FeedSea=zeros(1,length(N));


    for i=1:length(N)
        
        n=N(i);
        
        %%% No Policy
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=1;
        
        CASE=1;
        GUESS=[];

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


        CASE=2;
                
       if true
            GUESS(CASE).I=Is_No;
            GUESS(CASE).W=Ws_No;
            GUESS(CASE).X=Xs_No;
            GUESS(CASE).L=Ls_No;
            GUESS(CASE).P=Ps_No;
            GUESS(CASE).U=0;
        end
        
        

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_N_NoPolicy(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_N_NoPolicy(i)./ts_No(end);
        Rotation_Length_Ratio_N_NoPolicy(i)= ts(end)./ts_No(end);
               
        if OBJ==1
        Convergence_Infections_W=((Is(end)-Is(1))*H);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H);
        Difference_Infections_N_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000);    
            
        elseif OBJ==2
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_N_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_N_NoPolicy(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_N_NoPolicy(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        end
        
        
        if OBJ==1
        Difference_Infections_N_NoPolicy_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W);               
        elseif OBJ==2
        Difference_Infections_N_NoPolicy_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);       
        end
        
        
        Profits_N_NoPolicy(i)=Profits(end);
        
        if OBJ==1
        Average_Cost_N_NoPolicy(i)=Profits(end);    
        elseif OBJ==2
        Average_Cost_N_NoPolicy(i)=Profits(end)./ts(end);
        end
        
        end
        
        %%% Minimum Rotation Lenght
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=2;
        

        CASE=1;
        GUESS=[];

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


        CASE=2;
       if true
            GUESS(CASE).I=Is;
            GUESS(CASE).W=Ws;
            GUESS(CASE).X=Xs;
            GUESS(CASE).L=Ls;
            GUESS(CASE).P=Ps;
            GUESS(CASE).U=Us;
        end
        

[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_N_RotLen(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_N_RotLen(i)./ts_No(end);      
        Rotation_Length_Ratio_N_RotLen(i)= ts(end)./ts_No(end);

        if OBJ==1
        Convergence_Infections_W=((Is(end)-Is(1))*H);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H);     
        Difference_Infections_N_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;

        
        Convergence_Infections_W=((Is(end)-Is(1))*5000);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000);    
        
        elseif OBJ==2
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_N_RotLen(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H)*Nb_Rotations_NoFeed_DuringConvergence;     
        Difference_Infections_N_RotLen(i)=Convergence_Infections_WO-Convergence_Infections_W;

        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_N_RotLen(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000)*Nb_Rotations_NoFeed_DuringConvergence;    
        end
        
        if OBJ==1
        Difference_Infections_N_RotLen_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W);    
        elseif OBJ==2
        Difference_Infections_N_RotLen_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);
        end
        
        
        Profits_N_RotLen(i)=Profits(end);
       
        if OBJ==1
        Average_Cost_N_RotLen(i)=round((Average_Cost_N_NoPolicy(i)-Profits(end)),0)./((round(Difference_Infections_N_NoPolicy_5000(i),0)-round(Difference_Infections_N_RotLen_5000(i),0)));   
        elseif OBJ==2
        Average_Cost_N_RotLen(i)=round((Average_Cost_N_NoPolicy(i)-Profits(end)./ts(end)),0)./((round(Difference_Infections_N_NoPolicy_5000(i),0)-round(Difference_Infections_N_RotLen_5000(i),0)));
        end
        
        
        if isnan(Average_Cost_N_RotLen(i))==1 
           Average_Cost_N_RotLen(i)=0   ;
        end
        if isinf(Average_Cost_N_RotLen(i))==1
           Average_Cost_N_RotLen(i)=0   ;
        end
        
        
        
        end
        
        %%% Limiting Feeding Season
        %%%%%%%%%%%%%%%%%%%
        if true
        POLICY=3;
        

        CASE=1;
        GUESS=[];

[ts_No, Topt_No, Is_No, Ws_No, Xs_No, Ns_No, Ls_No, Ps_No, Bs_No, Omegas_No, Profits_No, psiWs_No, psiXs_No, psiNs_No, alphaNs_No, Ratios_No, Ths_No, ks_No, Results_No] = ...
   SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS);


        CASE=2;
        
        if true
            GUESS(CASE).I=Is_No;
            GUESS(CASE).W=Ws_No;
            GUESS(CASE).X=Xs_No;
            GUESS(CASE).L=Ls_No;
            GUESS(CASE).P=Ps_No;
            GUESS(CASE).U=0;
        end


[ts, Topt, Us, Costs, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, ks_I, ks_E, Results] = ...
    SchistoAquaculture_Feed(Topt_No,Ws_No(end),Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,POLICY,GUESS);


        
        Convergence_Time_N_FeedSea(i)=((1-Ws_No(end)/Ws_No(1))/(1-Ws(end)/Ws(1)))*ts(end);
        Nb_Rotations_NoFeed_DuringConvergence=Convergence_Time_N_FeedSea(i)./ts_No(end);
        Rotation_Length_Ratio_N_FeedSea(i)= ts(end)./ts_No(end);
        
        if OBJ==1
        Convergence_Infections_W=((Is(end)-Is(1))*H);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H);
        Difference_Infections_N_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000);
        
                    
        elseif OBJ==2
        Convergence_Infections_W=((Is(end)-Is(1))*H)*Convergence_Time_N_FeedSea(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*H)*Nb_Rotations_NoFeed_DuringConvergence;
        Difference_Infections_N_FeedSea(i)=Convergence_Infections_WO-Convergence_Infections_W;
        
        Convergence_Infections_W=((Is(end)-Is(1))*5000)*Convergence_Time_N_FeedSea(i);
        Convergence_Infections_WO=((Is_No(1)-Is_No(end))*5000)*Nb_Rotations_NoFeed_DuringConvergence;
        end
        
        
        
        if OBJ==1
        Difference_Infections_N_FeedSea_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W);    
        elseif OBJ==2
        Difference_Infections_N_FeedSea_5000(i)=(Convergence_Infections_WO-Convergence_Infections_W)./ts(end);
        end
        
        
        Profits_N_FeedSea(i)=Profits(end)./ts(end);
        
        if OBJ==1
        Average_Cost_N_FeedSea(i)=round((Average_Cost_N_NoPolicy(i)-Profits(end)),0)./((round(Difference_Infections_N_NoPolicy_5000(i),0)-round(Difference_Infections_N_FeedSea_5000(i),0)));    
        elseif OBJ==2
        Average_Cost_N_FeedSea(i)=round((Average_Cost_N_NoPolicy(i)-Profits(end)./ts(end)),0)./((round(Difference_Infections_N_NoPolicy_5000(i),0)-round(Difference_Infections_N_FeedSea_5000(i),0)));
        end
        
       if isnan(Average_Cost_N_FeedSea(i))==1 
           Average_Cost_N_FeedSea(i)=0   ;
       end
       if isinf(Average_Cost_N_FeedSea(i))==1
           Average_Cost_N_FeedSea(i)=0   ;
        end
        
        
         
        
        end
        
    end
    
[beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()


%%%% Figure
if true
fig=figure
subplot(231)
yyaxis left
plot(N,Convergence_Time_N_NoPolicy,'LineWidth',3); hold on
title({'No Policy','(A)'},'FontSize', 16)
ylabel({'Number of Years','Until Convergence'}, 'FontSize', 14)
%plot([alphaU alphaU],[0 40], 'Color', 'k')
plot([n n],[0 5], 'Color', 'k')
%
yyaxis right
plot(N,Rotation_Length_Ratio_N_NoPolicy,'LineWidth',3); hold on
%ylim([0.45 1.05])
ylim([0.25 1.05])

subplot(232)
yyaxis left
plot(N,Convergence_Time_N_RotLen,'LineWidth',3); hold on
title({' Minimum Lenght','(B)'},'FontSize', 16)
plot([n n],[0 5], 'Color', 'k')
%plot([alphaU alphaU],[0 40], 'Color', 'k')
%
yyaxis right
plot(N,Rotation_Length_Ratio_N_RotLen,'LineWidth',3); hold on
%ylim([0.45 1.05])
ylim([0.25 1.05])

subplot(233)
yyaxis left
plot(N,Convergence_Time_N_FeedSea,'LineWidth',3); hold on
title({'Limiting Feeding','(C)'},'FontSize', 16)
plot([n n],[0 5], 'Color', 'k')
%plot([alphaU alphaU],[0 40], 'Color', 'k')

%
yyaxis right
plot(N,Rotation_Length_Ratio_N_FeedSea,'LineWidth',3); hold on
%ylim([0.45 1.05])
ylim([0.25 1.05])
ylabel({'Ratio Feed/No Feed of','the Rotation Length'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Length of Rotation w/ Feed}}{\textrm{Lenght of Rotation w/o Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

subplot(234)
yyaxis left
plot(N,Difference_Infections_N_NoPolicy,'LineWidth',3); hold on
title('(D)', 'FontSize', 16)
ylabel({'Additional Human','Cases per 1M People'}, 'FontSize', 14)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
plot([n n],[0 6500], 'Color', 'k')
ylim([0 6500])
%
yyaxis right
%plot(AlphaU,Average_Profits_NoPolicy,'LineWidth',3);
ylim([-20 80])
%ylim([0 375])

subplot(235)
yyaxis left
plot(N,Difference_Infections_N_RotLen,'LineWidth',3); hold on
title('(E) ', 'FontSize', 16)
%plot([alphaU alphaU],[0 3000], 'Color', 'k')
plot([n n],[0 6500], 'Color', 'k')
ylim([0 6500])
%
yyaxis right
%plot(AlphaU,Average_Profits_RotLen,'LineWidth',3); hold on
plot(N,Average_Cost_N_RotLen,'LineWidth',3)
ylim([-20 80])
%ylim([0 375])
%
subplot(236)
yyaxis left
plot(N,Difference_Infections_N_FeedSea,'LineWidth',3); hold on
title('(F) ', 'FontSize', 16)
plot([n n],[0 6500], 'Color', 'k')
ylim([0 6500])
%plot([alphaU alphaU],[0 25000], 'Color', 'k')
%
yyaxis right
%plot(AlphaU,Average_Profits_FeedSea,'LineWidth',3); hold on
plot(N,Average_Cost_N_FeedSea,'LineWidth',3)
ylim([-20 80])
%ylim([0 375])
ylabel({'Average Cost of', 'an Averted Case'}, 'FontSize', 14)
%ylabel({'$\frac{\textrm{Rotation Length Feed}}{\textrm{Rotation Length No Feed}}$'},'Interpreter','latex', 'FontSize', 12);
    

    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,{'Exponent of Holling''s Type III Functional Response '}, 'FontSize', 14);
 %saveas(gcf,'Holling_TypeIII.png'); hold off    
end



    
end


end



if false

%%% Predation Figures
Wp=linspace(0,1);
alphaNp=3;
Thp=25;
Np= 0.3896;
Up=4.5;
n=1;


psiW_wo=(alphaNp.*epsilon.*Wp.^n)./(1+alphaNp.*epsilon.*Thp.*Np.^n);

psiW_w=(alphaNp.*epsilon.*Wp.^n)./(1+alphaNp.*epsilon.*Thp.*Np.^n  + eta.*alphaU.*Thp.*Up.^n);


figure
plot(Wp,psiW_wo); hold on
plot(Wp,psiW_w);

end

