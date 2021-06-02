 function [ts, Us, Is, Ws, Xs, Ns, Ls, Ps, Omegas, Profits, psiWs, psiXs, Avoided_HCosts, alphaNs, Ths, ks, Results] = ...
    SchistoAquaculture_Feed(Topt_Health,Wno_T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,epsilon,n,cU,alphaU,ThU,CASE,OBJ,POLICY,GUESS)

if CASE==2 %% CASE== 2 : Aquaculture with feeding
    
toms t Tstar;
p=tomPhase('p', t, 0, Tstar, Nset,[],'gauss'); %gauss is what has been working (May 10, 2021) %'gauss' 'cheb' or 'fem1s' or 'fem1'
setPhase(p);
tomStates I W X L P
tomControls U

%Terminal Conditions
if POLICY==1 %No Policy
 cterm=final({Tstar>=1/365}); %L>=140 
 cterm2=final({}); 
 cterm3=final({});  %Tstar<=1
elseif POLICY==2 %Standard Rotation Length
 cterm=final({Tstar>=1/365}); 
 cterm2=final({Tstar>=Topt_Health}); 
 cterm3=final({}); 
%cterm3=final({Tstar<=1}); 
elseif POLICY==3 %Limiting feeding season
cterm=final({Tstar>=1/365}); %; 
cterm2={collocate(( U <= 1e3*(t<((Topt_Health)*(1/2))) ))};
cterm3=final({}); 
%cterm3=final({Tstar<=1}); 
elseif POLICY==4 % Tentative Policy penalizing change in Udot
    %Forced to reach W_T
cterm=final({Tstar>=1/365}); 
cterm2=final({});% cterm2={collocate(( U <= 1e3*(t<((Topt_Health)*(1/2))) ))};
cterm3=final({}); %cterm3=final({W<=Wno_T}); 
end


%%% Functions used in ODEs
N=X+W; %Total population of snails (N) equals the sum of the non-infected (X) and infected snails (W)
B=(aP.*L.^(bP)).*0.00001; % Mean Prawn body size (need to "*0.00001" to get in grams)
Omega=B.*P; % Biomass (i.e. mean prawn body size x number of prawns); needs to be in grams because "g" and "omega" are for grams
Ratio=B./((aN.*12.^bN)*.001);  %(need to "*0.001" to get in grams) %Biomass ratio prawn to snail; I chose 8mm snail (medium); goes from 4mm (small) to 8mm (medium) to 12mm (large); see Hoover et al.
alphaN=aM.*log(Ratio); %attack rate
Th=((th.*Ratio).^(-1)*365); %Handling time estimated from Sokolow et al. (13; see Hoover et al)

psiW=(alphaN.*epsilon.*W.^n)./(1+alphaN.*epsilon.*Th.*N.^n + alphaU.*ThU.*U.^n);%per-capita prawn predation rate on infected snails 
psiX=(alphaN.*epsilon.*X.^n)./(1+alphaN.*epsilon.*Th.*N.^n + alphaU.*ThU.*U.^n);%per-capita prawn predation rate on healthy snails 

k=(max(kMax,alphaN.*N) + alphaU.*sqrt(U))./(1+g.*Omega);

%k=(max(kMax,alphaN.*N) + alphaU.*U.^(0.2))./(1+g.*Omega);

 
%ODEs
ode = collocate({
    dot(I) == beta.*W.*(1-I) - gamma.*I ; %I dot
    dot(W) == lambda.*I.*X - delta.*W - psiW.*P;%    %W dot
    dot(X) == f.*X.*(1-N) - lambda.*I.*X - psiX.*P ;
    dot(L) == k.*(Linf-L);
    dot(P) == -P(muP.*B.^(-d) + omega.*Omega)});

cbb={icollocate(0 <= I)
     icollocate(0 <= W)
     icollocate(0 <= X)
     icollocate(0 <= L)
     icollocate(0 <= P)
     icollocate(0 <= Tstar)
     icollocate(0 <= U ) 
 initial(I == x0ic(1))
 initial(W == x0ic(2))
 initial(X == x0ic(3))
 initial(L == x0ic(4))
 initial(P == x0ic(5))}; 



 %%% Initial guess
    if isempty(GUESS)==1

    x0_guess = { icollocate({I == x0ic(1);
        W==x0ic(2)
        X==x0ic(3)
        L==x0ic(4)
        P==x0ic(5)
        U==0})}; 
    elseif isempty(GUESS)~=1  

     x0_guess = { icollocate({I == GUESS(CASE).I;
        W==GUESS(CASE).W;
        X==GUESS(CASE).X;
        L==GUESS(CASE).L;
        P==GUESS(CASE).P;});
        collocate(U==GUESS(CASE).U)};
    end

     
options=struct;
    if OBJ==1
        if POLICY==1
        options.name='Supplemental Feeding & One-Shot Horizon & No Policy';
        elseif POLICY==2
          options.name='Supplemental Feeding & One-Shot Horizon & Rotation Lenght';  
        elseif POLICY==3
        options.name='Supplemental Feeding & One-Shot Horizon & Feeding Season';    
        end
    elseif OBJ==2
        options.name='Supplemental Feeding & Infinite Horizon';
    elseif OBJ==3
        options.name='Societal Optimum & Feed & Infinite Horizon';
    elseif OBJ==4
        options.name='Societal Optimum & Feed & Infinite Horizon';
    elseif OBJ==5
        options.name='Health Optimum & Feed & One-Shot Horizon';
    elseif OBJ==6
        options.name='Health Optimum & Feed & Infinite Horizon';
    end
    
    options.solver = 'knitro'; % choose solver ('knitro' is other good one)
    
    Omega=Omega.*0.001; %in kilograms %multiplied by "0.001" to get kilograms because price is per kg
    
    %Penalizing Change in U
    if POLICY==4
        POLICY_Cost=1e6.*dot(U);
    else
        POLICY_Cost=0;
    end
    %Feeding cost parameter and cost function
    cU_quad=200; 
    feed= integrate(exp(-r.*t).*(cU_quad.*U.^2 + cU.*U));
    %feed= integrate(exp(-r.*t).*cU_quad.*U);
    
if OBJ==1 %Private: Single-Rotation Harvest 
objective= - exp(-r.*Tstar).*price.*final(Omega) + feed + cP.*initial(P); 
elseif OBJ==2 %Private: Infinite-horizon Harvest
objective= - (1./(exp(r.*Tstar)-1)).*(price.*final(Omega) - feed - cP.*initial(P)) + cP.*initial(P) + feed ;
elseif OBJ==3 %Health: Single-Rotation Harvest
objective=  cP.*initial(P) - cI.*initial(I) +  exp(-r.*Tstar).*cI.*final(I) + feed ; 
elseif OBJ==4 %Health: Infinite-Horizon Harvest
objective= - (1./(exp(r.*Tstar)-1)).*(cI.*(initial(I)-final(I)) - cP.*initial(P) - feed ) - cI.*initial(I) + cP.*initial(P) + feed ;
elseif OBJ==5
objective= - exp(-r.*Tstar).*price.*final(Omega) + cP.*initial(P) - cI.*initial(I) +  exp(-r.*Tstar).*cI.*final(I) + feed ; 
elseif OBJ==6 %Societal Optimum Infinite-horizon
objective= - (1./(exp(r.*Tstar)-1)).*(cI.*(initial(I)-final(I)) - cP.*initial(P) - feed  + price.*final(Omega)) - cI.*initial(I) + cP.*initial(P) + feed ; 
end


constr={ode,cbb,cterm,cterm2,cterm3};

    [solution,result]= ezsolve(objective,constr,x0_guess,options);

    [solution,result]= ezsolve(objective,constr,solution,options);
    
    %%% Changing solver if exitflag~=0
    counter=0;
    while result.ExitFlag~=0  && counter==0  % limits the lopp size

        options.solver = 'snopt';
    [solution,result]= ezsolve(objective,constr,x0_guess,options);
    %[solution,result]= ezsolve(objective,constr,solution,options);
    
        if result.ExitFlag~=0
            
    
        options.solver = 'npsol';
    [solution,result]= ezsolve(objective,constr,x0_guess,options);
    %[solution,result]= ezsolve(objective,constr,solution,options);
        end
        
        if result.ExitFlag~=0
            
        options.solver = 'knitro';
    [solution,result]= ezsolve(objective,constr,x0_guess,options);
     %[solution,result]= ezsolve(objective,constr,solution,options);
        end
        
        counter=counter+1;
    end
    
  if true  
       % [solution,result]= ezsolve(objective,constr,x0_guess,options);
    
    while result.ExitFlag~=0 && counter==1  % limits the lopp size

        options.solver = 'snopt';
    [solution,result]= ezsolve(objective,constr,solution,options);
    
        if result.ExitFlag~=0
            
    
        options.solver = 'npsol';
    [solution,result]= ezsolve(objective,constr,solution,options);
        end
        
        if result.ExitFlag~=0
            
        options.solver = 'knitro';
 
    [solution,result]= ezsolve(objective,constr,solution,options);
        end
        
        counter=counter+1;
    end
  end
  
  
    if result.ExitFlag==0 
        
        Results=result; %%% not sure about this
        ts=subs(icollocate(t), solution);
        Is=subs(icollocate(I), solution);
        Ws=subs(icollocate(W), solution);
        Xs=subs(icollocate(X), solution);
        Ns=subs(icollocate(W+X), solution);
        Ls=subs(icollocate(L), solution);
        Ps=subs(icollocate(P), solution);
        Us=subs(icollocate(U), solution); 


        
    else
     Result=-999999; % did not work

    end
    
        
        Bs= (aP.*Ls.^(bP)).*0.00001; %in grams
        Omegas=(Bs.*Ps).*0.001; %in kilograms      
        Ratios=Bs./((aN.*12.^bN)*.001); %Biomass ratio prawn to snail; I chose 8mm snail (medium); goes from 4mm (small) to 8mm (medium) to 12mm (large); see Hoover et al.
        alphaNs=aM.*log(Ratios); %attack rate
        Ths=((th.*Ratios).^(-1)*365); %Handling time estimated from Sokolow et al. (13; see Hoover et al)
               
        psiWs=(alphaNs.*epsilon.*Ws.^n)./(1+alphaNs.*epsilon.*Ths.*Ns.^n + alphaU.*Ths.*Us.^n);%per-capita prawn predation rate on snails 
        psiXs=(alphaNs.*epsilon.*Xs.^n)./(1+alphaNs.*epsilon.*Ths.*Ns.^n + alphaU.*Ths.*Us.^n);%per-capita prawn predation rate on snails 
        ks=max(kMax,alphaNs.*Ns)+alphaU.*Us.^(0.5); %Effective growth rate

        Avoided_HCosts= cI.*Is(1) - exp(-r.*ts).*cI.*Is;
        


        %Type of interpolation
        fitApprox = fittype('pchipinterp'); %pchipinterp
        %%% %%Interpolation of NPV
       % gCosts = fit(ts,exp(-r.*ts).*(cU_quad.*Us.^2 + cU.*Us),fitApprox); 
        gCosts = fit(ts,exp(-r.*ts).*(cU_quad.*Us.^2),fitApprox); 
        %gCosts = fit(ts,exp(-r.*ts).*(cU_quad.*Us),fitApprox); 
        Costs=integrate(gCosts,linspace(ts(1),ts(end),Nset+1),0);

        for i=1:Nset+1
        Profits=exp(-r.*ts).*(price.*Omegas)-cP.*Ps(1) - Costs(i);
        end
        

end



end