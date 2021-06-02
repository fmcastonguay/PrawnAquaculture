function [ts, Topt, Is, Ws, Xs, Ns, Ls, Ps, Bs, Omegas, Profits, psiWs, psiXs, psiNs, alphaNs, Ratios, Ths, ks, Results] = ...
    SchistoAquaculture_NoFeed(T,Nset,x0ic,beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,alphaY,Y,alphaU,ThY,ThU,eta,CASE,OBJ,GUESS)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE=1
if CASE==1 %% CASE== 1 : Aquaculture without feeding
    
toms t Tstar;
p=tomPhase('p', t, 0, Tstar, Nset,[],'gauss'); % Gauss has been working well (May 10, 2021); %'gauss' 'cheb' or 'fem1s' or 'fem1'
setPhase(p);
tomStates I W X L P

%%%Things that don't change with cases
state = [I; W; X; L; P];
combined = [I; W; X; L; P; Tstar]; 

%Terminal Conditions
cterm2=final({Tstar>=1/365}); %L>=140
%cterm=final({Tstar<=100});
cterm=final({}); %Tstar<=1

%%% Functions used in ODEs
N=X+W; %Total population of snails (N) equals the sum of the non-infected (X) and infected snails (W)
B=(aP.*L.^(bP)).*0.00001; % Mean Prawn body size (need to "*0.00001" to get in grams)
Omega=B.*P; % Biomass (i.e. mean prawn body size x number of prawns); needs to be in grams because "g" and "omega" are for grams
Ratio=B./((aN.*12.^bN)*.001);  %(need to "*0.001" to get in grams) %Biomass ratio prawn to snail; I chose 8mm snail (medium); goes from 4mm (small) to 8mm (medium) to 12mm (large); see Hoover et al.
alphaN=aM.*log(Ratio); %attack rate
Th=((th.*Ratio).^(-1)*365); %Handling time estimated from Sokolow et al. (13; see Hoover et al)

%epsilon=0.01;

psiW=(alphaN.*epsilon.*W.^n)./(1+alphaN.*epsilon.*Th.*N.^n);%per-capita prawn predation rate on infected snails 
psiX=(alphaN.*epsilon.*X.^n)./(1+alphaN.*epsilon.*Th.*N.^n);%per-capita prawn predation rate on healthy snails 


k=max(kMax,alphaN.*N)./(1+g.*Omega);%Modified Brody growth coefficient



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
     icollocate(0/365 <= Tstar)
 initial(I == x0ic(1))
 initial(W == x0ic(2))
 initial(X == x0ic(3))
 initial(L == x0ic(4))
 initial(P == x0ic(5))}; 




 %%% Initial guess
    if isempty(GUESS)==1
        %Creating initial guess
    x0_guess = { icollocate({I == x0ic(1);
    W==x0ic(2)
    X==x0ic(3)
    L==x0ic(4)
    P==x0ic(5)})};


    elseif isempty(GUESS)~=1
        

    x0_guess = { icollocate({I == GUESS(CASE).I;
        W==GUESS(CASE).W
        N==GUESS(CASE).X
        L==GUESS(CASE).L
        P==GUESS(CASE).P})};
        
    
    end

     
options=struct;
    if OBJ==1
        options.name='No Feed & One-Shot Horizon';
    elseif OBJ==2
        options.name='No Feed & Infinite Horizon';
    elseif OBJ==3
        options.name='Societal Optimum & One-Shot Horizon';
    elseif OBJ==4
        options.name='Societal Optimum & Infinite Horizon';
    elseif OBJ==5
        options.name='Health Optimum & One-Shot Horizon';
    elseif OBJ==6
        options.name='Health Optimum & Infinite Horizon';
    end
    

    options.solver = 'knitro'; % choose solver ('knitro' is other good one)

    
    Omega=Omega.*0.001; %in kilograms

  
if OBJ==1 %One-shot horizon
objective= - exp(-r.*Tstar).*price.*final(Omega) + cP.*initial(P);
elseif OBJ==2 %Infinite horizon
objective= - (1./(exp(r.*Tstar)-1)).*(price.*final(Omega) - cP.*initial(P)) + cP.*initial(P);
elseif OBJ==3 %Health Objective Single-Rotation
objective= cP.*initial(P) - cI.*initial(I) +  exp(-r.*Tstar).*cI.*final(I) ;
elseif OBJ==4 %Health Objective Infinite-Horizon
objective= - (1./(exp(r.*Tstar)-1)).*(cI.*(initial(I)-final(I))  - cP.*initial(P)) - cI.*initial(I) + cP.*initial(P) ;
elseif OBJ==5 %Benevolant Myopic Planner
objective= - exp(-r.*Tstar).*price.*final(Omega) + cP.*initial(P) - cI.*initial(I) +  exp(-r.*Tstar).*cI.*final(I) ;
elseif OBJ==6 %Benevolant Planner
objective= - (1./(exp(r.*Tstar)-1)).*(cI.*(initial(I)-final(I)) + price.*final(Omega) - cP.*initial(P)) - cI.*initial(I) + cP.*initial(P) ;
end



constr={ode,cbb,cterm,cterm2};

    [solution,result]= ezsolve(objective,constr,x0_guess,options);
    %This gives the solution of the first try as a guess for the second
    [solution,result]= ezsolve(objective,constr,solution,options);
    
    %%% Changing solver if exitflag~=0
    counter=0;
    while result.ExitFlag~=0 && counter<1  % limits the lopp size

        options.solver = 'snopt';
    [solution,result]= ezsolve(objective,constr,x0_guess,options);
        
        if result.ExitFlag~=0
            
    
        options.solver = 'npsol';
    [solution,result]= ezsolve(objective,constr,x0_guess,options);
        end
        
        if result.ExitFlag~=0
            
        options.solver = 'knitro';
    [solution,result]= ezsolve(objective,constr,x0_guess,options);
        end
        
        counter=counter+1;
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


        
    else
     Result=-999999; % did not work

    end
    
        Topt=ts(end)*365;
        
        Bs= (aP.*Ls.^(bP)).*0.00001; %in grams
        Omegas=(Bs.*Ps).*0.001; %in kilograms
        



        Ratios=Bs./((aN.*12.^bN)*.001); %Biomass ratio prawn to snail; I chose 8mm snail (medium); goes from 4mm (small) to 8mm (medium) to 12mm (large); see Hoover et al.
        alphaNs=aM.*log(Ratios); %attack rate
        Ths=((th.*Ratios).^(-1)*365); %Handling time estimated from Sokolow et al. (13; see Hoover et al)

        psiWs=(alphaNs.*epsilon.*Ws.^n)./(1+alphaNs.*epsilon.*Ths.*Ns.^n);%per-capita prawn predation rate on snails 
        psiXs=(alphaNs.*epsilon.*Xs.^n)./(1+alphaNs.*epsilon.*Ths.*Ns.^n);%per-capita prawn predation rate on snails 
        
        
        % - exp(-r.*Tstar).*price.*final(Omega) + cP.*initial(P) - cI.*initial(I) +  exp(-r.*Tstar).*cI.*final(I)
        psiNs=cI.*Is(1) - exp(-r.*ts).*cI.*Is;
         

        ks=max(kMax, alphaNs.*Ns);

        %Type of interpolation
        fitApprox = fittype('pchipinterp'); %pchipinterp

        %Interpolation of NPV
        gRevenues = fit(ts,exp(-r.*ts).*(price.*Omegas),fitApprox); 
        Revenues=integrate(gRevenues,linspace(ts(1),ts(end),Nset+2),0);

        %Profits=Revenues-cP.*Ps(1);
        Profits=exp(-r.*ts).*(price.*Omegas)-cP.*Ps(1);

end



end