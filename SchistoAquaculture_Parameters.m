
function [beta,lambda,gamma,Linf,delta,r,f,cP,cI,price,kMax,g,d,aP,bP,muP,omega,aN,bN,aM,th,K,epsilon,n,cU,alphaY,Y,alphaU,ThY,ThU,eta]= SchistoAquaculture_Parameters()


%%%Parameters

%Village/Water parameters parameters
params.Area=1000; %1000 m2 of water
params.K=50*params.Area; % Carrying capacity of snails in the 1000 m2 of water
params.H=5000; % Human population

Area=params.Area; K=params.K; H=params.H;

%Aquaculture/Prawn Parameters
params.aP= 0.07329961; %0.00244; %  % % Value used in Hoover, but they cite this ; % Allometric parameters for prawns length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii males)
params.bP=3.5502; % (3.4617 - 3.6386)  Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, males)
params.g=3.5e-06; % Density dependent growth reduction (could be adjusted)
params.muP=2.21; %2.21/12;  %(DEPENDS ON TIME) Natural dyearly prawn mortality rate (M. volenhovenii males) from Nwosu & Wolfi 2006 %(0.0061.*365); %daily prawn mortality rate
params.d=-0.382; %size-dependent mortality scaling coefficient
params.omega=5.5e-09;%5.5e-09;% Density dependent mortality factor
params.Linf=213.63; % Maximum length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii males)
params.kMax=  3.7960;% 4.452;% 3.19  % (DEPENDS ON TIME) Growth rate (mm/day) of M. rosenbergii; 3.19/365; % Growth rate (mm/year) of M. vollenhovenii (0.0104*365) (mm/day) ; % Growth rate of prawns (Hoover; macrobrachium rosenbergii); 0.371*12 for M. Rosenbergii



Linf=params.Linf; kMax=params.kMax; g=params.g;  aP=params.aP; bP=params.bP; muP=params.muP; d=params.d; omega=params.omega; 

%Predation Parameters
params.aN=0.187178454; % Allometric parameters for snails' length-weight relationship
params.bN=2.536764792; % Idem as above;
params.aM = 0.9050; %Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
params.th = 0.38561;     % Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
params.n=2; %Exponent of Hollingtype III functional response
params.epsilon=0.1; %(epsilon=Prawn predation attack rate penalty associated withsearching for prey in wildlife rather than laboratory conditions) range=[0.01,1]

aN=params.aN; bN=params.bN; aM=params.aM; th=params.th; n=params.n; epsilon=params.epsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Feed and non-snail prey
%%% This is to test the snail-dependent somatic growth
B_min=(aP.*40.^(bP)).*0.00001;
B_max=(aP.*Linf.^(bP)).*0.00001;
N_min=0.0042;

Ratio_min=B_min./((aN.*12.^bN)*.001);
Ratio_max=B_max./((aN.*12.^bN)*.001);
alphaN_min=aM.*log(Ratio_min);

Y=1;
alphaY=kMax;% - alphaN_min*N_min;
alphaU=alphaY./2;

ThY=((th.*Ratio_min).^(-1)*365);
ThU=((th.*Ratio_min).^(-1)*365);
eta=1; %500





%Snail parameters 
params.f=0.16*Area*365/K;% (DEPENDS ON TIME)(0.16*365*200)/10000; %Snail Growth Rate; (0.16 per day*365 days*200 m^2)/10000 (carrying capacity) (Sokolow et al. reduced transmission)
params.delta=6/K;%(DEPENDS ON TIME) %death rate of the disease (i.e. snails, 2months) in the environment *** See Kariuki, H. Curtis, et al. "Divergent effects of Schistosoma haematobium exposure on intermediate-host snail species Bulinus nasutus and Bulinus globosus from coastal Kenya." The American journal of tropical medicine and hygiene 96.4 (2017): 850-855. for snail life span; see Sokolow et al. (2015) for prawn density per site (i.e. 250

f=params.f; delta=params.delta; 

% Epidemiologivcal Parameters
params.beta = 1/200*20*2.077893*2*5/4/100;%/100;%1/200;  %contact rate ; daily infection probability from snail to man ; Sokolow et al.(2015)
params.lambda=4e-05*3.5*2.5/2/20/2/5*4*100;%*100;%1/555; %shedding rate ; per capita snail infection probability ; Sokolow et al. *375 (gives an R0 of 3.5)
params.gamma=(1/(3.3))/70;%/H; %(DEPENDS ON TIME) %death rate of the disease in human; i.e. adult worm lifespan in human host ; Sokolow et al.(2015) (see worm lifespan (3.3 years) and mean adult worm (70 worms) per human host)

params.q=.8;   %Efficiency of MDA control

beta=params.beta; lambda=params.lambda; gamma=params.gamma; q=params.q;

%Economic parameters
params.r=.07; % (DEPENDS ON TIME) Discount rate
params.cI=520670; % cost of infection to adults \\\\ (from Lo et al. 2016) 
params.cP=0.1; % Cost of juvenile prawns (Hoover et al)
params.cU= 109.601*0; %0.729;  % %%250   % Cost of feeding based off of Dasgupta and Tidwell 2003: 254$/metric ton, need 4315 kg/hectare/year and we have 0.1 hectare

params.c31=H.*.105; %(assuming 5000 people) cost of MDA treatment to children  (from Lo et al. 2015, 2016)
params.F=1000*0; % Fixed Cost
params.price=12;%/1e+05; % Price of prawns (12$/kg); Omega is given in Micrograms Hoover et al.; estimates from Dasgupta and Tidwell

r=params.r; cI=params.cI; cP=params.cP; cU=params.cU ; c31=params.c31 ; F=params.F; price=params.price;


%Parameters for which I have no data: cost of feed

%R naught
if true
efF= [0 beta ; 0 0];
V= [-gamma  0; lambda -delta ];
FV=-efF*inv(V);
R0=eigs(FV,1,'lr');

% Rnaught = (lambda*beta)/(gamma*delta)
end


end





