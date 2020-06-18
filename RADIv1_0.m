%% 1-D Reaction-Advection-Diffusion-Irrigation (RADI) Diagenetic Sediment Module
%% Source code by O. Sulpis and M. Wilhelmus
%% Optimised for MATLAB by M.P. Humphreys [v20, March 2020]
%% Uses: CO2SYS and calc_pco2

disp("RADIv20 is running the following experiment:")
disp(Station)
if rerun==1
    disp("it is a rerun: you can stop the run at any time to visualize the evolution of the system")
    disp("years after start of simulation:")
elseif rerun==0
    disp("it is not a rerun: you can stop the run at any time to visualize the evolution of the system")
    disp("when you stop, set rerun to '1' to continue the analysis")
    disp("years after start of simulation")
else
    disp('initial conditions loaded')
    disp("years after start of simulation")
end

tStart = tic;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% carbonate system initialization this is used only at the first time step to initialize the calc_pco2 program
CO2SYS_data = CO2SYS(dtalkw*1e6/rho_sw,dtCO2w*1e6/rho_sw,1,2,S,T,T,P*10,P*10,dSiw*1e6/rho_sw,dtPO4w*1e6/rho_sw,1,10,1);
k1(1,1:ndepths) = CO2SYS_data(1,67);           %carbonic acid first dissociation constant
k2(1,1:ndepths) = CO2SYS_data(1,68);           %carbonic acid second dissociation constant
k1p(1,1:ndepths) = CO2SYS_data(1,75);         %phosphate constant 1
k2p(1,1:ndepths) = CO2SYS_data(1,76);         %phosphate constant 2
k3p(1,1:ndepths) = CO2SYS_data(1,77);         %phosphate constant 3
kb(1,1:ndepths) = CO2SYS_data(1,72);            %boron constant 
kw(1,1:ndepths) = CO2SYS_data(1,71);           %water dissociation constants
ksi(1,1:ndepths) = CO2SYS_data(1,78);           %silica constants
bt(1,1:ndepths) = CO2SYS_data(1,79);             %[umol/kg] total boron 
omegaC = CO2SYS_data(1,30);                         %calcite saturation state
omegaA = CO2SYS_data(1,31);                          %aragonite saturation state
co3 = CO2SYS_data(1,22) .* 10^-6;                   %[mol/kg] CO3 
hco3 = CO2SYS_data(1,21) .* 10^-6;                 %[mol/kg] HCO3
ph = CO2SYS_data(1,37) ;                                   %pH on the total scale
Ca_ini = dCaw ./ rho_sw;                                 %[mol/kg] Ca concentration 
fg(1,1:ndepths)=dtalkw./rho_sw-hco3-2*co3; %sum of all alkalinity species that are not carbon
kspc = (co3 .* Ca_ini) ./ omegaC;                %[mol2/kg2] calcite in situ solubility
kspa = (co3 .* Ca_ini) ./ omegaA;                %[mol2/kg2] aragonite in situ solubility
ff(1,1:ndepths) = 1;                                        %random parameter needed for calc_pco2
H(1,1:ndepths) = 10^-ph;                             %[mol/kg] H concentration first guess
clear co3 hco3 ph omegaC omegaA Ca_ini CO2SYS_data   
sit = 120 * 10^-6;                                              %[mol/kg] convert silica concentration
bt = bt .* 10^-6;                                                %[mol/kg] convert boron concentration

%% temperature dependent "free solution" diffusion coefficients
D_dO2=0.034862+0.001409*T; %[m2/a] oxygen diffusion coefficient from Li and Gregory
D_dtalk=0.015169+0.000793*T;         %[m2/a] approximted to bicarbonate diffusion coefficient from Hulse et al (2018)
D_dtCO2=0.015169+0.000793*T;         %[m2/a] approximted to bicarbonate diffusion coefficient from Hulse et al (2018)
D_dtNO3=0.030842+0.001226*T;        %[m2/a] nitrate diffusion coefficient from Li and Gregory (1974)
D_dtSO4=0.015768+0.000788*T;        %[m2/a] sulfate diffusion coefficient from Li and Gregory (1974)
D_dtPO4=0.011291+0.000559*T;        %[m2/a] phosphate diffusion coefficient from Li and Gregory (1974)
D_dtNH4=0.030905+0.001226*T;        %[m2/a] ammonium diffusion coefficient from Li and Gregory (1974)
D_dtH2S=0.030748+0.000964*T;        %[m2/a] hydrogen sulfide diffusion coefficient from the UNISENSE table by Ramsing and Gundersen
D_dMn=0.0086+0.001525*T;           %[m2/a] manganese diffusion coefficient from Li and Gregory (1974)
D_dFe=0.0108+0.001478*T;           %[m2/a] iron diffusion coefficient from Li and Gregory (1974)
D_dCa=0.0107+0.001677*T;         %[m2/a] calcium diffusion coefficient from Li and Gregory (1974)

%% bioturbation (for solids)
D_bio_0=1e-4*0.0232*(Foc*1e2)^0.85; %[m2/a] surf bioturb coeff, Archer et al (2002)
lambda_b = 0.08;
D_bio=D_bio_0*exp(-(depths./lambda_b).^2).*((dO2w/1e-3)/((dO2w/1e-3)+20)); %[m2/a] bioturb coeff, Archer et al (2002)

%% irrigation (for solutes)
alpha_0=11*(atan((5*Foc*1e2-400)/400)/pi+0.5)-0.9...
    +20*((dO2w/1e-3)/((dO2w/1e-3)+10))*exp(-(dO2w/1e-3)/10)*Foc*1e2/(Foc*1e2+30);    %[/a] from Archer et al (2002)
lambda_i=0.05;
alpha=alpha_0.*exp(-(depths/lambda_i).^2);                                                                                   %[/a] Archer et al (2002) the depth of 5 cm was changed

%% depth-dependent porosity and diffusion coefficient loss
% % Use differences - values should then be divided by z_res?
% delta_phi = [0 diff(phi)]; % depth-dependent porosity loss
% delta_phiS = [0 diff(phiS)]; % depth-dependent solid fraction gain
% delta_tort2 = [0 diff(tort.^2)]; % depth-dependent tortuosity gain
% delta_D_bio_i = [0 diff(D_bio)]; % [m/a]
% Use derivative equations instead! all checked vs finite differences
delta_phi = -phiBeta.*(phi0 - phiInf).*exp(-phiBeta*depths);
% delta_phi(1) = 0; % don't do this
delta_phiS = -delta_phi;
delta_tort2 = -2*delta_phi./phi; % not used in Julia
delta_D_bio = -2*depths.*D_bio/lambda_b^2; % not used in Julia

% biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis (1977)
Peh=w.*z_res./(2*D_bio);      %one half the cell Peclet number (Eq. 97 in Boudreau 1996)
% when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
sigma=1./tanh(Peh)-1./(Peh);  %Eq. 96 in Boudreau 1996

%% organic matter degradation parameters

 KdO2=0.003; %[mol/m3] Monod constant from Soetaert et al. 1996 (GCA)
 KindO2=0.01; %[mol/m3] Monod inhibition constant from Soetaert et al. 1996 (GCA)
 KdtNO3=0.03; %[mol/m3] Monod constant from Soetaert et al. 1996 (GCA)
 KindtNO3=0.005; %[mol/m3] Monod inhibition constant from Soetaert et al. 1996 (GCA)
 KpMnO2=42.4; %[mol/m3] Monod constant from Van Cappellen and Wang 1996
 KinpMnO2=KpMnO2; %[mol/m3] Monod inhibition constant from Van Cappellen and Wang 1996
 KpFeOH3=265; %[mol/m3] Monod constant from Van Cappellen and Wang 1996 
 KinpFeOH3=KpFeOH3; %[mol/m3] Monod inhibition constant from Van Cappellen and Wang 1996
 KdtSO4=1.6; %[mol/m3] Monod constant from Van Cappellen and Wang 1996
 KindtSO4=KdtSO4; %[mol/m3] Monod inhibition constant from Van Cappellen and Wang 1996

kslow_0=1e-4 * (Foc*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
lambda_slow=1;     %[m] tuned parameter
kslow=kslow_0*exp(-depths./lambda_slow);    %[/a] from Archer et al (2002)

kfast_0=1e-2 * (Foc*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
lambda_fast=0.03;     %[m] tuned parameter
kfast=kfast_0*exp(-depths./lambda_fast);    %[/a] from Archer et al (2002)

%% redox reaction first order constants
kMnox=1e6; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
kFeox=1e6; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
kNHox=1e4; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
kSox=3e5; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% R.A.D.I. main loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerun == 0 %concentrations set to zero for solids, bottom water values for not-too-sensitive solutes
    dO2(1,1:ndepths)=0;        
    dtalk(1,1:ndepths)=0;            
    dtCO2(1,1:ndepths)=0;            
    dtNO3(1,1:ndepths)=0;      
    dtSO4(1,1:ndepths)=0;      
    dtPO4(1,1:ndepths)=0;      
    dtNH4(1,1:ndepths)=0;      
    dtH2S(1,1:ndepths)=0;      
    dFe(1,1:ndepths)=0;      
    dMn(1,1:ndepths)=0;      
    proc(1,1:ndepths)=0;
    psoc(1,1:ndepths)=0;             
    pfoc(1,1:ndepths)=0;        
    pFeOH3(1,1:ndepths)=0;    
    pMnO2(1,1:ndepths)=0;    
    dCa(1,1:ndepths)=0;            
    
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime:t_length;  %we will only keep the variables every year
    plot_number(1)=1;
    
elseif rerun==1 %if it is a rerun, initial conditions are concentrations from last time step
    dO2=dO2f(idx-1,:);            %[mol/m3]
    dtalk=dtalkf(idx-1,:);            %[mol/m3]
    dtCO2=dtCO2f(idx-1,:);      %[mol/m3]
    dtNO3=dtNO3f(idx-1,:);            %[mol/m3]
    dtSO4=dtSO4f(idx-1,:);            %[mol/m3]
    dtPO4=dtPO4f(idx-1,:);            %[mol/m3]
    dtNH4=dtNH4f(idx-1,:);            %[mol/m3]
    dtH2S=dtH2Sf(idx-1,:);            %[mol/m3]
    dFe=dFef(idx-1,:);            %[mol/m3]
    dMn=dMnf(idx-1,:);            %[mol/m3]
    proc=proc(idx-1,:);                %[mol/m3]
    psoc=psoc(idx-1,:);                %[mol/m3]
    pfoc=pfoc(idx-1,:);                %[mol/m3]
    pFeOH3=pFeOH3f(idx-1,:);                %[mol/m3]
    pMnO2=pMnO2f(idx-1,:);                %[mol/m3]
    dCa=dCaf(idx-1,:);            %[mol/m3]

    plot_number=0:t_length/stoptime:t_length;  %we will only keep the variables every year
    i=plot_number(idx-1);
    
else
    % initial condition for solutes: bottom-water value
    dO2=dO2ic;                %[mol/m3]
    dtalk=dtalkic;                %[mol/m3]
    dtCO2=dtCO2ic;                %[mol/m3]
    dtNO3=dtNO3ic;            %[mol/m3]
    dtSO4=dtSO4ic;            %[mol/m3]
    dtPO4=dtPO4ic;            %[mol/m3]
    dtNH4=dtNH4ic;            %[mol/m3]
    dtH2S=dtH2Sic;            %[mol/m3]
    dFe=dFeic;            %[mol/m3]
    dMn=dMnic;            %[mol/m3]
    dCa=dCaic;            %[mol/m3]
    proc=procic;
    psoc=psocic;
    pfoc=pfocic;
    pFeOH3=pFeOH3ic;
    pMnO2=pMnO2icic;
    
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime:t_length;  %we will only keep the variables every year
    plot_number(1)=1;    
end

%% short-cut transport variables
APPW=w - delta_D_bio - delta_phiS.* D_bio./ phiS;
DFF=(tort2.* delta_phi./ phi - delta_tort2)./ (tort2.^2);
DBF=phiS.*D_bio+D_bio.*(-delta_phi);
TR=(2*z_res.* (tort.^2)./ dbl);

%% Prepare for timestep calculations // MPH [v20]

% Set indices for depth-varying reactions
j = 2:ndepths-1;
jp1 = j + 1;
jm1 = j - 1;
z_res2 = z_res.^2;

% Subsample variables that do not change from timestep to timestep
alpha_j = alpha(j);
z_res_j = z_res(j);
z_res2_j = z_res2(j);
DFF_j = DFF(j);
tort2_j = tort2(j);
u_j = u(j);
APPW_j = APPW(j);
sigma_j = sigma(j);
D_bio_j = D_bio(j);
phi_j = phi(j);
phiS_j = phiS(j);

% % Precalculations - these don't actually speed things up!
% uDFF_TA = u_j - D_TA*DFF_j;
% Dtort2_TA = D_TA./tort2_j;

%% Begin timestep loop
% Start with some O2 and OC
dO2(:) = dO2w;
dtalk(:) = dtalkw;
dtCO2(:) = dtCO2w;
dtNO3(:)=dtNO3w;            %[mol/m3]
dtSO4(:)=dtSO4w;            %[mol/m3]
dtPO4(:)=dtPO4w;            %[mol/m3]
dtNH4(:)=dtNH4w;            %[mol/m3]
dtH2S(:)=dtH2Sw;            %[mol/m3]
dFe(:)=dFew;            %[mol/m3]
dMn(:)=dMnw;            %[mol/m3]
dCa(:)=dCaw;
proc(:) = 3e4;
psoc(:) = 3e3;
pfoc(:) = 3e2;
pFeOH3(:)= 0;
pMnO2(:) = 0;

% Preallocate saving arrays
%dO2f = NaN(ndepths, t_length);
%dtalkf = NaN(ndepths, t_length);
%dtCO2f = NaN(ndepths, t_length);
%dtNO3f = NaN(ndepths, t_length);
%dtSO4f = NaN(ndepths, t_length);
%dtPO4f = NaN(ndepths, t_length);
%dtNH4f = NaN(ndepths, t_length);
%dtH2Sf = NaN(ndepths, t_length);
%dFef = NaN(ndepths, t_length);
%dMnf = NaN(ndepths, t_length);
%dCaf = NaN(ndepths, t_length);
%procf = NaN(ndepths, t_length);
%psocf = NaN(ndepths, t_length);
%pfocf = NaN(ndepths, t_length);
%pFeOH3f = NaN(ndepths, t_length);
%pMnO2f = NaN(ndepths, t_length);
dO2f = NaN(ndepths, stoptime+1);
dtalkf = NaN(ndepths, stoptime+1);
dtCO2f = NaN(ndepths, stoptime+1);
dtNO3f = NaN(ndepths, stoptime+1);
dtSO4f = NaN(ndepths, stoptime+1);
dtPO4f = NaN(ndepths, stoptime+1);
dtNH4f = NaN(ndepths, stoptime+1);
dtH2Sf = NaN(ndepths, stoptime+1);
dFef = NaN(ndepths, stoptime+1);
dMnf = NaN(ndepths, stoptime+1);
dCaf = NaN(ndepths, stoptime+1);
procf = NaN(ndepths, stoptime+1);
psocf = NaN(ndepths, stoptime+1);
pfocf = NaN(ndepths, stoptime+1);
pFeOH3f = NaN(ndepths, stoptime+1);
pMnO2f = NaN(ndepths, stoptime+1);

for i=i:t_length-1

%     disp(i)
    
    %F_O2i=D_O2*phi(1)*(O2(:,1)-O2w)./5e-3;
    %F_DICi=D_DIC*phi(1)*(DIC(:,1)-DICw)./5e-3;
    %F_TAi=D_TA*phi(1)*(TA(:,1)-TAw)./5e-3;
    
       
    %% Organic matter respiration pathways
    fdO2=dO2./(KdO2+dO2);                   %from the code of Couture et al. (EST 2010), following Boudreau (1996)
    fdtNO3=dtNO3./(KdtNO3+dtNO3).*(KindO2./(KindO2+dO2));
    fpMnO2=pMnO2./(pMnO2+KpMnO2).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fpFeOH3=pFeOH3./(pFeOH3+KpFeOH3).*(KinpMnO2./(pMnO2+KinpMnO2)).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fdtSO4=dtSO4./(dtSO4+KdtSO4).*(KinpFeOH3./(pFeOH3+KinpFeOH3)).*(KinpMnO2./(pMnO2+KinpMnO2)).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fdCH4=(KindtSO4./(dtSO4+KindtSO4)).*(KinpFeOH3./(pFeOH3+KinpFeOH3)).*(KinpMnO2./(pMnO2+KinpMnO2)).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fox=fdO2+fdtNO3+fpMnO2+fpFeOH3+fdtSO4+fdCH4;
    
    %% Redox reaction rates
    Rs_o2 = psoc.*kslow.*fdO2; %degradation of slow-decay organic matter by oxygen
    Rf_o2 = pfoc.*kfast.*fdO2; %degradation of fast-decay organic matter by oxygen
    Rs_no3 = psoc.*kslow.*fdtNO3; %...by nitrate
    Rf_no3 = pfoc.*kfast.*fdtNO3;
    Rs_mno2 = psoc.*kslow.*fpMnO2; %... by manganese oxide
    Rf_mno2 = pfoc.*kfast.*fpMnO2;
    Rs_feoh3 = psoc.*kslow.*fpFeOH3; %by iron oxide
    Rf_feoh3 = pfoc.*kfast.*fpFeOH3;
    Rs_so4 = psoc.*kslow.*fdtSO4; %... by sulfate
    Rf_so4 = pfoc.*kfast.*fdtSO4;
    Rs_ch4 = psoc.*kslow.*fdCH4; %... by itself
    Rf_ch4 = pfoc.*kfast.*fdCH4;
    Rs_tot = psoc.*kslow.*fox; % total degradation rate of slow-decay organic matter
    Rf_tot = pfoc.*kfast.*fox;  % total degradation rate of fast-decay organic matter
    RFeox = kFeox.*dFe.*dO2; % oxidation of dissolved iron
    RMnox = kMnox.*dMn.*dO2; % oxidation of dissolved manganese 
    RSox = kSox.*dtH2S.*dO2; % oxidation of hydrogen sulfide
    RNHox = kNHox.*dtNH4.*dO2; % oxidation of ammonia
    
    %% calc_pco2: time efficient carbonate system solver    
    DIC_molPerKg = dtCO2 ./ rho_sw;     %convert DIC to [mol/kg]
    TA_molPerKg = dtalk ./ rho_sw;         %convert TA to [mol/kg]
    PO4_molPerKg = dtPO4 ./ rho_sw;   %convert PO4 to [mol/kg]
    
    [H] = calc_pCO2(DIC_molPerKg,PO4_molPerKg,sit,bt,TA_molPerKg,ff,k1,k2,k1p,k2p,k3p,kb,kw,ksi,fg,H);    %[mol/kg] H concentration
    
    co3_molPerKg = (DIC_molPerKg .* k1 .* k2)./ (H.^2 + (k1 .* H) + (k1 .* k2)); %[mol/kg] CO3 concentration
    co3 = co3_molPerKg .* rho_sw; %[mol m^-3] CO3 concentraiton
    
    %% CaCO3 reactions
    OmegaC = dCa.*co3./ (kspc.* rho_sw.^2); %[no unit] calcite saturation state
    OmegaA = dCa.*co3./ (kspa.* rho_sw.^2); %[no unit] aragonite saturation state
    
    %% Calculate all reactions (14 species, units: [mol/m3/a])
    % This section ~2x faster by not putting all the reactions into a
    % single matrix but keeping as separate vectors // MPH
    TotR_dO2 = - phiS./phi.*(Rs_o2 + Rf_o2) - 0.25.*RFeox - 0.5.*RMnox - 2.*RSox - 2.*RNHox;
    TotR_dtalk = + phiS./ phi.*(Rs_o2.*(RN./RC - RP./RC) + Rf_o2.*(RN./RC - RP./RC) + Rs_no3.*(0.8+RN./RC - RP./RC)...
        + Rf_no3.*(0.8+RN./RC - RP./RC) + Rs_mno2.*(4+RN./RC - RP./RC) + Rf_mno2.*(4+RN./RC - RP./RC)...
        + Rs_feoh3.*(8+RN./RC - RP./RC) + Rf_feoh3.*(0.8+RN./RC - RP./RC) + Rs_so4.*(1+RN./RC - RP./RC)...
        + Rf_so4.*(1+RN./RC - RP./RC) + Rs_ch4.*(RN./RC - RP./RC) + Rf_ch4.*(RN./RC - RP./RC)) - 2.* RFeox...
            - 2.* RMnox - 2.* RSox - 2.* RNHox;
    TotR_dtCO2 = phiS./phi.*(Rs_o2 + Rf_o2 + Rs_no3 + Rf_no3 + Rs_mno2 + Rf_mno2 + Rs_feoh3 + Rf_feoh3...
        + Rs_so4 + Rf_so4 + Rs_ch4.*0.5 + Rf_ch4.*0.5);
    TotR_dtNO3 = - phiS./phi.*0.8.*(Rs_no3 + Rf_no3) + RNHox; 
    TotR_dtSO4 = - phiS./phi.*0.5.*(Rs_so4 + Rf_so4) + RSox; 
    TotR_dtPO4 = phiS./phi.*(RP./RC).*(Rs_tot + Rf_tot);
    TotR_dtNH4 = phiS./phi.*(RN./RC).*(Rs_tot + Rf_tot) - RNHox;
    TotR_dtH2S = phiS./phi.*0.5.*(Rs_so4 + Rf_so4) - RSox;    
    TotR_dFe = phiS./phi.*4.*(Rs_feoh3 + Rf_feoh3) - RFeox;
    TotR_dMn = phiS./phi.*2.*(Rs_mno2 + Rf_mno2) - RMnox;
    TotR_psoc = - Rs_tot;
    TotR_pfoc = - Rf_tot;
    TotR_pFeOH3 = - 4.*(Rs_feoh3 + Rf_feoh3) + phi./phiS.*RFeox;
    TotR_pMnO2 = - 2.*(Rs_mno2 + Rf_mno2) + phi./phiS.*RMnox;
    
    %% top boundary condition: prescribed solid fluxes and diffusive boundary layer control on solutes
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    
    dO2_1 = dO2(1) + interval * ( D_dO2 / tort2(1) * (2*dO2(2) - 2*dO2(1) + TR(1) * (dO2w - dO2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dO2.*DFF(1)) * -1 * TR(1) * ( dO2w - dO2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dO2w - dO2(1) ) ... %irrigation
        + TotR_dO2(1)); %reaction
    
     dtalk_1 = dtalk(1) + interval * ( D_dtalk / tort2(1) * (2*dtalk(2) - 2*dtalk(1) + TR(1) * (dtalkw - dtalk(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtalk.*DFF(1)) * -1 * TR(1) * ( dtalkw - dtalk(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtalkw - dtalk(1) ) ... %irrigation
        + TotR_dtalk(1)); %reaction
    
    dtCO2_1 = dtCO2(1) + interval * ( D_dtCO2 / tort2(1) * (2*dtCO2(2) - 2*dtCO2(1) + TR(1) * (dtCO2w - dtCO2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtCO2.*DFF(1)) * -1 * TR(1) * ( dtCO2w - dtCO2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtCO2w - dtCO2(1) ) ... %irrigation
        + TotR_dtCO2(1)); %reaction
    
    dtNO3_1 = dtNO3(1) + interval * ( D_dtNO3 / tort2(1) * (2*dtNO3(2) - 2*dtNO3(1) + TR(1) * (dtNO3w - dtNO3(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtNO3.*DFF(1)) * -1 * TR(1) * ( dtNO3w - dtNO3(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtNO3w - dtNO3(1) ) ... %irrigation
        + TotR_dtNO3(1)); %reaction
    
    dtSO4_1 = dtSO4(1) + interval * ( D_dtSO4 / tort2(1) * (2*dtSO4(2) - 2*dtSO4(1) + TR(1) * (dtSO4w - dtSO4(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtSO4.*DFF(1)) * -1 * TR(1) * ( dtSO4w - dtSO4(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtSO4w - dtSO4(1) ) ... %irrigation
        + TotR_dtSO4(1)); %reaction
    
    dtPO4_1 = dtPO4(1) + interval * ( D_dtPO4 / tort2(1) * (2*dtPO4(2) - 2*dtPO4(1) + TR(1) * (dtPO4w - dtPO4(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtPO4.*DFF(1)) * -1 * TR(1) * ( dtPO4w - dtPO4(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtPO4w - dtPO4(1) ) ... %irrigation
        + TotR_dtPO4(1)); %reaction
    
    dtNH4_1 = dtNH4(1) + interval * ( D_dtNH4 / tort2(1) * (2*dtNH4(2) - 2*dtNH4(1) + TR(1) * (dtNH4w - dtNH4(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtNH4.*DFF(1)) * -1 * TR(1) * ( dtNH4w - dtNH4(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtNH4w - dtNH4(1) ) ... %irrigation
        + TotR_dtNH4(1)); %reaction
    
    dtH2S_1 = dtH2S(1) + interval * ( D_dtH2S / tort2(1) * (2*dtH2S(2) - 2*dtH2S(1) + TR(1) * (dtH2Sw - dtH2S(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtH2S.*DFF(1)) * -1 * TR(1) * ( dtH2Sw - dtH2S(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtH2Sw - dtH2S(1) ) ... %irrigation
        + TotR_dtH2S(1)); %reaction
    
    dFe_1 = dFe(1) + interval * ( D_dFe / tort2(1) * (2*dFe(2) - 2*dFe(1) + TR(1) * (dFew - dFe(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dFe.*DFF(1)) * -1 * TR(1) * ( dFew - dFe(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dFew - dFe(1) ) ... %irrigation
        + TotR_dFe(1)); %reaction
    
    dMn_1 = dMn(1) + interval * ( D_dMn / tort2(1) * (2*dMn(2) - 2*dMn(1) + TR(1) * (dMnw - dMn(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dMn.*DFF(1)) * -1 * TR(1) * ( dMnw - dMn(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dMnw - dMn(1) ) ... %irrigation
        + TotR_dMn(1)); %reaction

    dCa_1 = dCa(1) + interval * ( D_dCa / tort2(1) * (2*dCa(2) - 2*dCa(1) + TR(1) * (dCaw - dCa(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dCa.*DFF(1)) * -1 * TR(1) * ( dCaw - dCa(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dCaw - dCa(1) ) ... %irrigation
        + 0); %reaction
    
    proc_1 = proc(1) + interval * ( D_bio(1) * ( 2 * proc(2) - 2 * proc(1) +... %diffusion
        2 * z_res(1) * (Froc - phiS(1) * w(1) * proc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*proc(2) + 2*sigma(1)*proc(1) - ... %advection
        (1 + sigma(1))*(proc(2)+2*z_res(1)/D_bio(1)*(Froc/phiS(1)-w(1)*proc(1))))/(2*z_res(1))); ... %advection
            
    psoc_1 = psoc(1) + interval * ( D_bio(1) * ( 2 * psoc(2) - 2 * psoc(1) +... %diffusion
        2 * z_res(1) * (Fsoc - phiS(1) * w(1) * psoc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*psoc(2) + 2*sigma(1)*psoc(1) - ... %advection
        (1 + sigma(1))*(psoc(2)+2*z_res(1)/D_bio(1)*(Fsoc/phiS(1)-w(1)*psoc(1))))/(2*z_res(1)) ... %advection
        +TotR_psoc(1)); %reaction
    
    pfoc_1 = pfoc(1) + interval * ( D_bio(1) * ( 2 * pfoc(2) - 2 * pfoc(1) +... %diffusion
        2 * z_res(1) * (Ffoc - phiS(1) * w(1) * pfoc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pfoc(2) + 2*sigma(1)*pfoc(1) - ... %advection
        (1 + sigma(1))*(pfoc(2)+2*z_res(1)/D_bio(1)*(Ffoc/phiS(1)-w(1)*pfoc(1))))/(2*z_res(1)) ... %advection
        +TotR_pfoc(1)); %reaction
    
     pFeOH3_1 = pFeOH3(1) + interval * ( D_bio(1) * ( 2 * pFeOH3(2) - 2 * pFeOH3(1) +... %diffusion
        2 * z_res(1) * (FFeOH3 - phiS(1) * w(1) * pFeOH3(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pFeOH3(2) + 2*sigma(1)*pFeOH3(1) - ... %advection
        (1 + sigma(1))*(pFeOH3(2)+2*z_res(1)/D_bio(1)*(FFeOH3/phiS(1)-w(1)*pFeOH3(1))))/(2*z_res(1)) ... %advection
        +TotR_pFeOH3(1)); %reaction

    pMnO2_1 = pMnO2(1) + interval * ( D_bio(1) * ( 2 * pMnO2(2) - 2 * pMnO2(1) +... %diffusion
        2 * z_res(1) * (FMnO2 - phiS(1) * w(1) * pMnO2(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pMnO2(2) + 2*sigma(1)*pMnO2(1) - ... %advection
        (1 + sigma(1))*(pMnO2(2)+2*z_res(1)/D_bio(1)*(FMnO2/phiS(1)-w(1)*pMnO2(1))))/(2*z_res(1)) ... %advection
        +TotR_pMnO2(1)); %reaction
 
    %% bottom boundary condition: gradients disappear
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    dO2_z = dO2(ndepths) + interval * (D_dO2 / tort2(ndepths) * 2 * ((dO2(ndepths-1) - dO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dO2w - dO2(ndepths)) ... %irrigation
        +TotR_dO2(ndepths));

    dtalk_z = dtalk(ndepths) + interval * (D_dtalk / tort2(ndepths) * 2 * ((dtalk(ndepths-1) - dtalk(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtalkw - dtalk(ndepths)) ... %irrigation
        +TotR_dtalk(ndepths));
    
    dtCO2_z = dtCO2(ndepths) + interval * (D_dtCO2 / tort2(ndepths) * 2 * ((dtCO2(ndepths-1) - dtCO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtCO2w - dtCO2(ndepths)) ... %irrigation
        +TotR_dtCO2(ndepths));
    
    dtNO3_z = dtNO3(ndepths) + interval * (D_dtNO3 / tort2(ndepths) * 2 * ((dtNO3(ndepths-1) - dtNO3(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtNO3w - dtNO3(ndepths)) ... %irrigation
        +TotR_dtNO3(ndepths));
    
    dtSO4_z = dtSO4(ndepths) + interval * (D_dtSO4 / tort2(ndepths) * 2 * ((dtSO4(ndepths-1) - dtSO4(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtSO4w - dtSO4(ndepths)) ... %irrigation
        +TotR_dtSO4(ndepths));
    
    dtPO4_z = dtPO4(ndepths) + interval * (D_dtPO4 / tort2(ndepths) * 2 * ((dtPO4(ndepths-1) - dtPO4(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtPO4w - dtPO4(ndepths)) ... %irrigation
        +TotR_dtPO4(ndepths));
    
    dtNH4_z = dtNH4(ndepths) + interval * (D_dtNH4 / tort2(ndepths) * 2 * ((dtNH4(ndepths-1) - dtNH4(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtNH4w - dtNH4(ndepths)) ... %irrigation
        +TotR_dtNH4(ndepths));
    
    dtH2S_z = dtH2S(ndepths) + interval * (D_dtH2S / tort2(ndepths) * 2 * ((dtH2S(ndepths-1) - dtH2S(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtH2Sw - dtH2S(ndepths)) ... %irrigation
        +TotR_dtH2S(ndepths));
    
    dFe_z = dFe(ndepths) + interval * (D_dFe / tort2(ndepths) * 2 * ((dFe(ndepths-1) - dFe(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dFew - dFe(ndepths)) ... %irrigation
        +TotR_dFe(ndepths));
    
    dMn_z = dMn(ndepths) + interval * (D_dMn / tort2(ndepths) * 2 * ((dMn(ndepths-1) - dMn(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dMnw - dMn(ndepths)) ... %irrigation
        +TotR_dMn(ndepths));

    dCa_z = dCa(ndepths) + interval * (D_dCa / tort2(ndepths) * 2 * ((dCa(ndepths-1) - dCa(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dCaw - dCa(ndepths)) ... %irrigation
        +0);

    proc_z = proc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (proc(ndepths-1) - proc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*proc(ndepths-1) + sigma(ndepths)*proc(ndepths))/z_res(ndepths)); %advection
    
    psoc_z = psoc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (psoc(ndepths-1) - psoc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*psoc(ndepths-1) + sigma(ndepths)*psoc(ndepths))/z_res(ndepths)... %advection
        +TotR_psoc(ndepths));
    
    pfoc_z = pfoc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pfoc(ndepths-1) - pfoc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pfoc(ndepths-1) + sigma(ndepths)*pfoc(ndepths))/z_res(ndepths)... %advection
        +TotR_pfoc(ndepths));
    
    pFeOH3_z = pFeOH3(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pFeOH3(ndepths-1) - pFeOH3(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pFeOH3(ndepths-1) + sigma(ndepths)*pFeOH3(ndepths))/z_res(ndepths)... %advection
        +TotR_pFeOH3(ndepths));
    
    pMnO2_z = pMnO2(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pMnO2(ndepths-1) - pMnO2(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pMnO2(ndepths-1) + sigma(ndepths)*pMnO2(ndepths))/z_res(ndepths)... %advection
        +TotR_pMnO2(ndepths));
    
    
    %% all other depths
    % ndepths=100 seems to be the sweet spot where loop and logical
    % approaches take about the same time as each other. For greater
    % ndepths, logical is faster. Indices are defined once, before the
    % loop begins. // MPH [v20]
    
    % Oxygen
    dO2_j = dO2(j);
    dO2_jp1 = dO2(jp1);
    dO2_jm1 = dO2(jm1);
    dO2(j) = dO2_j + interval*(TotR_dO2(j) - ...
        (u_j - D_dO2*DFF_j).*(dO2_jp1 - dO2_jm1)./(2*z_res_j) + ...
        (D_dO2./tort2_j).*((dO2_jp1 - 2*dO2_j + dO2_jm1)./z_res2_j) + ...
        alpha_j.*(dO2w - dO2_j));
    
    % Total alkalinity
    dtalk_j = dtalk(j);
    dtalk_jp1 = dtalk(jp1);
    dtalk_jm1 = dtalk(jm1);
    dtalk(j) = dtalk_j + interval*(TotR_dtalk(j) - ...
        (u_j - D_dtalk*DFF_j).*(dtalk_jp1 - dtalk_jm1)./(2*z_res_j) + ...
        (D_dtalk./tort2_j).*((dtalk_jp1 - 2*dtalk_j + dtalk_jm1)./z_res2_j) + ...
        alpha_j.*(dtalkw - dtalk_j));

    % Dissolved inorganic carbon
    dtCO2_j = dtCO2(j);
    dtCO2_jp1 = dtCO2(jp1);
    dtCO2_jm1 = dtCO2(jm1);
    dtCO2(j) = dtCO2_j + interval*(TotR_dtCO2(j) - ...
        (u_j - D_dtCO2*DFF_j).*(dtCO2_jp1 - dtCO2_jm1)./(2*z_res_j) + ...
        (D_dtCO2./tort2_j).*((dtCO2_jp1 - 2*dtCO2_j + dtCO2_jm1)./z_res2_j) + ...
        alpha_j.*(dtCO2w - dtCO2_j));

    % Nitrate
    dtNO3_j = dtNO3(j);
    dtNO3_jp1 = dtNO3(jp1);
    dtNO3_jm1 = dtNO3(jm1);
    dtNO3(j) = dtNO3_j + interval*(TotR_dtNO3(j) - ...
        (u_j - D_dtNO3*DFF_j).*(dtNO3_jp1 - dtNO3_jm1)./(2*z_res_j) + ...
        (D_dtNO3./tort2_j).*((dtNO3_jp1 - 2*dtNO3_j + dtNO3_jm1)./z_res2_j) + ...
        alpha_j.*(dtNO3w - dtNO3_j));
    
    % Sulfate
    dtSO4_j = dtSO4(j);
    dtSO4_jp1 = dtSO4(jp1);
    dtSO4_jm1 = dtSO4(jm1);
    dtSO4(j) = dtSO4_j + interval*(TotR_dtSO4(j) - ...
        (u_j - D_dtSO4*DFF_j).*(dtSO4_jp1 - dtSO4_jm1)./(2*z_res_j) + ...
        (D_dtSO4./tort2_j).*((dtSO4_jp1 - 2*dtSO4_j + dtSO4_jm1)./z_res2_j) + ...
        alpha_j.*(dtSO4w - dtSO4_j));
    
    % Phosphate
    dtPO4_j = dtPO4(j);
    dtPO4_jp1 = dtPO4(jp1);
    dtPO4_jm1 = dtPO4(jm1);
    dtPO4(j) = dtPO4_j + interval*(TotR_dtPO4(j) - ...
        (u_j - D_dtPO4*DFF_j).*(dtPO4_jp1 - dtPO4_jm1)./(2*z_res_j) + ...
        (D_dtPO4./tort2_j).*((dtPO4_jp1 - 2*dtPO4_j + dtPO4_jm1)./z_res2_j) + ...
        alpha_j.*(dtPO4w - dtPO4_j));
    
    % Ammonia
    dtNH4_j = dtNH4(j);
    dtNH4_jp1 = dtNH4(jp1);
    dtNH4_jm1 = dtNH4(jm1);
    dtNH4(j) = dtNH4_j + interval*(TotR_dtNH4(j) - ...
        (u_j - D_dtNH4*DFF_j).*(dtNH4_jp1 - dtNH4_jm1)./(2*z_res_j) + ...
        (D_dtNH4./tort2_j).*((dtNH4_jp1 - 2*dtNH4_j + dtNH4_jm1)./z_res2_j) + ...
        alpha_j.*(dtNH4w - dtNH4_j));
    
    % Hydrogen sulfide
    dtH2S_j = dtH2S(j);
    dtH2S_jp1 = dtH2S(jp1);
    dtH2S_jm1 = dtH2S(jm1);
    dtH2S(j) = dtH2S_j + interval*(TotR_dtH2S(j) - ...
        (u_j - D_dtH2S*DFF_j).*(dtH2S_jp1 - dtH2S_jm1)./(2*z_res_j) + ...
        (D_dtH2S./tort2_j).*((dtH2S_jp1 - 2*dtH2S_j + dtH2S_jm1)./z_res2_j) + ...
        alpha_j.*(dtH2Sw - dtH2S_j));
    
    % Dissolved iron
    dFe_j = dFe(j);
    dFe_jp1 = dFe(jp1);
    dFe_jm1 = dFe(jm1);
    dFe(j) = dFe_j + interval*(TotR_dFe(j) - ...
        (u_j - D_dFe*DFF_j).*(dFe_jp1 - dFe_jm1)./(2*z_res_j) + ...
        (D_dFe./tort2_j).*((dFe_jp1 - 2*dFe_j + dFe_jm1)./z_res2_j) + ...
        alpha_j.*(dFew - dFe_j));
    
    %Dissolved manganese
    dMn_j = dMn(j);
    dMn_jp1 = dMn(jp1);
    dMn_jm1 = dMn(jm1);
    dMn(j) = dMn_j + interval*(TotR_dMn(j) - ...
        (u_j - D_dMn*DFF_j).*(dMn_jp1 - dMn_jm1)./(2*z_res_j) + ...
        (D_dMn./tort2_j).*((dMn_jp1 - 2*dMn_j + dMn_jm1)./z_res2_j) + ...
        alpha_j.*(dMnw - dMn_j));
    
    % Dissolved calcium
    dCa_j = dCa(j);
    dCa_jp1 = dCa(jp1);
    dCa_jm1 = dCa(jm1);
    dCa(j) = dCa_j + interval*( - ...
        (u_j - D_dCa*DFF_j).*(dCa_jp1 - dCa_jm1)./(2*z_res_j) + ...
        (D_dCa./tort2_j).*((dCa_jp1 - 2*dCa_j + dCa_jm1)./z_res2_j) + ...
        alpha_j.*(dCaw - dCa_j));
    
    % Refractory organic carbon
    proc_j = proc(j);
    proc_jp1 = proc(jp1);
    proc_jm1 = proc(jm1);
    proc(j) = proc_j + interval*(... 
        - APPW_j.*(((1 - sigma_j).*proc_jp1 + ...
        2*sigma_j.*proc_j - ...
        (1 + sigma_j).*proc_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((proc_jp1 - 2*proc_j + ...
        proc_jm1)./z_res2_j));

    % Slow decay organic carbon
    psoc_j = psoc(j);
    psoc_jp1 = psoc(jp1);
    psoc_jm1 = psoc(jm1);
    psoc(j) = psoc_j + interval*(TotR_psoc(j) - ...
        APPW_j.*(((1 - sigma_j).*psoc_jp1 + ...
        2*sigma_j.*psoc_j - ...
        (1 + sigma_j).*psoc_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((psoc_jp1 - 2*psoc_j + ...
        psoc_jm1)./z_res2_j));
    
    % Fast decay organic carbon
    pfoc_j = pfoc(j);
    pfoc_jp1 = pfoc(jp1);
    pfoc_jm1 = pfoc(jm1);
    pfoc(j) = pfoc_j + interval*(TotR_pfoc(j) - ...
        APPW_j.*(((1 - sigma_j).*pfoc_jp1 + ...
        2*sigma_j.*pfoc_j - ...
        (1 + sigma_j).*pfoc_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pfoc_jp1 - 2*pfoc_j + ...
        pfoc_jm1)./z_res2_j));
    
            % Iron oxide
    pFeOH3_j = pFeOH3(j);
    pFeOH3_jp1 = pFeOH3(jp1);
    pFeOH3_jm1 = pFeOH3(jm1);
    pFeOH3(j) = pFeOH3_j + interval*(TotR_pFeOH3(j) - ...
        APPW_j.*(((1 - sigma_j).*pFeOH3_jp1 + ...
        2*sigma_j.*pFeOH3_j - ...
        (1 + sigma_j).*pFeOH3_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pFeOH3_jp1 - 2*pFeOH3_j + ...
        pFeOH3_jm1)./z_res2_j));    
    
            % Manganese oxide
    pMnO2_j = pMnO2(j);
    pMnO2_jp1 = pMnO2(jp1);
    pMnO2_jm1 = pMnO2(jm1);
    pMnO2(j) = pMnO2_j + interval*(TotR_pMnO2(j) - ...
        APPW_j.*(((1 - sigma_j).*pMnO2_jp1 + ...
        2*sigma_j.*pMnO2_j - ...
        (1 + sigma_j).*pMnO2_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pMnO2_jp1 - 2*pMnO2_j + ...
        pMnO2_jm1)./z_res2_j));    
    
    %% Set top and bottom conditions in arrays
    % Doing this here means that the correct values are used in calculating
    % diffusion and advection at all other depths // MPH [v20]
    dO2(1) = dO2_1;
    dtalk(1) = dtalk_1; 
    dtCO2(1) = dtCO2_1; 
    dtNO3(1) = dtNO3_1;
    dtSO4(1) = dtSO4_1;
    dtPO4(1) = dtPO4_1;
    dtNH4(1) = dtNH4_1;
    dtH2S(1) = dtH2S_1;
    dFe(1) = dFe_1; 
    dMn(1) = dMn_1;
    dCa(1) = dCa_1;
    proc(1) = proc_1;
    psoc(1) = psoc_1;
    pfoc(1) = pfoc_1; 
    pFeOH3(1) = pFeOH3_1;
    pMnO2(1) = pMnO2_1; 
    
    dO2(ndepths) = dO2_z;
    dtalk(ndepths) = dtalk_z; 
    dtCO2(ndepths) = dtCO2_z; 
    dtNO3(ndepths) = dtNO3_z;
    dtSO4(ndepths) = dtSO4_z;
    dtPO4(ndepths) = dtPO4_z;
    dtNH4(ndepths) = dtNH4_z;
    dtH2S(ndepths) = dtH2S_z;
    dFe(ndepths) = dFe_z; 
    dCa(ndepths) = dCa_z;
    dMn(ndepths) = dMn_z;
    proc(ndepths) = proc_z;
    psoc(ndepths) = psoc_z;
    pfoc(ndepths) = pfoc_z; 
    pFeOH3(ndepths) = pFeOH3_z;
    pMnO2(ndepths) = pMnO2_z; 
    
    %% set very small or negative concentration to zero
    dO2(dO2<0)=0;
    dtalk(dtalk<0)=0;
    dtCO2(dtCO2<0)=0;
    dtNO3(dtNO3<0)=0;
    dtSO4(dtSO4<0)=0; 
    dtPO4(dtPO4<0)=0;
    dtNH4(dtNH4<0)=0; 
    dtH2S(dtH2S<0)=0; 
    dFe(dFe<0)=0; 
    dMn(dMn<0)=0;    
    dCa(dCa<0)=0;
    proc(proc<0)=0;
    psoc(psoc<0)=0;
    pfoc(pfoc<0)=0;
    pFeOH3(pFeOH3<0)=0;
    pMnO2(pMnO2<0)=0;    
    
    %% save data every step
    %dO2f(:, i+1) = dO2;
    %dtalkf(:, i+1) = dtalk;
    %dtCO2f(:, i+1) = dtCO2; 
    %dtNO3f(:, i+1) = dtNO3;
    %dtSO4f(:, i+1) = dtSO4;
    %dtPO4f(:, i+1) = dtPO4;
    %dtNH4f(:, i+1) = dtNH4;
    %dtH2Sf(:, i+1) = dtH2S;
    %dFef(:, i+1) = dFe; 
    %dMnf(:, i+1) = dMn;
    %dCaf(:, i+1) = dCa;
    %procf(:, i+1) = proc;
    %psocf(:, i+1) = psoc;
    %pfocf(:, i+1) = pfoc; 
    %pFeOH3f(:, i+1) = pFeOH3;
    %pMnO2f(:, i+1) = pMnO2; 
    
     if i == plot_number(idx)
       disp(plot_number(idx)*interval)
       dO2f(:, idx) = dO2;
       dtalkf(:, idx) = dtalk;
       dtCO2f(:, idx) = dtCO2; 
       dtNO3f(:, idx) = dtNO3;
       dtSO4f(:, idx) = dtSO4;
       dtPO4f(:, idx) = dtPO4;
       dtNH4f(:, idx) = dtNH4;
       dtH2Sf(:, idx) = dtH2S;
       dFef(:, idx) = dFe; 
       dMnf(:, idx) = dMn;
       dCaf(:, idx) = dCa;
       procf(:, idx) = proc;
       psocf(:, idx) = psoc;
       pfocf(:, idx) = pfoc; 
       pFeOH3f(:, idx) = pFeOH3;
       pMnO2f(:, idx) = pMnO2; 
      idx=idx+1;
     end  
    end

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%RADIplot
