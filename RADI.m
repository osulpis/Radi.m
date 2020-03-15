%% 1-D Reaction-Advection-Diffusion-Irrigation (RADI) Diagenetic Sediment Module
%% Source code by O. Sulpis and M. Wilhelmus
%% Optimised for MATLAB by M.P. Humphreys [March 2020]
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
CO2SYSr = CO2SYS(TAw*1e6/rho_sw,DICw*1e6/rho_sw,1,2,S,T,T,P*10,P*10,Siw*1e6/rho_sw,PO4w*1e6/rho_sw,1,10,1);
k1(1,1:z_length) = CO2SYSr(1,67);           %carbonic acid first dissociation constant
k2(1,1:z_length) = CO2SYSr(1,68);           %carbonic acid second dissociation constant
k1p(1,1:z_length) = CO2SYSr(1,75);         %phosphate constant 1
k2p(1,1:z_length) = CO2SYSr(1,76);         %phosphate constant 2
k3p(1,1:z_length) = CO2SYSr(1,77);         %phosphate constant 3
kb(1,1:z_length) = CO2SYSr(1,72);            %boron constant 
kw(1,1:z_length) = CO2SYSr(1,71);           %water dissociation constants
ksi(1,1:z_length) = CO2SYSr(1,78);           %silica constants
bt(1,1:z_length) = CO2SYSr(1,79);             %[umol/kg] total boron 
omegaC = CO2SYSr(1,30);                         %calcite saturation state
omegaA = CO2SYSr(1,31);                          %aragonite saturation state
co3 = CO2SYSr(1,22) .* 10^-6;                   %[mol/kg] CO3 
hco3 = CO2SYSr(1,21) .* 10^-6;                 %[mol/kg] HCO3
ph = CO2SYSr(1,37) ;                                   %pH on the total scale
CALK = Caw ./ rho_sw;                                 %[mol/kg] Ca concentration 
fg(1,1:z_length)=TAw./rho_sw-hco3-2*co3; %sum of all alkalinity species that are not carbon
kspc = (co3 .* CALK) ./ omegaC;                %[mol2/kg2] calcite in situ solubility
kspa = (co3 .* CALK) ./ omegaA;                %[mol2/kg2] aragonite in situ solubility
ff(1,1:z_length) = 1;                                        %random parameter needed for calc_pco2
H(1,1:z_length) = 10^-ph;                             %[mol/kg] H concentration first guess
clear co3 hco3 ph omegaC omegaA CALK CO2SYSr   
sit = 120 * 10^-6;                                              %[mol/kg] convert silica concentration
bt = bt .* 10^-6;                                                %[mol/kg] convert boron concentration

%% temperature dependent "free solution" diffusion coefficients
D_TA=0.015169+0.000793*T;           %[m2/a] approximated to bicarbonate diffusion coefficient from Hulse et al (2018)
D_DIC=0.015169+0.000793*T;         %[m2/a] approximted to bicarbonate diffusion coefficient from Hulse et al (2018)
D_Ca=0.0107+0.001677*T;               %[m2/a] calcium diffusion coefficient from Li and Gregory
D_O2=0.034862+0.001409*T;           %[m2/a] oxygen diffusion coefficient from Li and Gregiry
D_NO3=0.030842+0.001226*T;        %[m2/a] nitrate diffusion coefficient from Li and Gregory (1974)
D_SO4=0.015768+0.000788*T;        %[m2/a] sulfate diffusion coefficient from Li and Gregory (1974)
D_PO4=0.011291+0.000559*T;        %[m2/a] phosphate diffusion coefficient from Li and Gregory (1974)
D_NH4=0.030905+0.001226*T;        %[m2/a] ammonium diffusion coefficient from Li and Gregory (1974)
D_H2S=0.030748+0.000964*T;        %[m2/a] hydrogen sulfide diffusion coefficient from the UNISENSE table by Ramsing and Gundersen
D_Mn2p=0.0086+0.001525*T;           %[m2/a] manganese diffusion coefficient from Li and Gregory (1974)
D_Fe2p=0.0108+0.001478*T;           %[m2/a] iron diffusion coefficient from Li and Gregory (1974)

%% bioturbation (for solids)
D_bio_0=1e-4*0.0232*(Foc*1e2)^0.85;                                              %[m2/a] surf bioturb coeff, Archer et al (2002)
D_bio=D_bio_0*exp(-(z./0.08).^2).*((O2w/1e-3)/((O2w/1e-3)+20)); %[m2/a] bioturb coeff, Archer et al (2002)

%% irrigation (for solutes)
alpha_0=11*(atan((5*Foc*1e2-400)/400)/3.1416+0.5)-0.9...
    +20*((O2w/1e-3)/((O2w/1e-3)+10))*exp(-(O2w/1e-3)/10)*Foc*1e2/(Foc*1e2+30);    %[/a] from Archer et al (2002)
alpha=alpha_0.*exp(-(z/0.05).^2);                                                                                   %[/a] Archer et al (2002) the depth of 5 cm was changed

%% depth-dependent porosity and diffusion coefficient loss
delta_phi = [0 diff(phi)]; % depth-dependent porosity loss
delta_phiS = [0 diff(phiS)]; % depth-dependent solid fraction gain
delta_tort2 = [0 diff(tort.^2)]; % depth-dependent tortuosity gain
delta_D_bio = [0 diff(D_bio)]; % [m/a]

% biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis (1977)
Peh=w.*z_res./(2*D_bio);      %one half the cell Peclet number (Eq. 97 in Boudreau 1996)
% when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
sigma=1./tanh(Peh)-1./(Peh);  %Eq. 96 in Boudreau 1996

%% organic matter degradation parameters
 KO2=0.003;            %[mol/m3] Monod constant from Soetaert et al. 1996 (GCA)
 KinO2=0.01;           %[mol/m3] Monod inhibition constant from Soetaert et al. 1996 (GCA)
 KNO3=0.03;           %[mol/m3] Monod constant from Soetaert et al. 1996 (GCA)
 KinNO3=0.005;       %[mol/m3] Monod inhibition constant from Soetaert et al. 1996 (GCA)
 KMnO2=42.4;            %[mol/m3] Monod constant from Van Cappellen and Wang 1996
 KFeOH3=265;         %[mol/m3] Monod constant from Van Cappellen and Wang 1996 
 KSO4=1.6;                 %[mol/m3] Monod constant from Van Cappellen and Wang 1996
 KinMnO2=KMnO2;
 KinFeOH3=KFeOH3;
 KinSO4=KSO4;
 
klabile=93750*D_bio_0*exp(-z./0.03);                 %[/a] from Archer et al (2002), consumption scale of 3mm
krefractory=80.25*D_bio_0*exp(-z./1);         %[/a] from Archer et al (2002)

%% redox reaction first order constants
slowfac=1;                           %slowing factor for redox reactions to avoid numerical issues
kmnox=1e6./slowfac;                %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
kfeox=1e6./slowfac;                   %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
knhox=1e4./slowfac;                  %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
ksox=3e5./slowfac;                    %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)

%% adsorption parameters
KNA=1.78;                    %[no unit] adsorption coefficient for ammonia in W Atlantic sediments from Mackin and Aller (LO, 1984)
KPA=500;                     %[no unit] adsorption coefficient for phosphate in oxic marine sediments from Krom and Berner (LO, 1980)
phiN=phi+phiS.*KNA;  %ammonia adsorption correction
phiP=phi+phiS.*KPA;  %phosphate adsorption correction

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% R.A.D.I. main loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerun == 0 %concentrations set to zero for solids, bottom water values for not-too-sensitive solutes
    TA(1,1:z_length)=TAw;
    DIC(1,1:z_length)=DICw;
    OmegaC=NaN(1,z_length);
    O2(1,1:z_length)=0;            
    Ca(1,1:z_length)=Caw;
    NO3(1,1:z_length)=0;
    SO4(1,1:z_length)=0;
    PO4(1,1:z_length)=PO4w;
    NH4(1,1:z_length)=0;
    H2S(1,1:z_length)=0;
    Mn2p(1,1:z_length)=0;
    Fe2p(1,1:z_length)=0;
    % initial condition for solids
    Calcite(1,1:z_length)=0;                    
    Aragonite(1,1:z_length)=0;          
    OC_labile(1,1:z_length)=0;             
    OC_refractory(1,1:z_length)=0;   
    MnO2(1,1:z_length)=0;                       
    FeOH3(1,1:z_length)=0;                     
    Clay(1,1:z_length)=0;   
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/t_end:t_length;  %we will only keep the variables every year
    plot_number(1)=1;
    
elseif rerun==1 %if it is a rerun, initial conditions are concentrations from last time step
    TA=TAf(idx-1,:);            %[mol/m3]
    DIC=DICf(idx-1,:);         %[mol/m3]
    OmegaC=OmegaCf(idx-1,:);              %[mol/m3]
    O2=O2f(idx-1,:);            %[mol/m3]
    Ca=Caf(idx-1,:);            %[mol/m3]
    NO3=NO3f(idx-1,:);      %[mol/m3]
    SO4=SO4f(idx-1,:);      %[mol/m3]
    PO4=PO4f(idx-1,:);      %[mol/m3]
    NH4=NH4f(idx-1,:);       %[mol/m3]
    H2S=H2Sf(idx-1,:);       %[mol/m3]
    Mn2p=Mn2pf(idx-1,:);   %[mol/m3]
    Fe2p=Fe2pf(idx-1,:);    %[mol/m3]
    % initial condition for solids
    Calcite=Calcitef(idx-1,:);                           %[mol/m3]
    Aragonite=Aragonitef(idx-1,:);          
    OC_labile=OC_labilef(idx-1,:);                %[mol/m3]
    OC_refractory=OC_refractoryf(idx-1,:);   %[mol/m3]
    MnO2=MnO2f(idx-1,:);                               %[mol/m3]
    FeOH3=FeOH3f(idx-1,:);                           %[mol/m3]
    Clay=Clayf(idx-1,:);                                      %[mol/m3]
    plot_number=0:t_length/t_end:t_length;  %we will only keep the variables every year
    i=plot_number(idx-1);
    
else
    % initial condition for solutes: bottom-water value
    TA=TAic;          %[mol/m3]
    DIC=DICic;       %[mol/m3]
    OmegaC=OmegaCic;   %[mol/m3]
    O2=O2ic;                %[mol/m3]
    Ca=Caic;          %[mol/m3]
    NO3=NO3ic;             %[mol/m3]
    SO4=SO4ic;             %[mol/m3]
    PO4=PO4ic;     %[mol/m3]
    NH4=NH4ic;              %[mol/m3]
    H2S=H2Sic;              %[mol/m3]
    Mn2p=Mn2pic;            %[mol/m3]
    Fe2p=Fe2pic;            %[mol/m3]
    % initial condition for solids
    Calcite=Calciteic;                    
    Aragonite=Aragoniteic;          
    OC_labile=OC_labileic;             
    OC_refractory=OC_refractoryic;   
    MnO2=MnO2ic;                       
    FeOH3=FeOH3ic;                     
    Clay=Clayic;   
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/t_end:t_length;  %we will only keep the variables every year
    plot_number(1)=1;
    
end

%% short-cut transport variables
APPW=w - delta_D_bio - delta_phiS.* D_bio./ phiS;
DFF=(tort2.* delta_phi./ phi - delta_tort2)./ (tort2.^2);
DBF=phiS.*D_bio+D_bio.*(-delta_phi);
TR=(2*z_res.* (tort.^2)./ dbl);

%% Prepare for timestep calculations // MPH [v20]

% Set indices for depth-varying reactions
j = 2:z_length-1;
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

% Preallocate dissolution rate vectors
Rd_aragonite = zeros(1, z_length);
Rd_calcite = zeros(1, z_length);
Rp_calcite = zeros(1, z_length);

% % Precalculations - these don't actually speed things up!
% uDFF_TA = u_j - D_TA*DFF_j;
% Dtort2_TA = D_TA./tort2_j;

%% Begin timestep loop
for i=i:t_length-1

    %F_O2i=D_O2*phi(1)*(O2(:,1)-O2w)./5e-3;
    %F_DICi=D_DIC*phi(1)*(DIC(:,1)-DICw)./5e-3;
    %F_TAi=D_TA*phi(1)*(TA(:,1)-TAw)./5e-3;
    
    %% calc_pco2: time efficient carbonate system solver    
    DIC_molPerKg = DIC ./ rho_sw;     %convert DIC to [mol/kg]
    TA_molPerKg = TA ./ rho_sw;         %convert TA to [mol/kg]
    PO4_molPerKg = PO4 ./ rho_sw;   %convert PO4 to [mol/kg]
    
    [H] = calc_pCO2(DIC_molPerKg,PO4_molPerKg,sit,bt,TA_molPerKg,ff,k1,k2,k1p,k2p,k3p,kb,kw,ksi,fg,H);    %[mol/kg] H concentration
    
    co3_molPerKg = (DIC_molPerKg .* k1 .* k2)./ (H.^2 + (k1 .* H) + (k1 .* k2)); %[mol/kg] CO3 concentration
    co3 = co3_molPerKg .* rho_sw; %[mol m^-3] CO3 concentraiton
    
    %% CaCO3 reactions
    OmegaC = Ca.*co3./ (kspc.* rho_sw.^2); %[no unit] calcite saturation state
    OmegaA = Ca.*co3./ (kspa.* rho_sw.^2); %[no unit] aragonite saturation state

%     %Aragonite dissolution from Dong et al. (2018) GCA
%     Rd_aragonite=Aragonite.*2.*(1-OmegaA).^1.4;
    % This Rd_aragonite not used so now commented out // MPH [v20]
    
    % Aragonite dissolution rate from Dong et al. (2019) EPSL
    % Converted from loops to logicals // MPH [v20]
    Rd_aragonite(:) = 0;
    % 0.835 is the OmegaA value at which both laws are equal
    J = OmegaA > 0.835 & OmegaA <= 1;
    if any(J)
        Rd_aragonite(J) = (Aragonite(J)*100.09*4.* ...
            (0.0019*(1 - OmegaA(J)).^0.13))/25;
    end % if
    J = OmegaA <= 0.835;
    if any(J)
        Rd_aragonite(J) = (Aragonite(J)*100.09*4.* ...
            (0.021*(1 - OmegaA(J)).^1.46))/25;
    end % if
    
    %Aragonite precipitation
        %To be implemented
        %OmegaA(OmegaA>1)=1;
    
    %Calcite dissolution from Naviaux et al. (2019) Marine Chemistry
    %data for forams with SSA of 4 m2/g
    % Converted from loops to logicals // MPH [v20]
    Rd_calcite(:) = 0;
     % 0.8275 is the OmegaC value at which both laws are equal
    J = OmegaC > 0.8275 & OmegaC <= 1;
    if any(J)
        Rd_calcite(J) = (Calcite(J)*100.09*4.* ...
            (0.00158*(1 - OmegaC(J)).^0.11))/25;
    end % if
    J = OmegaC <= 0.8275;
    if any(J)
        Rd_calcite(J) = (Calcite(J)*100.09*4.* ...
            (5*(1 - OmegaC(J)).^4.7))/25;
    end % if

    %Calcite precipitation rate from Zuddas and Mucci, GCA (1998)
    %normalized to the same surface area than for dissolution (4m2/g)
    % Converted from loops to logicals // MPH [v20]
    Rp_calcite(:) = 0;
    J = Calcite < 23500 & OmegaC > 1;
    if any(J)
%         Rp_calcite(J) = (Calcite(J)*100.09*4*1.63.* ...
%             (OmegaC(J) - 1).^1.76)/100;
        Rp_calcite(J) = 1.63*(OmegaC(J) - 1).^1.76;
    end
   
    %% Organic matter respiration pathways
    fO2=O2./(KO2+O2);                   %from the code of Couture et al. (EST 2010), following Boudreau (1996)
    fNO3=NO3./(KNO3+NO3).*(KinO2./(KinO2+O2));
    fMnO2=MnO2./(MnO2+KMnO2).*(KinNO3./(KinNO3+NO3)).*(KinO2./(KinO2+O2));
    fFe3=FeOH3./(FeOH3+KFeOH3).*(KinMnO2./(MnO2+KinMnO2)).*(KinNO3./(KinNO3+NO3)).*(KinO2./(KinO2+O2));
    fSO4=SO4./(SO4+KSO4).*(KinFeOH3./(FeOH3+KinFeOH3)).*(KinMnO2./(MnO2+KinMnO2)).*(KinNO3./(KinNO3+NO3)).*(KinO2./(KinO2+O2));
    fCH4=(KinSO4./(SO4+KinSO4)).*(KinFeOH3./(FeOH3+KinFeOH3)).*(KinMnO2./(MnO2+KinMnO2)).*(KinNO3./(KinNO3+NO3)).*(KinO2./(KinO2+O2));
    fOx=fO2+fNO3+fMnO2+fFe3+fSO4+fCH4;
    fNH=fOx-fO2;
    
    %% Redox reaction rates
    Rg = OC_labile.*klabile + OC_refractory.*krefractory;
    Rfeox = kfeox.* Fe2p.* O2;
    Rmnox = kmnox.* Mn2p.* O2;
    Rsox = ksox.* H2S.* O2;
    Rnhox = knhox.* NH4.* O2;  
    
    %% Calculate all reactions (19 species, units: [mol/m3/a])
    % This section ~2x faster by not putting all the reactions into a
    % single matrix but keeping as separate vectors // MPH
    TotR_TA = phiS./phi.*Rg.*((RN/RC - RP/RC).*fO2 + ...
        (0.8 + RN/RC - RP/RC).*fNO3 + (4 + RN/RC - RP/RC).*fMnO2 + ...
        (8 + RN/RC - RP/RC).*fFe3 + (1 + RN/RC - RP/RC).*fSO4 + ...
        (RN/RC - RP/RC).*fCH4) + phiS./phi.*2.* ...
        (Rd_calcite + Rd_aragonite - Rp_calcite) - ...
        2*Rfeox - 2*Rmnox - 2*Rsox - 2*Rnhox;
    TotR_DIC = phiS./phi.*Rg.* ...
        (fO2 + fNO3 + fMnO2 + fFe3 + fSO4 + fCH4/2) + ...
        phiS./phi.*(Rd_calcite + Rd_aragonite - Rp_calcite);
    TotR_calcite = -Rd_calcite + Rp_calcite;
    TotR_aragonite = -Rd_aragonite;
    TotR_O2 = -phiS./phi.*Rg.*fO2 - 2*Rsox - 2*Rnhox - 0.5*Rmnox - ...
        0.25*Rfeox;
    TotR_OC_labile = -OC_labile.*klabile.*fOx;
    TotR_OC_refractory = -OC_refractory.*krefractory.*fOx;
    TotR_Ca = phiS./phi.*(Rd_calcite + Rd_aragonite - Rp_calcite);
    TotR_MnO2 = -2*Rg.*fMnO2 + phi./phiS.*Rmnox;
    TotR_FeOH3 = -4*Rg.*fFe3 + phi./phiS.*Rfeox;
	TotR_NO3 = -phiS./phi.*Rg.*0.8.*fNO3 + Rnhox;
	TotR_SO4 = -phiS./phi.*Rg.*fSO4/2 + Rsox;
	TotR_PO4 = phiS.*(RP/RC).*Rg.*fOx;
	TotR_NH4 = phiS.*(RN/RC).*Rg.*fOx - Rnhox;
	TotR_H2S = phiS./phi.*Rg.*fSO4/2 - Rsox;
    TotR_Mn2p = phiS./phi.*2.*Rg.*fMnO2 - Rmnox;
    TotR_Fe2p = phiS./phi.*4.*Rg.*fFe3 - Rfeox;
    
    %% top boundary condition: prescribed solid fluxes and diffusive boundary layer control on solutes
    TA(1) = TA(1) + t_res * ( D_TA(1) / tort2(1) * (2*TA(2) - 2*TA(1) + TR(1) * (TAw - TA(1))) / (z_res(1).^2) ... %diffusion
        - u(1) *  -1 * TR(1) * (TAw - TA(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * (TAw - TA(1)) ... %irrigation
        + TotR_TA(1)); %reaction

    DIC(1) = DIC(1) + t_res * ( D_DIC(1) / tort2(1) * (2*DIC(2) - 2*DIC(1) + TR(1) * (DICw - DIC(1))) / (z_res(1).^2) ... %diffusion
        - u(1) *  -1 * TR(1) * (DICw - DIC(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * (DICw - DIC(1)) ... %irrigation
        + TotR_DIC(1)); %reaction

    Calcite(1) = Calcite(1) + t_res * (D_bio(1) * ( 2 * Calcite(2) - 2 * Calcite(1) +... %diffusion
        2 * z_res(1) * (Fcalcite - phiS(1) * w(1) * Calcite(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
        - w(1) * -1 * (Fcalcite - phiS(1) * w(1) * Calcite(1)) / (D_bio(1) * phiS(1))... %advection
        + TotR_calcite(1)); %reaction

    Aragonite(1) = Aragonite(1) + t_res * (D_bio(1) * ( 2 * Aragonite(2) - 2 * Aragonite(1) +... %diffusion
        2 * z_res(1) * (Faragonite - phiS(1) * w(1) * Aragonite(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
        - w(1) * -1 * (Faragonite - phiS(1) * w(1) * Aragonite(1)) / (D_bio(1) * phiS(1))... %advection
        + TotR_aragonite(1)); %reaction

    O2(1) = O2(1) + t_res * ( D_O2 / tort2(1) * (2*O2(2) - 2*O2(1) + TR(1) * (O2w - O2(1))) / (z_res(1)^2) ... %diffusion
        - u(1) * -1 * TR(1) * ( O2w - O2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( O2w - O2(1) ) ... %irrigation
        + TotR_O2(1)); %reaction

    OC_labile(1) = OC_labile(1) + t_res * (D_bio(1) * ( 2 * OC_labile(2) - 2 * OC_labile(1) +... %diffusion
        2 * z_res(1) * (Flabile - phiS(1) * w(1) * OC_labile(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
        -w(1) * -1 * (Flabile - phiS(1) * w(1) * OC_labile(1)) / (D_bio(1) * phiS(1))... %advection
        + TotR_OC_labile(1)); %reaction

    OC_refractory(1) = OC_refractory(1) + t_res * (D_bio(1) * ( 2 * OC_refractory(2) - 2 * OC_refractory(1) +... %diffusion
        2 * z_res(1) * (2*Frefractory - phiS(1) * w(1) * OC_refractory(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        -w(1) * -1 * (2*Frefractory - phiS(1) * w(1) * OC_refractory(1)) / (D_bio(1) * phiS(1))... %advection
        +TotR_OC_refractory(1)); %reaction

    Ca(1) = Ca(1) + t_res * ( D_Ca(1) / tort2(1) * (2*Ca(2) - 2*Ca(1) + TR(1) * (Caw - Ca(1))) / (z_res(1).^2) ... %diffusion
        - u(1) *  -1 * TR(1) * (Caw - Ca(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * (Caw - Ca(1)) ... %irrigation
        + TotR_Ca(1)); %reaction

    MnO2(1) = MnO2(1) + t_res * (D_bio(1) * ( 2 * MnO2(2) - 2 * MnO2(1) +... %diffusion
        2 * z_res(1) * (Fmno2 - phiS(1) * w(1) * MnO2(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
        - w(1) * -1 * (Fmno2 - phiS(1) * w(1) * MnO2(1)) / (D_bio(1) * phiS(1))... %advection
        + TotR_MnO2(1)); %reaction                                                                                  

    FeOH3(1) = FeOH3(1) + t_res * (D_bio(1) * ( 2 * FeOH3(2) - 2 * FeOH3(1) +... %diffusion
        2 * z_res(1) * (Ffeoh3 - phiS(1) * w(1) * FeOH3(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
        - w(1) * -1 * (Ffeoh3 - phiS(1) * w(1) * FeOH3(1)) / (D_bio(1) * phiS(1))... %advection
        + TotR_FeOH3(1)); %reaction

    NO3(1) = NO3(1) + t_res * ( D_NO3(1) / tort2(1) * (2*NO3(2) - 2*NO3(1) + TR(1) * (NO3w - NO3(1))) / (z_res(1).^2) ... %diffusion
        - u(1) *  -1 * TR(1) * (NO3w - NO3(1)) / (2*z_res(1)) ...  %advection
        + alpha(1) * (NO3w - NO3(1))... %irrigation
        + TotR_NO3(1));  %reaction

    SO4(1) = SO4(1) + t_res * ( D_SO4(1) / tort2(1) * (2*SO4(2) - 2*SO4(1) + TR(1) * (SO4w - SO4(1))) / (z_res(1).^2) ...  %diffusion
        - u(1) *  -1 * TR(1) * (SO4w - SO4(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * (SO4w - SO4(1))... %irrigation
        + TotR_SO4(1)); %reaction

    PO4(1) = PO4(1) + t_res * ( (phi(1) * D_PO4 / tort2(1) + phiS(1) * D_bio(1) * KPA) * ( 2 * PO4(2) - 2 * PO4(1)... %diffusion
        + TR(1) * (PO4w - PO4(1) )) / (z_res(1).^2) ... %diffusion
        - (phi(1) * u(1) + phiS(1) * w(1) * KPA - phi(1) * D_PO4 * DFF(1) - KPA * DBF(1)) * -1 * TR(1)... %advection
        * (PO4w - PO4(1)) / (2*z_res(1)) ... %advection
        + phi(1) * alpha(1) * (PO4w - PO4(1)) ... %irrigation
        +TotR_PO4(1))./ (phi(z_length)+phiS(z_length)*KPA); %reaction

    NH4(1) = NH4(1) + t_res * ( (phi(1) * D_NH4 / tort2(1) + phiS(1) * D_bio(1) * KNA) * ( 2 * NH4(2) - 2 * NH4(1)... %diffusion
        + TR(1) * (NH4w - NH4(1) )) / (z_res(1).^2) ...  %diffusion
        - (phi(1) * u(1) + phiS(1) * w(1) * KNA - phi(1) * D_NH4 * DFF(1) - KNA * DBF(1)) * -1 * TR(1)... %advection
        * (NH4w - NH4(1)) / (2*z_res(1)) ... %advection
        + phi(1) * alpha(1) * (NH4w - NH4(1)) ...  %irrigation
        + TotR_NH4(1)) / (phi(z_length)+phiS(z_length)*KNA); %reaction

    H2S(1) = H2S(1) + t_res * ( D_H2S(1) / tort2(1) * (2*H2S(2) - 2*H2S(1) + TR(1) * (H2Sw - H2S(1))) / (z_res(1).^2) ...  %diffusion
        - u(1) *  -1 * TR(1) * (H2Sw - H2S(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * (H2Sw - H2S(1)) ... %irrigation
        + ( phiS(1) / phi(1) ) * (OC_labile(1) * klabile(1) + OC_refractory(1) * krefractory(1) ) * fSO4(1) / 2 ...  %reaction
        +TotR_H2S(1)); %reaction

    Mn2p(1) = Mn2p(1) + t_res * ( D_Mn2p(1) / tort2(1) * (2*Mn2p(2) - 2*Mn2p(1) + TR(1) * (Mn2pw - Mn2p(1))) / (z_res(1).^2) ...  %diffusion
        - u(1) *  -1 * TR(1) * (Mn2pw - Mn2p(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * (Mn2pw - Mn2p(1)) ...   %irrigation
        + ( phiS(1) / phi(1) ) * 2 * (OC_labile(1) * klabile(1) + OC_refractory(1) * krefractory(1) ) * fMnO2(1)... %reaction
        +TotR_Mn2p(1)); %reaction

    Fe2p(1) = Fe2p(1) + t_res * ( D_Fe2p(1) / tort2(1) * (2*Fe2p(2) - 2*Fe2p(1) + TR(1) * (Fe2pw - Fe2p(1))) / (z_res(1).^2) ... %diffusion
        - u(1) *  -1 * TR(1) * (Fe2pw - Fe2p(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * (Fe2pw - Fe2p(1)) ... %irrigation
        + TotR_Fe2p(1)); %reaction

    Clay(1) = Clay(1) + t_res * (D_bio(1) * ( 2 * Clay(2) - 2 * Clay(1) +... %diffusion
        2 * z_res(1) * (Fclay - phiS(1) * w(1) * Clay(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
        - w(1) * -1 * (Fclay - phiS(1) * w(1) * Clay(1)) / (D_bio(1) * phiS(1))); %advection
      
    %% bottom boundary condition: gradients disappear
    TA(z_length) = TA(z_length) + t_res * (D_TA / tort2(z_length) * 2 * ((TA(z_length-1) - TA(z_length)) / z_res(z_length).^2) ...   %diffusion
        +alpha(z_length) * (TAw - TA(z_length)) ... %irrigation
        +TotR_TA(z_length));

    DIC(z_length) = DIC(z_length) + t_res * (D_DIC / tort2(z_length) * 2 * ((DIC(z_length-1) - DIC(z_length)) / z_res(z_length).^2) ...   %diffusion
        +alpha(z_length) * (DICw - DIC(z_length)) ...  %irrigation
        +TotR_DIC(z_length));

    Calcite(z_length) = Calcite(z_length) + t_res * (D_bio(z_length) * 2 * ( (Calcite(z_length-1) - Calcite(z_length))  / z_res(z_length).^2)...%diffusion
        - APPW(z_length) * (-sigma(z_length)*Calcite(z_length-1) + sigma(z_length)*Calcite(z_length))/z_res(z_length)...  %advection
        +TotR_calcite(z_length));

    Aragonite(z_length) = Aragonite(z_length) + t_res * (D_bio(z_length) * 2 * ( (Aragonite(z_length-1) - Aragonite(z_length))  / z_res(z_length).^2)...%diffusion
        - APPW(z_length) * (-sigma(z_length)*Aragonite(z_length-1) + sigma(z_length)*Aragonite(z_length))/z_res(z_length)... %advection
        +TotR_aragonite(z_length));

    O2(z_length) = O2(z_length) + t_res * (D_O2 / tort2(z_length) * 2 * ((O2(z_length-1) - O2(z_length)) / z_res(z_length).^2) ...  %diffusion
        + alpha(z_length) * (O2w - O2(z_length)) ... %irrigation
        +TotR_O2(z_length));

    OC_labile(z_length) = OC_labile(z_length) + t_res * (D_bio(z_length) * 2 * ( (OC_labile(z_length-1) - OC_labile(z_length)) / z_res(z_length).^2)... %diffusion
        - APPW(z_length) * (-sigma(z_length)*OC_labile(z_length-1) + sigma(z_length)*OC_labile(z_length))/z_res(z_length)... %advection
        +TotR_OC_labile(z_length));

    OC_refractory(z_length) = OC_refractory(z_length) + t_res * (D_bio(z_length) * 2 * ( (OC_refractory(z_length-1) - OC_refractory(z_length)) / z_res(z_length).^2)... %diffusion
        - APPW(z_length) * (-sigma(z_length)*OC_refractory(z_length-1) + sigma(z_length)*OC_refractory(z_length))/z_res(z_length)... %advection
        +TotR_OC_refractory(z_length));

    Ca(z_length) = Ca(z_length) + t_res * (D_Ca / tort2(z_length) * 2 * ((Ca(z_length-1) - Ca(z_length)) / z_res(z_length).^2) ... %diffusion
        +alpha(z_length) * (Caw - Ca(z_length)) ... %irrigation
        +TotR_Ca(z_length));

    MnO2(z_length) = MnO2(z_length) + t_res * (D_bio(z_length) * 2 * ( (MnO2(z_length-1) - MnO2(z_length))  / z_res(z_length).^2)... %diffusion
        - APPW(z_length) * (-sigma(z_length) * MnO2(z_length-1) + sigma(z_length)*MnO2(z_length))/z_res(z_length)...  %advection
        +TotR_MnO2(z_length));

    FeOH3(z_length) = FeOH3(z_length) + t_res * (D_bio(z_length) * 2 * ( (FeOH3(z_length-1) - FeOH3(z_length))  / z_res(z_length).^2)...%diffusion
        - APPW(z_length) * (-sigma(z_length)*FeOH3(z_length-1) + sigma(z_length)*FeOH3(z_length))/z_res(z_length)...%advection
        +TotR_FeOH3(z_length));

    NO3(z_length) = NO3(z_length) + t_res * (D_NO3 / tort2(z_length) * 2 * ((NO3(z_length-1) - NO3(z_length)) / z_res(z_length).^2) ...  %diffusion
        + alpha(z_length) * (NO3w - NO3(z_length)) ...  %irrigation
        +TotR_NO3(z_length));

    SO4(z_length) = SO4(z_length) + t_res * (D_SO4 / tort2(z_length) * 2 * ((SO4(z_length-1) - SO4(z_length)) / z_res(z_length).^2) ...  %diffusion
        + alpha(z_length) * (SO4w - SO4(z_length)) ... %irrigation
        +TotR_SO4(z_length));

    PO4(z_length) = PO4(z_length) + t_res * ( (phi(z_length) * D_PO4 / tort2(z_length) + phiS(z_length) * D_bio(z_length) * KPA) ... %diffusion
        * 2 * ((PO4(z_length-1) - PO4(z_length)) / z_res(z_length).^2) ...  %diffusion
        + phi(z_length) * alpha(z_length) * (PO4w - PO4(z_length)) ... %irrigation
        +TotR_PO4(z_length)) / (phi(z_length)+phiS(z_length)*KPA); %reaction

    NH4(z_length) = NH4(z_length) + t_res * ( (phi(z_length) * D_NH4 / tort2(z_length) + phiS(z_length) * D_bio(z_length) * KNA) ... %diffusion
        * 2 * ((NH4(z_length-1) - NH4(z_length)) / z_res(z_length).^2) ...     %diffusion
        + phi(z_length) * alpha(z_length) * (NH4w - NH4(z_length)) ...  %irrigation
        +TotR_NH4(z_length)) / (phi(z_length)+phiS(z_length)*KNA);  %reaction

    H2S(z_length) = H2S(z_length) + t_res * (D_H2S / tort2(z_length) * 2 * ((H2S(z_length-1) - H2S(z_length)) / z_res(z_length).^2) ...  %diffusion
        + alpha(z_length) * (H2Sw - H2S(z_length)) ... %irrigation
        +TotR_H2S(z_length));

    Mn2p(z_length) = Mn2p(z_length) + t_res * (D_Mn2p / tort2(z_length) * 2 * ((Mn2p(z_length-1) - Mn2p(z_length)) / z_res(z_length).^2) ...  %diffusion
        + alpha(z_length) * (Mn2pw - Mn2p(z_length)) ... %irrigation
        +TotR_Mn2p(z_length));

    Fe2p(z_length) = Fe2p(z_length) + t_res * (D_Fe2p / tort2(z_length) * 2 * ((Fe2p(z_length-1) - Fe2p(z_length)) / z_res(z_length).^2) ...  %diffusion
        + alpha(z_length) * (Fe2pw - Fe2p(z_length)) ...   %irrigation
        +TotR_Fe2p(z_length));

    Clay(z_length) = Clay(z_length) + t_res * (D_bio(z_length) * 2 * ((Clay(z_length-1) - Clay(z_length))  / z_res(z_length).^2) ...%diffusion
    - APPW(z_length) * (-sigma(z_length)*Clay(z_length-1) + sigma(z_length)*Clay(z_length))/z_res(z_length));... %advection    
    
    %% all other depths
    % z_length=100 seems to be the sweet spot where loop and logical
    % approaches take about the same time as each other. For greater
    % z_length, logical is faster. Indices are defined once, before the
    % loop begins. // MPH [v20]
    
    % Total alkalinity
    TA_j = TA(j);
    TA_jp1 = TA(jp1);
    TA_jm1 = TA(jm1);
    TA(j) = TA_j + t_res*(TotR_TA(j) - ...
        (u_j - D_TA*DFF_j).*(TA_jp1 - TA_jm1)./(2*z_res_j) + ...
        (D_TA./tort2_j).*((TA_jp1 - 2*TA_j + TA_jm1)./z_res2_j) + ...
        alpha_j.*(TAw - TA_j));
    % ^For some reason, it's actually (very slightly) slower to put
    % (u_j - D_TA*DFF_j) and (D_TA./tort2_j) into variables that are
    % evaluated before the loop. Also, putting the TA calculation into a
    % function is ~10x slower. // MPH
   
    % Dissolved inorganic carbon
    DIC_j = DIC(j);
    DIC_jp1 = DIC(jp1);
    DIC_jm1 = DIC(jm1);
    DIC(j) = DIC_j + t_res*(TotR_DIC(j) - ...
        (u_j - D_DIC*DFF_j).*(DIC_jp1 - DIC_jm1)./(2*z_res_j) + ...
        (D_DIC./tort2_j).*((DIC_jp1 - 2*DIC_j + DIC_jm1)./z_res2_j) + ...
        alpha_j.*(DICw - DIC_j));
   
    % Calcite
    Calcite_j = Calcite(j);
    Calcite_jp1 = Calcite(jp1);
    Calcite_jm1 = Calcite(jm1);
    Calcite(j) = Calcite_j + t_res*(TotR_calcite(j) - ...
        APPW_j.*(((1 - sigma_j).*Calcite_jp1 + ...
        2*sigma_j.*Calcite_j - (1 + sigma_j).*Calcite_jm1)./ ...
        (2*z_res_j)) + D_bio_j.* ...
        ((Calcite_jp1 - 2*Calcite_j + Calcite_jm1)./z_res2_j));

    % Aragonite
    Aragonite_j = Aragonite(j);
    Aragonite_jp1 = Aragonite(jp1);
    Aragonite_jm1 = Aragonite(jm1);
    Aragonite(j) = Aragonite_j + t_res*(TotR_aragonite(j) - ...
        APPW_j.*(((1 - sigma_j).*Aragonite_jp1 + ...
        2*sigma_j.*Aragonite_j - (1 + sigma_j).*Aragonite_jm1)./ ...
        (2*z_res_j)) + D_bio_j.* ...
        ((Aragonite_jp1 - 2*Aragonite_j + Aragonite_jm1)./z_res2_j));

    % Oxygen
    O2_j = O2(j);
    O2_jp1 = O2(jp1);
    O2_jm1 = O2(jm1);
    O2(j) = O2_j + t_res*(TotR_O2(j) - ...
        (u_j - D_O2*DFF_j).*(O2_jp1 - O2_jm1)./(2*z_res_j) + ...
        (D_O2./tort2_j).*((O2_jp1 - 2*O2_j + O2_jm1)./z_res2_j) + ...
        alpha_j.*(O2w - O2_j));
    
    % Labile organic carbon
    OC_labile_j = OC_labile(j);
    OC_labile_jp1 = OC_labile(jp1);
    OC_labile_jm1 = OC_labile(jm1);
    OC_labile(j) = OC_labile_j + t_res*(TotR_OC_labile(j) - ...
        APPW_j.*(((1 - sigma_j).*OC_labile_jp1 + ...
        2*sigma_j.*OC_labile_j - (1 + sigma_j).*OC_labile_jm1)./ ...
        (2*z_res_j)) + D_bio_j.* ...
        ((OC_labile_jp1 - 2*OC_labile_j + OC_labile_jm1)./z_res2_j));

    % Refractory organic carbon
    OC_refractory_j = OC_refractory(j);
    OC_refractory_jp1 = OC_refractory(jp1);
    OC_refractory_jm1 = OC_refractory(jm1);
    OC_refractory(j) = OC_refractory_j + t_res*(TotR_OC_refractory(j) - ...
        APPW_j.*(((1 - sigma_j).*OC_refractory_jp1 + ...
        2*sigma_j.*OC_refractory_j - ...
        (1 + sigma_j).*OC_refractory_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((OC_refractory_jp1 - 2*OC_refractory_j + ...
        OC_refractory_jm1)./z_res2_j));

    % Calcium
    Ca_j = Ca(j);
    Ca_jp1 = Ca(jp1);
    Ca_jm1 = Ca(jm1);
    Ca(j) = Ca_j + t_res*(TotR_Ca(j) - ...
        (u_j - D_Ca*DFF_j).*(Ca_jp1 - Ca_jm1)./(2*z_res_j) + ...
        (D_Ca./tort2_j).*((Ca_jp1 - 2*Ca_j + Ca_jm1)./(z_res2_j)) + ...
        alpha_j.*(Caw - Ca_j));

    % Manganese dioxide
    MnO2_j = MnO2(j);
    MnO2_jp1 = MnO2(jp1);
    MnO2_jm1 = MnO2(jm1);
    MnO2(j) = MnO2_j + t_res*(TotR_MnO2(j) - ...
        APPW_j.*(((1 - sigma_j).*MnO2_jp1 + 2*sigma_j.*MnO2_j - ...
        (1 + sigma_j).*MnO2_jm1)./(2*z_res_j)) + D_bio_j.* ...
        ((MnO2_jp1 - 2*MnO2_j + MnO2_jm1)./z_res2_j));

    % Iron(III) hydroxide
    FeOH3_j = FeOH3(j);
    FeOH3_jp1 = FeOH3(jp1);
    FeOH3_jm1 = FeOH3(jm1);
    FeOH3(j) = FeOH3_j + t_res*(TotR_FeOH3(j) - ...
        APPW_j.*(((1 - sigma_j).*FeOH3_jp1 + 2*sigma_j.*FeOH3_j - ...
        (1 + sigma_j).*FeOH3_jm1)./(2*z_res_j)) + D_bio_j.* ...
        ((FeOH3_jp1 - 2*FeOH3_j + FeOH3_jm1)./z_res2_j));

    % Nitrate
    NO3_j = NO3(j);
    NO3_jp1 = NO3(jp1);
    NO3_jm1 = NO3(jm1);
    NO3(j) = NO3_j + t_res*(TotR_NO3(j) - ...
       (u_j - D_NO3*DFF_j).*(NO3_jp1 - NO3_jm1)./(2*z_res_j) + ...
       (D_NO3./tort2_j).*((NO3_jp1 - 2*NO3_j + NO3_jm1)./z_res2_j) + ...
       alpha_j.*(NO3w - NO3_j));
    
    % Sulfate
    SO4_j = SO4(j);
    SO4_jp1 = SO4(jp1);
    SO4_jm1 = SO4(jm1);
    SO4(j) = SO4_j + t_res*(TotR_SO4(j) - ...
        (u_j - D_SO4*DFF_j).*(SO4_jp1 - SO4_jm1)./(2*z_res_j) + ...
        (D_SO4./tort2_j).*((SO4_jp1 - 2*SO4_j + SO4_jm1)./z_res2_j) + ...
        alpha_j.*(SO4w - SO4_j));

    % Phosphate
    PO4_j = PO4(j);
    PO4_jp1 = PO4(jp1);
    PO4_jm1 = PO4(jm1);
    PO4(j) = PO4_j + t_res*(TotR_PO4(j) - ...
        (u_j - D_PO4*DFF_j).*(PO4_jp1 - PO4_jm1)./(2*z_res_j) + ...
        (D_PO4./tort2_j).*((PO4_jp1 - 2*PO4_j + PO4_jm1)./z_res2_j) + ...
        alpha_j.*(PO4w - PO4_j))./(phi_j + phiS_j*KPA);
    % ^ Switched final division to elementwise not matrix // MPH [v20]

    % Ammonium
    NH4_j = NH4(j);
    NH4_jp1 = NH4(jp1);
    NH4_jm1 = NH4(jm1);
    NH4(j) = NH4_j + t_res*(TotR_NH4(j) - ...
        (u_j - D_NH4*DFF_j).*(NH4_jp1 - NH4_jm1)./(2*z_res_j) + ...
        (D_NH4./tort2_j).*((NH4_jp1 - 2*NH4_j + NH4_jm1)./z_res2_j) + ...
        alpha_j.*(NH4w - NH4_j))./(phi_j + phiS_j*KNA);
    % ^ Switched final division to elementwise not matrix // MPH [v20]

    % Hydrogen sulfide
    H2S_j = H2S(j);
    H2S_jp1 = H2S(jp1);
    H2S_jm1 = H2S(jm1);
    H2S(j) = H2S_j + t_res*(TotR_H2S(j) - ...
       (u_j - D_H2S*DFF_j).*(H2S_jp1 - H2S_jm1)./(2*z_res_j) + ...
       (D_H2S./tort2_j).*((H2S_jp1 - 2*H2S_j + H2S_jm1)./z_res2_j) + ...
       alpha_j.*(H2Sw - H2S_j));

    % Mn2p
    Mn2p_j = Mn2p(j);
    Mn2p_jp1 = Mn2p(jp1);
    Mn2p_jm1 = Mn2p(jm1);
    Mn2p(j) = Mn2p_j + t_res*(TotR_Mn2p(j) - ...
        (u_j - D_Mn2p*DFF_j).*(Mn2p_jp1 - Mn2p_jm1)./(2*z_res_j) + ...
        (D_Mn2p./tort2_j).*((Mn2p_jp1 - 2*Mn2p_j + Mn2p_jm1)./z_res2_j) ...
        + alpha_j.*(Mn2pw - Mn2p_j));
    
    % Fe2p
    Fe2p_j = Fe2p(j);
    Fe2p_jp1 = Fe2p(jp1);
    Fe2p_jm1 = Fe2p(jm1);
    Fe2p(j) = Fe2p_j + t_res*(TotR_Fe2p(j) - ...
        (u_j - D_Fe2p*DFF_j).*(Fe2p_jp1 - Fe2p_jm1)./(2*z_res_j) + ...
        (D_Fe2p./tort2_j).*((Fe2p_jp1 - 2*Fe2p_j + Fe2p_jm1)./z_res2_j) ...
        + alpha_j.*(Fe2pw - Fe2p_j));
    
    % Clay
    Clay_j = Clay(j);
    Clay_jp1 = Clay(jp1);
    Clay_jm1 = Clay(jm1);
    Clay(j) = Clay_j + t_res*(-APPW_j.*(((1 - sigma_j).*Clay_jp1 + ...
        2*sigma_j.*Clay_j - (1 + sigma_j).*Clay_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((Clay_jp1 - 2*Clay_j + Clay_jm1)./z_res2_j));

    %% set very small or negative concentration to zero
    TA(TA<0)=0;
    DIC(DIC<0)=0;
    O2(O2<0)=0;
    Ca(Ca<0)=0; % missing =0 added // MPH [v20]
    NO3(NO3<0)=0;
    SO4(SO4<0)=0;
    PO4(PO4<0)=0;
    NH4(NH4<0)=0;
    H2S(H2S<0)=0;
    Mn2p(Mn2p<0)=0;
    Fe2p(Fe2p<0)=0;
    Calcite(Calcite<0)=0;
    Aragonite(Aragonite<0)=0;
    OC_labile(OC_labile<0)=0;
    OC_refractory(OC_refractory<0)=0;
    MnO2(MnO2<0)=0;
    FeOH3(FeOH3<0)=0;
    
    %% save data every year
    if i == plot_number(idx)
        disp(plot_number(idx)*t_res)
        TAf(idx,:)=TA;
        DICf(idx,:)=DIC;
        OmegaCf(idx,:)=OmegaC;
        CO3f(idx,:)=co3;
        pHf(idx,:)=-log10(H);
        OC_labilef(idx,:)=OC_labile;
        OC_refractoryf(idx,:)=OC_refractory;
        O2f(idx,:)=O2;
        Calcitef(idx,:)=Calcite;
        Aragonitef(idx,:)=Aragonite;
        Caf(idx,:)=Ca;
        MnO2f(idx,:)=MnO2;
        FeOH3f(idx,:)=FeOH3;
        NO3f(idx,:)=NO3;
        SO4f(idx,:)=SO4;
        PO4f(idx,:)=PO4;
        NH4f(idx,:)=NH4;
        H2Sf(idx,:)=H2S;
        Mn2pf(idx,:)=Mn2p;
        Fe2pf(idx,:)=Fe2p;
        Clayf(idx,:)=Clay;
        idx=idx+1;
    end

    if i==217%any(isnan(OC_labile))
        disp(i)
        break
    end % if
end

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
