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

%% temperature dependent "free solution" diffusion coefficients
D_dO2=0.034862+0.001409*T; %[m2/a] oxygen diffusion coefficient from Li and Gregory
D_dtCO2=0.015169+0.000793*T;         %[m2/a] approximted to bicarbonate diffusion coefficient from Hulse et al (2018)

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

kslow_0=1e-4 * (Foc*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
lambda_slow=1;     %[m] tuned parameter
kslow=kslow_0*exp(-depths./lambda_slow);    %[/a] from Archer et al (2002)

kfast_0=1e-2 * (Foc*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
lambda_fast=0.03;     %[m] tuned parameter
kfast=kfast_0*exp(-depths./lambda_fast);    %[/a] from Archer et al (2002)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% R.A.D.I. main loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerun == 0 %concentrations set to zero for solids, bottom water values for not-too-sensitive solutes
    dO2(1,1:ndepths)=0;            
    dtCO2(1,1:ndepths)=0;            
    proc(1,1:ndepths)=0;
    psoc(1,1:ndepths)=0;             
    pfoc(1,1:ndepths)=0;             
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime:t_length;  %we will only keep the variables every year
    plot_number(1)=1;
    
elseif rerun==1 %if it is a rerun, initial conditions are concentrations from last time step
    dO2=dO2f(idx-1,:);            %[mol/m3]
    dtCO2(1,1:ndepths)=0;            
    proc=proc(idx-1,:);                %[mol/m3]
    psoc=psoc(idx-1,:);                %[mol/m3]
    pfoc=pfoc(idx-1,:);                %[mol/m3]
    plot_number=0:t_length/stoptime:t_length;  %we will only keep the variables every year
    i=plot_number(idx-1);
    
else
    % initial condition for solutes: bottom-water value
    dO2=dO2ic;                %[mol/m3]
    dtCO2=dtCO2ic;                %[mol/m3]
    proc=proc_ic;
    psoc=psoc_ic;
    pfoc=pfoc_ic;
    
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
dO2(:) = dO2w*2/3;
dtCO2(:) = dtCO2w;
proc(:) = 3e4;
psoc(:) = 3e3;
pfoc(:) = 3e2;

% Preallocate saving arrays
dO2f = NaN(ndepths, t_length);
dtCO2f = NaN(ndepths, t_length);
procf = NaN(ndepths, t_length);
psocf = NaN(ndepths, t_length);
pfocf = NaN(ndepths, t_length);

for i=i:t_length-1

%     disp(i)
    
    %F_O2i=D_O2*phi(1)*(O2(:,1)-O2w)./5e-3;
    %F_DICi=D_DIC*phi(1)*(DIC(:,1)-DICw)./5e-3;
    %F_TAi=D_TA*phi(1)*(TA(:,1)-TAw)./5e-3;
    
    %% Redox reaction rates
    Rs = psoc.*kslow;
    Rf = pfoc.*kfast;
    
    %% Calculate all reactions (19 species, units: [mol/m3/a])
    % This section ~2x faster by not putting all the reactions into a
    % single matrix but keeping as separate vectors // MPH
    TotR_dO2 = -phiS./phi.*(Rs + Rf);
    TotR_dtCO2 = phiS./phi.*(Rs + Rf);
    TotR_psoc=-Rs;
    TotR_pfoc=-Rf;
    
    %% top boundary condition: prescribed solid fluxes and diffusive boundary layer control on solutes
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    
    dO2_1 = dO2(1) + interval * ( D_dO2 / tort2(1) * (2*dO2(2) - 2*dO2(1) + TR(1) * (dO2w - dO2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dO2.*DFF(1)) * -1 * TR(1) * ( dO2w - dO2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dO2w - dO2(1) ) ... %irrigation
        + TotR_dO2(1)); %reaction
    
    dtCO2_1 = dtCO2(1) + interval * ( D_dtCO2 / tort2(1) * (2*dtCO2(2) - 2*dtCO2(1) + TR(1) * (dtCO2w - dtCO2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtCO2.*DFF(1)) * -1 * TR(1) * ( dtCO2w - dtCO2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtCO2w - dtCO2(1) ) ... %irrigation
        + TotR_dtCO2(1)); %reaction
    
    proc_1 = proc(1) + interval * (D_bio(1) * ( 2 * proc(2) - 2 * proc(1) +... %diffusion
        2 * z_res(1) * (Froc - phiS(1) * w(1) * proc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        + (delta_D_bio(1) + D_bio(1) / phiS(1) * - delta_phiS(1) - w(1)) * -1 * (Froc - phiS(1) * w(1) * proc(1)) / (D_bio(1) * phiS(1))); %advection

    psoc_1 = psoc(1) + interval * (D_bio(1) * ( 2 * psoc(2) - 2 * psoc(1) +... %diffusion
        2 * z_res(1) * (Fsoc - phiS(1) * w(1) * psoc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        + (delta_D_bio(1) + D_bio(1) / phiS(1) * - delta_phiS(1) - w(1)) * -1 * (Fsoc - phiS(1) * w(1) * psoc(1)) / (D_bio(1) * phiS(1))... %advection
        +TotR_psoc(1)); %reaction

    pfoc_1 = pfoc(1) + interval * (D_bio(1) * ( 2 * pfoc(2) - 2 * pfoc(1) +... %diffusion
        2 * z_res(1) * (Ffoc - phiS(1) * w(1) * pfoc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        + (delta_D_bio(1) + D_bio(1) / phiS(1) * - delta_phiS(1) - w(1)) * -1 * (Ffoc - phiS(1) * w(1) * pfoc(1)) / (D_bio(1) * phiS(1))... %advection
        +TotR_pfoc(1)); %reaction        
      
    %% bottom boundary condition: gradients disappear
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    dO2_z = dO2(ndepths) + interval * (D_dO2 / tort2(ndepths) * 2 * ((dO2(ndepths-1) - dO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dO2w - dO2(ndepths)) ... %irrigation
        +TotR_dO2(ndepths));

    dtCO2_z = dtCO2(ndepths) + interval * (D_dtCO2 / tort2(ndepths) * 2 * ((dtCO2(ndepths-1) - dtCO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtCO2w - dtCO2(ndepths)) ... %irrigation
        +TotR_dtCO2(ndepths));
    
    proc_z = proc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (proc(ndepths-1) - proc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*proc(ndepths-1) + sigma(ndepths)*proc(ndepths))/z_res(ndepths)); %advection
    
    psoc_z = psoc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (psoc(ndepths-1) - psoc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*psoc(ndepths-1) + sigma(ndepths)*psoc(ndepths))/z_res(ndepths)... %advection
        +TotR_psoc(ndepths));
    
    pfoc_z = pfoc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pfoc(ndepths-1) - pfoc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pfoc(ndepths-1) + sigma(ndepths)*pfoc(ndepths))/z_res(ndepths)... %advection
        +TotR_pfoc(ndepths));
    
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

    % Dissolved inorganic carbon
    dtCO2_j = dtCO2(j);
    dtCO2_jp1 = dtCO2(jp1);
    dtCO2_jm1 = dtCO2(jm1);
    dtCO2(j) = dtCO2_j + interval*(TotR_dtCO2(j) - ...
        (u_j - D_dtCO2*DFF_j).*(dtCO2_jp1 - dtCO2_jm1)./(2*z_res_j) + ...
        (D_dtCO2./tort2_j).*((dtCO2_jp1 - 2*dtCO2_j + dtCO2_jm1)./z_res2_j) + ...
        alpha_j.*(dtCO2w - dtCO2_j));
    
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
    
    %% Set top and bottom conditions in arrays
    % Doing this here means that the correct values are used in calculating
    % diffusion and advection at all other depths // MPH [v20]
    dO2(1) = dO2_1;
    dtCO2(1) = dtCO2_1;
    proc(1) = proc_1;
    psoc(1) = psoc_1;
    pfoc(1) = pfoc_1;
    dO2(ndepths) = dO2_z;
    dtCO2(ndepths) = dtCO2_z;
    proc(ndepths) = proc_z;
    psoc(ndepths) = psoc_z;
    pfoc(ndepths) = pfoc_z;
    
    %% set very small or negative concentration to zero
    dO2(dO2<0)=0;
    dtCO2(dtCO2<0)=0;
    proc(proc<0)=0;
    psoc(psoc<0)=0;
    pfoc(pfoc<0)=0;
    
    %% save data every step
    dO2f(:, i+1) = dO2;
    dtCO2f(:,i+1) = dtCO2;
    procf(:, i+1) = proc;
    psocf(:, i+1) = psoc;
    pfocf(:, i+1) = pfoc;
    
end

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

RADIplot
