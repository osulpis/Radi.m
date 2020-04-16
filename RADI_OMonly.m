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
D_dO2=0.034862+0.001409*T; %[m2/a] oxygen diffusion coefficient from Li and Gregiry

%% bioturbation (for solids)
D_bio_0=1e-4*0.0232*(FoC*1e2)^0.85; %[m2/a] surf bioturb coeff, Archer et al (2002)
lambda_b = 0.08;
D_bio=D_bio_0*exp(-(depths./lambda_b).^2).*((dO2w/1e-3)/((dO2w/1e-3)+20)); %[m2/a] bioturb coeff, Archer et al (2002)

%% irrigation (for solutes)
alpha_0=11*(atan((5*FoC*1e2-400)/400)/pi+0.5)-0.9...
    +20*((dO2w/1e-3)/((dO2w/1e-3)+10))*exp(-(dO2w/1e-3)/10)*FoC*1e2/(FoC*1e2+30);    %[/a] from Archer et al (2002)
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

kslow_0=1e-4 * (FoC*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
lambda_slow=1;     %[m] tuned parameter
kslow=kslow_0*exp(-depths./lambda_slow);    %[/a] from Archer et al (2002)

kfast_0=1e-2 * (FoC*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
lambda_fast=0.03;     %[m] tuned parameter
kfast=kfast_0*exp(-depths./lambda_fast);    %[/a] from Archer et al (2002)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% R.A.D.I. main loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerun == 0 %concentrations set to zero for solids, bottom water values for not-too-sensitive solutes
    dO2(1,1:ndepths)=0;            
    proC(1,1:ndepths)=0;
    psoC(1,1:ndepths)=0;             
    pfoC(1,1:ndepths)=0;             
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime:t_length;  %we will only keep the variables every year
    plot_number(1)=1;
    
elseif rerun==1 %if it is a rerun, initial conditions are concentrations from last time step
    dO2=dO2f(idx-1,:);            %[mol/m3]
    proC=proC(idx-1,:);                %[mol/m3]
    psoC=psoC(idx-1,:);                %[mol/m3]
    pfoC=pfoC(idx-1,:);                %[mol/m3]
    plot_number=0:t_length/stoptime:t_length;  %we will only keep the variables every year
    i=plot_number(idx-1);
    
else
    % initial condition for solutes: bottom-water value
    dO2=dO2ic;                %[mol/m3]
    proC=proC_ic;
    psoC=psoC_ic;
    pfoC=pfoC_ic;
    
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
proC(:) = 3e4;
psoC(:) = 3e3;
pfoC(:) = 3e2;

% Preallocate saving arrays
dO2f = NaN(ndepths, t_length);
proCf = NaN(ndepths, t_length);
psoCf = NaN(ndepths, t_length);
pfoCf = NaN(ndepths, t_length);

for i=i:t_length-1

%     disp(i)
    
    %F_O2i=D_O2*phi(1)*(O2(:,1)-O2w)./5e-3;
    %F_DICi=D_DIC*phi(1)*(DIC(:,1)-DICw)./5e-3;
    %F_TAi=D_TA*phi(1)*(TA(:,1)-TAw)./5e-3;
    
    %% Redox reaction rates
    Rs = psoC.*kslow;
    Rf = pfoC.*kfast;
    
    %% Calculate all reactions (19 species, units: [mol/m3/a])
    % This section ~2x faster by not putting all the reactions into a
    % single matrix but keeping as separate vectors // MPH
    TotR_dO2 = -phiS./phi.*(Rs + Rf);
    TotR_psoC=-Rs;
    TotR_pfoC=-Rf;
    
    %% top boundary condition: prescribed solid fluxes and diffusive boundary layer control on solutes
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    
    %delete// OS
    %O2_1 = O2(1) + interval * ( D_O2 / tort2(1) * (2*O2(2) - 2*O2(1) + TR(1) * (O2w - O2(1))) / (z_res(1)^2) ... %diffusion
    %    - u(1) * -1 * TR(1) * ( O2w - O2(1)) / (2*z_res(1)) ... %advection
    %   + alpha(1) * ( O2w - O2(1) ) ... %irrigation
    %   + TotR_O2(1)); %reaction
    
    %implement nonzero delta_phi at the interface // OS
    dO2_1 = dO2(1) + interval * ( D_dO2 / tort2(1) * (2*dO2(2) - 2*dO2(1) + TR(1) * (dO2w - dO2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dO2.*DFF(1)) * -1 * TR(1) * ( dO2w - dO2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dO2w - dO2(1) ) ... %irrigation
        + TotR_dO2(1)); %reaction
    
    %delete//OS
    %OC_1 = OC(1) + interval * (D_bio(1) * ( 2 * OC(2) - 2 * OC(1) +... %diffusion
    %    2 * z_res(1) * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
    %   -w(1) * -1 * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1))... %advection
    %   +TotR_OC(1)); %reaction

    %implement nonzero sigma at the interface // OS
    proC_1 = proC(1) + interval * (D_bio(1) * ( 2 * proC(2) - 2 * proC(1) +... %diffusion
        2 * z_res(1) * (FroC - phiS(1) * w(1) * proC(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        + (delta_D_bio(1) + D_bio(1) / phiS(1) * - delta_phiS(1) - w(1)) * -1 * (FroC - phiS(1) * w(1) * proC(1)) / (D_bio(1) * phiS(1))); %advection

    %implement nonzero sigma at the interface // OS
    psoC_1 = psoC(1) + interval * (D_bio(1) * ( 2 * psoC(2) - 2 * psoC(1) +... %diffusion
        2 * z_res(1) * (FsoC - phiS(1) * w(1) * psoC(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        + (delta_D_bio(1) + D_bio(1) / phiS(1) * - delta_phiS(1) - w(1)) * -1 * (FsoC - phiS(1) * w(1) * psoC(1)) / (D_bio(1) * phiS(1))... %advection
        +TotR_psoC(1)); %reaction

    %implement nonzero sigma at the interface // OS
    pfoC_1 = pfoC(1) + interval * (D_bio(1) * ( 2 * pfoC(2) - 2 * pfoC(1) +... %diffusion
        2 * z_res(1) * (FfoC - phiS(1) * w(1) * pfoC(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        + (delta_D_bio(1) + D_bio(1) / phiS(1) * - delta_phiS(1) - w(1)) * -1 * (FfoC - phiS(1) * w(1) * pfoC(1)) / (D_bio(1) * phiS(1))... %advection
        +TotR_pfoC(1)); %reaction    
    
%     if i == 1
%         disp(' ')
%         disp('original irrigative O2 term at top:')
%         disp(interval*alpha(1) * ( O2w - O2(1) ))
%         disp(' ')
%     end % if
%     if i == 1
%         disp(' ')
%         disp('original advective OC term at top:')
%         disp(interval*(-w(1) * -1 * (2*Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1))))
%         disp('corrected advective OC term at top:')
%         disp(interval*(-w(1) * -1 * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1))))
%         disp(' ')
%     end % if
%     if i == 1
%         disp(' ')
%         disp('original diffusive OC term at top:')
%         disp(interval*(D_bio(1) * ( 2 * OC(2) - 2 * OC(1) +... %diffusion
%             2 * z_res(1) * (2*Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1)) ...
%             ) / (z_res(1).^2)))
%         disp('corrected diffusive OC term at top:')
%         disp(interval*(D_bio(1) * ( 2 * OC(2) - 2 * OC(1) +... %diffusion
%             2 * z_res(1) * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1)) ...
%             ) / (z_res(1).^2)))
%         disp(' ')
%     end % if
      
    %% bottom boundary condition: gradients disappear
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    dO2_z = dO2(ndepths) + interval * (D_dO2 / tort2(ndepths) * 2 * ((dO2(ndepths-1) - dO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dO2w - dO2(ndepths)) ... %irrigation
        +TotR_dO2(ndepths));

    proC_z = proC(ndepths) + interval * (D_bio(ndepths) * 2 * ( (proC(ndepths-1) - proC(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*proC(ndepths-1) + sigma(ndepths)*proC(ndepths))/z_res(ndepths)); %advection
    
    psoC_z = psoC(ndepths) + interval * (D_bio(ndepths) * 2 * ( (psoC(ndepths-1) - psoC(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*psoC(ndepths-1) + sigma(ndepths)*psoC(ndepths))/z_res(ndepths)... %advection
        +TotR_psoC(ndepths));
    
    pfoC_z = pfoC(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pfoC(ndepths-1) - pfoC(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pfoC(ndepths-1) + sigma(ndepths)*pfoC(ndepths))/z_res(ndepths)... %advection
        +TotR_pfoC(ndepths));
    
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

    % Refractory organic carbon
    proC_j = proC(j);
    proC_jp1 = proC(jp1);
    proC_jm1 = proC(jm1);
    proC(j) = proC_j + interval*(... 
        - APPW_j.*(((1 - sigma_j).*proC_jp1 + ...
        2*sigma_j.*proC_j - ...
        (1 + sigma_j).*proC_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((proC_jp1 - 2*proC_j + ...
        proC_jm1)./z_res2_j));

        % Slow decay organic carbon
    psoC_j = psoC(j);
    psoC_jp1 = psoC(jp1);
    psoC_jm1 = psoC(jm1);
    psoC(j) = psoC_j + interval*(TotR_psoC(j) - ...
        APPW_j.*(((1 - sigma_j).*psoC_jp1 + ...
        2*sigma_j.*psoC_j - ...
        (1 + sigma_j).*psoC_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((psoC_jp1 - 2*psoC_j + ...
        psoC_jm1)./z_res2_j));
    
        % Fast decay organic carbon
    pfoC_j = pfoC(j);
    pfoC_jp1 = pfoC(jp1);
    pfoC_jm1 = pfoC(jm1);
    pfoC(j) = pfoC_j + interval*(TotR_pfoC(j) - ...
        APPW_j.*(((1 - sigma_j).*pfoC_jp1 + ...
        2*sigma_j.*pfoC_j - ...
        (1 + sigma_j).*pfoC_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pfoC_jp1 - 2*pfoC_j + ...
        pfoC_jm1)./z_res2_j));    
    
    %% Set top and bottom conditions in arrays
    % Doing this here means that the correct values are used in calculating
    % diffusion and advection at all other depths // MPH [v20]
    dO2(1) = dO2_1;
    proC(1) = proC_1;
    psoC(1) = psoC_1;
    pfoC(1) = pfoC_1;
    dO2(ndepths) = dO2_z;
    proC(ndepths) = proC_z;
    psoC(ndepths) = psoC_z;
    pfoC(ndepths) = pfoC_z;
    
    %% set very small or negative concentration to zero
    dO2(dO2<0)=0;
    proC(proC<0)=0;
    psoC(psoC<0)=0;
    pfoC(pfoC<0)=0;
    
    %% save data every step
    dO2f(:, i+1) = dO2;
    proCf(:, i+1) = proC;
    psoCf(:, i+1) = psoC;
    pfoCf(:, i+1) = pfoC;
    
end

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

RADIplot
