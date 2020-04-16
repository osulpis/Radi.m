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
D_O2=0.034862+0.001409*T; %[m2/a] oxygen diffusion coefficient from Li and Gregiry

%% bioturbation (for solids)
D_bio_0=1e-4*0.0232*(Foc*1e2)^0.85; %[m2/a] surf bioturb coeff, Archer et al (2002)
lambda_b = 0.08;
D_bio=D_bio_0*exp(-(z./lambda_b).^2).*((O2w/1e-3)/((O2w/1e-3)+20)); %[m2/a] bioturb coeff, Archer et al (2002)

%% irrigation (for solutes)
alpha_0=11*(atan((5*Foc*1e2-400)/400)/pi+0.5)-0.9...
    +20*((O2w/1e-3)/((O2w/1e-3)+10))*exp(-(O2w/1e-3)/10)*Foc*1e2/(Foc*1e2+30);    %[/a] from Archer et al (2002)
alpha=alpha_0.*exp(-(z/0.05).^2);                                                                                   %[/a] Archer et al (2002) the depth of 5 cm was changed

%% depth-dependent porosity and diffusion coefficient loss
% % Use differences - values should then be divided by z_res?
% delta_phi = [0 diff(phi)]; % depth-dependent porosity loss
% delta_phiS = [0 diff(phiS)]; % depth-dependent solid fraction gain
% delta_tort2 = [0 diff(tort.^2)]; % depth-dependent tortuosity gain
% delta_D_bio_i = [0 diff(D_bio)]; % [m/a]
% Use derivative equations instead! all checked vs finite differences
delta_phi = -phiBeta.*(phi0 - phiInf).*exp(-phiBeta*z);
% delta_phi(1) = 0; % don't do this
delta_phiS = -delta_phi;
delta_tort2 = -2*delta_phi./phi; % not used in Julia
delta_D_bio = -2*z.*D_bio/lambda_b^2; % not used in Julia

% biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis (1977)
Peh=w.*z_res./(2*D_bio);      %one half the cell Peclet number (Eq. 97 in Boudreau 1996)
% when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
sigma=1./tanh(Peh)-1./(Peh);  %Eq. 96 in Boudreau 1996

%% organic matter degradation parameters
krefractory=80.25*D_bio_0*exp(-z./1);         %[/a] from Archer et al (2002)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% R.A.D.I. main loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerun == 0 %concentrations set to zero for solids, bottom water values for not-too-sensitive solutes
    O2(1,1:z_length)=0;            
    OC(1,1:z_length)=0;             
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/t_end:t_length;  %we will only keep the variables every year
    plot_number(1)=1;
    
elseif rerun==1 %if it is a rerun, initial conditions are concentrations from last time step
    O2=O2f(idx-1,:);            %[mol/m3]
    OC=OC(idx-1,:);                %[mol/m3]
    plot_number=0:t_length/t_end:t_length;  %we will only keep the variables every year
    i=plot_number(idx-1);
    
else
    % initial condition for solutes: bottom-water value
    O2=O2ic;                %[mol/m3]
    OC=OC_ic;             
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

% % Precalculations - these don't actually speed things up!
% uDFF_TA = u_j - D_TA*DFF_j;
% Dtort2_TA = D_TA./tort2_j;

%% Begin timestep loop
% Start with some O2 and OC
O2(:) = O2w*2/3;
OC(:) = 1e4;

% Preallocate saving arrays
O2f = NaN(z_length, t_length);
OCf = NaN(z_length, t_length);
O2f(:, 1) = O2;
OCf(:, 1) = OC;

for i=i:t_length-1

%     disp(i)
    
    %F_O2i=D_O2*phi(1)*(O2(:,1)-O2w)./5e-3;
    %F_DICi=D_DIC*phi(1)*(DIC(:,1)-DICw)./5e-3;
    %F_TAi=D_TA*phi(1)*(TA(:,1)-TAw)./5e-3;
    
    %% Redox reaction rates
    Rg = OC.*krefractory;
    
    %% Calculate all reactions (19 species, units: [mol/m3/a])
    % This section ~2x faster by not putting all the reactions into a
    % single matrix but keeping as separate vectors // MPH
    TotR_O2 = -phiS./phi.*Rg;
    TotR_OC=-OC.*krefractory;
    
    %% top boundary condition: prescribed solid fluxes and diffusive boundary layer control on solutes
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    
    %delete// OS
    %O2_1 = O2(1) + t_res * ( D_O2 / tort2(1) * (2*O2(2) - 2*O2(1) + TR(1) * (O2w - O2(1))) / (z_res(1)^2) ... %diffusion
    %    - u(1) * -1 * TR(1) * ( O2w - O2(1)) / (2*z_res(1)) ... %advection
    %   + alpha(1) * ( O2w - O2(1) ) ... %irrigation
    %   + TotR_O2(1)); %reaction
    
    %implement nonzero delta_phi at the interface // OS
    O2_1 = O2(1) + t_res * ( D_O2 / tort2(1) * (2*O2(2) - 2*O2(1) + TR(1) * (O2w - O2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_O2.*DFF(1)) * -1 * TR(1) * ( O2w - O2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( O2w - O2(1) ) ... %irrigation
        + TotR_O2(1)); %reaction
    
    %delete//OS
    %OC_1 = OC(1) + t_res * (D_bio(1) * ( 2 * OC(2) - 2 * OC(1) +... %diffusion
    %    2 * z_res(1) * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
    %   -w(1) * -1 * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1))... %advection
    %   +TotR_OC(1)); %reaction

    %implement nonzero sigma at the interface // OS
    OC_1 = OC(1) + t_res * (D_bio(1) * ( 2 * OC(2) - 2 * OC(1) +... %diffusion
        2 * z_res(1) * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        + (delta_D_bio(1) + D_bio(1) / phiS(1) * - delta_phiS(1) - w(1)) * -1 * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1))... %advection
        +TotR_OC(1)); %reaction
    
%     if i == 1
%         disp(' ')
%         disp('original irrigative O2 term at top:')
%         disp(t_res*alpha(1) * ( O2w - O2(1) ))
%         disp(' ')
%     end % if
%     if i == 1
%         disp(' ')
%         disp('original advective OC term at top:')
%         disp(t_res*(-w(1) * -1 * (2*Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1))))
%         disp('corrected advective OC term at top:')
%         disp(t_res*(-w(1) * -1 * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1))))
%         disp(' ')
%     end % if
%     if i == 1
%         disp(' ')
%         disp('original diffusive OC term at top:')
%         disp(t_res*(D_bio(1) * ( 2 * OC(2) - 2 * OC(1) +... %diffusion
%             2 * z_res(1) * (2*Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1)) ...
%             ) / (z_res(1).^2)))
%         disp('corrected diffusive OC term at top:')
%         disp(t_res*(D_bio(1) * ( 2 * OC(2) - 2 * OC(1) +... %diffusion
%             2 * z_res(1) * (Foc - phiS(1) * w(1) * OC(1)) / (D_bio(1) * phiS(1)) ...
%             ) / (z_res(1).^2)))
%         disp(' ')
%     end % if
      
    %% bottom boundary condition: gradients disappear
    % Calculate here, but don't set in arrays yet, otherwise calculations
    % at other depths use values from the wrong timestep // MPH [v20]
    O2_z = O2(z_length) + t_res * (D_O2 / tort2(z_length) * 2 * ((O2(z_length-1) - O2(z_length)) / z_res(z_length).^2) ...  %diffusion
        + alpha(z_length) * (O2w - O2(z_length)) ... %irrigation
        +TotR_O2(z_length));

    OC_z = OC(z_length) + t_res * (D_bio(z_length) * 2 * ( (OC(z_length-1) - OC(z_length)) / z_res(z_length).^2)... %diffusion
        - APPW(z_length) * (-sigma(z_length)*OC(z_length-1) + sigma(z_length)*OC(z_length))/z_res(z_length)... %advection
        +TotR_OC(z_length));
    
    %% all other depths
    % z_length=100 seems to be the sweet spot where loop and logical
    % approaches take about the same time as each other. For greater
    % z_length, logical is faster. Indices are defined once, before the
    % loop begins. // MPH [v20]
    
    % Oxygen
    O2_j = O2(j);
    O2_jp1 = O2(jp1);
    O2_jm1 = O2(jm1);
    O2(j) = O2_j + t_res*(TotR_O2(j) - ...
        (u_j - D_O2*DFF_j).*(O2_jp1 - O2_jm1)./(2*z_res_j) + ...
        (D_O2./tort2_j).*((O2_jp1 - 2*O2_j + O2_jm1)./z_res2_j) + ...
        alpha_j.*(O2w - O2_j));

    % Refractory organic carbon
    OC_j = OC(j);
    OC_jp1 = OC(jp1);
    OC_jm1 = OC(jm1);
    OC(j) = OC_j + t_res*(TotR_OC(j) - ...
        APPW_j.*(((1 - sigma_j).*OC_jp1 + ...
        2*sigma_j.*OC_j - ...
        (1 + sigma_j).*OC_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((OC_jp1 - 2*OC_j + ...
        OC_jm1)./z_res2_j));

    %% Set top and bottom conditions in arrays
    % Doing this here means that the correct values are used in calculating
    % diffusion and advection at all other depths // MPH [v20]
    O2(1) = O2_1;
    OC(1) = OC_1;
    O2(z_length) = O2_z;
    OC(z_length) = OC_z;
    
    %% set very small or negative concentration to zero
    O2(O2<0)=0;
    OC(OC<0)=0;
    
    %% save data every step
    O2f(:, i+1) = O2;
    OCf(:, i+1) = OC;

end

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

RADIplot
