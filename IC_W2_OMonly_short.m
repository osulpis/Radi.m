%Station W2
%Hammon et al 1996 Deep-Sea Research
Station= "Hammond1996 - W2";

%% definition of the spatial domain with two different resolutions
z_end=20e-2;                               %[m] bottom sediment depth, should be a multiple of z_res
z_res=2e-2;                               %[m] depth resolution
z_length = 1 + z_end/z_res;                %[no unit] number of depth layers
z = linspace(0, z_end, z_length); %[m] depth
z_res = ones(size(z))*z_res; %[m] depth resolution

%% definition of the temporal domain
% t_end=20000;                             %[a] total timespan of the problem
t_end = 5/8760;
t_res=1/8760;                          %[a] time resolution (1/60000 is nine minutes, 1/8760 is one hour; 1/365.2 is a day)
% t_res = 1/8760;
t_length=t_end/t_res;                 %[no unit] number of time layers

%% bottom-water environmental conditions
T=1.4;                   %[C] temperature
SF_depth=4310; %[m] seafloor depth
S=34.69;              %[psu] salinity
rho_sw = gsw_rho(S,T,1); %[kg/m^3] in situ seawater density computed from GSW toolbox
P=rho_sw.*9.81.*SF_depth/1e5; %[bar] in situ pressure computed from GSW toolbox

%% bottom-water values of dissolved species
O2w=(159.7)*1e-6*rho_sw;                                   %[mol/m3] dissolved oxygen
PO4w=(2.39)*1e-6*rho_sw;                                   %[mol/m3] nitrate from GLODAP at sation location, bottom waters

%% depth-dependent porosity
phi = (0.85-0.74)*exp(-33.*z)+0.74;         %porosity profile (porewater bulk fraction) fitted from station7 mooring3 of cruise NBP98-2 by Sayles et al DSR 2001
phiS=1-phi;                                          %solid volume fraction
tort=(1-2*log(phi)).^0.5;                      %tortuosity from Boudreau (1996, GCA)
tort2=tort.^2;                                        %tortuosity squared

%% Redfield ratios
RC=1/(6.9*1e-3*PO4w./(1e-6*rho_sw)+6*1e-3);      %P:C computed as a function of SRP from Galbraith and Martiny PNAS 2015
RN=11;                    % value at 60 degS from Martiny et al. Nat G 2013 
RP=1;                          % Redfield ratio for P in the deep sea
M_OM=30.03+(RN/RC)*17.03+(RP/RC)*97.994; %[g of OM per mol of OC] Organic Matter molar mass

%% solid fluxes and solid initial conditions
Foc=1.000041991200773; %[mol/m2/a] flux of total organic carbon to the bottom 

Ftot=Foc.*M_OM;      %[g/m2/a] total sediment flux 
v0=(Ftot)/(2.65e6*phiS(1));                                             %[m/a] bulk burial velocity at sediment-water interface
vinf=v0*phiS(1)/phiS(1,z_length);                                    %[m/a]bulk burial velocity at the infinite depth
for j=1:z_length
    u(1,j)=vinf*phi(1,z_length)/phi(1,j);                               %[m/a] porewater burial velocity
    w(1,j)=vinf*phiS(1,z_length)/phiS(1,j);                         %[m/a] solid burial velocity
end

%% diffusive boundary layer
dbl=1e-3;            %[m] thickness at location taken from Sulpis et al 2018 PNAS
save('data/IC_W2_OMonly_short.mat')

rerun = 0;
