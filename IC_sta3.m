%Station W2
%Hammon et al 1996 Deep-Sea Research
clear all

Station= "Hammond1996 - W2";

%% Definition of the spatial domain
z_max=40e-2;     %[m] bottom sediment depth, should be a multiple of z_res
z_res=2e-2;     %[m] depth resolution
ndepths = 1 + z_max/z_res;     %[no unit] number of depth layers
depths = linspace(0, z_max, ndepths); %[m] depth
z_res = ones(size(depths))*z_res; %[m] depth resolution

%% Definition of the temporal domain
stoptime = 20000;       %[a] total timespan of the problem
interval=1/32000;          %[a] time resolution (1/60000 is nine minutes, 1/8760 is one hour; 1/365.2 is a day)
t_length=stoptime/interval;      %[no unit] number of time layers

%% Bottom-water environmental conditions
T=2.85;      %[C] temperature
S=34.6;   %[psu] salinity
P=1631;    %[dbar] pressure
rho_sw = gsw_rho(S,T,P);    %[kg/m^3] in situ seawater density computed from GSW toolbox

%% Bottom-water values of dissolved species
dO2w=(66)*1e-6*rho_sw; %[mol/m3] dissolved oxygen from CTD cast at sta3
dtalkw=(2400)*1e-6*rho_sw; %[mol/m3] dissolved oxygen from GLODAP at station location, bottom waters
dtCO2w=(2353)*1e-6*rho_sw; %[mol/m3] DIC from GLODAP at sation location, bottom waters
dt12CO2w=(2327.16843)*1e-6*rho_sw; %[mol/m3] DIC from GLODAP at sation location, bottom waters
dt13CO2w=(25.8315696)*1e-6*rho_sw; %[mol/m3] DIC from GLODAP at sation location, bottom waters
dtNO3w=(40)*1e-6*rho_sw; %[mol/m3] nitrate from GLODAP at sation location, bottom waters
dtSO4w=(29264.2*S/35)*1e-6*rho_sw; %[mol/m3] computer from salinity (Millero, 2013)
dtPO4w=(2.5)*1e-6*rho_sw; %[mol/m3] phosphate from GLODAP at sation location, bottom waters
dtNH4w=1e-6*rho_sw; %[mol/m3] assumed
dtH2Sw=0*1e-6*rho_sw; %[mol/m3] assumed
dFew=0.5*1e-9*rho_sw; %[mol/m3] typical for deep sea oxic bottom waters (Abadie et al., 2019)
dMnw=0.5*1e-9*rho_sw; %[mol/m3] typical for deep sea oxic bottom waters (Morton et al., 2019)
dtSiw=120*1e-6*rho_sw;  %[mol/m3] dissolved inorganic silica
dCaw=0.02163./40.087.*(S./1.80655)*rho_sw;  %[mol/m3] Ca, computed from salinity using Riley CG(1967)

%% depth-dependent porosity
phiBeta = 33;   %porosity attenuation coefficient
phiInf = 0.75;   %porosity at infinite depth
phi0 = 0.83;    %porosity at interface
phi = (phi0 - phiInf)*exp(-phiBeta*depths) + phiInf;   %porosity profile 
phiS=1-phi;   %solid volume fraction
tort=(1-2*log(phi)).^0.5;   %tortuosity from Boudreau (1996, GCA)
tort2=tort.^2;   %tortuosity squared

%% Redfield ratios
RC=106;     %Redfield ratio for C
RN=16;     %Redfield ratio for N
RP=1;    %Redfield ratio for P 
M_CH2O=30.031; %[g per mol] molar mass of CH2O
M_NH3=17.031; %[g per mol] molar mass of NH3
M_H3PO4=97.994; %[g per mol] molar mass of H3PO4
M_OM=M_CH2O+(RN/RC)*M_NH3+(RP/RC)*M_H3PO4; %[g of OM per mol of OC] Organic Matter molar mass

%% solid fluxes and solid initial conditions
Foc=0.12; %[mol/m2/a] flux of total organic carbon to the bottom 
Fo13c=1.2861e-3; %[mol/m2/a] flux of total organic carbon to the bottom 
Fo12c=118.7139e-3; %[mol/m2/a] flux of total organic carbon to the bottom 
Froc=Foc*0.03; %[mol/m2/a] flux of refractory organic carbon to the bottom 
Fro13c=Fo13c*0.03; %[mol/m2/a] flux of refractory organic carbon to the bottom 
Fro12c=Fo12c*0.03; %[mol/m2/a] flux of refractory organic carbon to the bottom 
Fsoc=Foc*0.27; %[mol/m2/a] flux of slow decay organic carbon to the bottom 
Fso13c=Fo13c*0.27; %[mol/m2/a] flux of slow decay organic carbon to the bottom 
Fso12c=Fo12c*0.27; %[mol/m2/a] flux of slow decay organic carbon to the bottom 
Ffoc=Foc*0.7; %[mol/m2/a] flux of fast decay organic carbon to the bottom 
Ffo13c=Fo13c*0.7; %[mol/m2/a] flux of fast decay organic carbon to the bottom 
Ffo12c=Fo12c*0.7; %[mol/m2/a] flux of fast decay organic carbon to the bottom 
FMnO2=0.00015; %typical for deep sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
FFeOH3=0.00015; %typical for deep sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
Fcalcite=0.12; %[mol/m2/a] flux of calcite to the seafloor 
F13calcite=1.31738e-3; %[mol/m2/a] flux of calcite to the seafloor 
F12calcite=118.68262e-3; %[mol/m2/a] flux of calcite to the seafloor 
Faragonite=0; %[mol/m2/a] flux of aragonite to the seafloor
Fclay=1/360.31; %[mol/m2/a] flux of clay to the bottom: 360.31 is the molar mass of montmorillonite, typical deep sea clay
%flux clay is computed to have a w value of 1.4 mm/a, corresponding to 14C ages

    fraction13foc=Ffo13c./(Ffo12c+Ffo13c);
    fraction12foc=Ffo12c./(Ffo12c+Ffo13c);
    fraction13soc=Fso13c./(Fso12c+Fso13c);
    fraction12soc=Fso12c./(Fso12c+Fso13c);
    fraction13roc=Fro13c./(Fro12c+Fro13c);
    fraction12roc=Fro12c./(Fro12c+Fro13c);
    fraction13calcite=F13calcite./(F12calcite+F13calcite);
    fraction12calcite=F12calcite./(F12calcite+F13calcite);

Ftot=Foc*M_OM+FMnO2*86.9368+FFeOH3*106.867+Fcalcite*100.0869+Faragonite*100.0869+Fclay*360.31; %[g/m2/a] total sediment flux 
v0=(Ftot)/(2.65e6*phiS(1));                                             %[m/a] bulk burial velocity at sediment-water interface
vinf=v0*phiS(1)/phiS(1,ndepths);                                    %[m/a]bulk burial velocity at the infinite depth
for j=1:ndepths
    u(1,j)=vinf*phi(1,ndepths)/phi(1,j);                               %[m/a] porewater burial velocity
    w(1,j)=vinf*phiS(1,ndepths)/phiS(1,j);                         %[m/a] solid burial velocity
end

%% diffusive boundary layer
dbl=1e-3;            %[m] thickness at location taken from Sulpis et al 2018 PNAS

rerun = 2;
load("ini.mat")
time_saved_resolution=1; %[a]