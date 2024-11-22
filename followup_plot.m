%% compute solid fraction 
load('data_st3.mat');
Xca=100*(pcalcitef.*100.09)./(pclayf.*258.17+pcalcitef.*100.09+paragonitef.*100.09+procf.*12+pfocf.*12+psocf.*12+pFeOH3f.*106.867);
Xoc=30*(pfocf.*12+psocf.*12+procf.*12)./(pclayf.*258.17+pcalcitef.*100.09+paragonitef.*100.09+procf.*12+pfocf.*12+psocf.*12+pFeOH3f.*106.867);
Xroc=30*(procf.*12)./(pclayf.*258.17+pcalcitef.*100.09+paragonitef.*100.09+procf.*12+pfocf.*12+psocf.*12+pFeOH3f.*106.867);
Xsoc=30*(psocf.*12)./(pclayf.*258.17+pcalcitef.*100.09+paragonitef.*100.09+procf.*12+pfocf.*12+psocf.*12+pFeOH3f.*106.867);
Xfoc=30*(pfocf.*12)./(pclayf.*258.17+pcalcitef.*100.09+paragonitef.*100.09+procf.*12+pfocf.*12+psocf.*12+pFeOH3f.*106.867);

%depth-time color plots

figure(1)

subplot 331
pcolor([1:10:idx.*10-10],-depths,dtCO2f(:,1:idx-1)./rho_sw.*1e6)
title('DIC (umol/kg)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat

subplot 332
pcolor([1:10:idx.*10-10],-depths,dtalkf(:,1:idx-1)./rho_sw.*1e6)
title('TAlk (umol/kg)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat

subplot 333
pcolor([1:10:idx.*10-10],-depths,dO2f(:,1:idx-1)./rho_sw.*1e6)
title('oxygen (umol/kg)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat

subplot 334
pcolor([1:10:idx.*10-10],-depths,dtNO3f(:,1:idx-1)./rho_sw.*1e6)
title('nitrate (umol/kg)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat

subplot 335
pcolor([1:10:idx.*10-10],-depths,Xca(:,1:idx-1))
title('calcite (dry wt %)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat

subplot 336
pcolor([1:10:idx.*10-10],-depths,Xoc(:,1:idx-1))
title('poc (dry wt %)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat

subplot 337
pcolor([1:10:idx.*10-10],-depths,dFef(:,1:idx-1)./rho_sw.*1e6)
title('Fe (umol/kg)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat

subplot 338
pcolor([1:10:idx.*10-10],-depths,dMnf(:,1:idx-1)./rho_sw.*1e6)
title('Mn (umol/kg)')
ylabel('depth (m)')
xlabel('time (a)')
colorbar
shading flat



%% depth profiles

figure(2)
clf

subplot 331
plot(dO2f(:,idx-1)./rho_sw.*1e6,-depths.*100)
title('O2')
ylabel('depth (cm)')
xlabel('concentration (umol/kg)')

subplot 332
plot(dtNO3f(:,idx-1)./rho_sw.*1e6,-depths.*100)
title('NO3')
hold on 
ylabel('depth (cm)')
xlabel('concentration (umol/kg)')
scatter(data_st3_NO3(:,2),-data_st3_NO3(:,1),'filled')

subplot 333
plot(dtalkf(:,idx-1)./rho_sw.*1e6,-depths.*100)
title('TAlk')
hold on 
ylabel('depth (cm)')
xlabel('concentration (umol/kg)')
scatter(data_st3_TAlk(:,2),-data_st3_TAlk(:,1),'filled')

subplot 334
plot(dtCO2f(:,idx-1)./rho_sw.*1e6,-depths.*100)
title('DIC')
hold on 
ylabel('depth (cm)')
xlabel('concentration (umol/kg)')
scatter(data_st3_DIC(:,2),-data_st3_DIC(:,1),'filled')

subplot 335
plot(dtSO4f(:,idx-1)./rho_sw.*1e3,-depths.*100)
title('SO4')
ylabel('depth (cm)')
xlabel('concentration (mmol/kg)')


subplot 336
plot(dFef(:,idx-1)./rho_sw.*1e6,-depths.*100)
title('Fe')
ylabel('depth (cm)')
xlabel('concentration (umol/kg)')

subplot 337
plot(Xoc(:,idx-1),-depths.*100)
title('X poc')
hold on 
scatter(data_st3_POC(:,1),-data_st3_PIC(:,1),'filled')
ylabel('depth (cm)')
xlabel('concentration (dry wt %)')

subplot 338
plot(Xca(:,idx-1),-depths.*100)
title('X calcite')
hold on 
scatter(data_st3_PIC(:,2),-data_st3_PIC(:,1),'filled')
ylabel('depth (cm)')
xlabel('concentration (dry wt %)')

%%
figure(3)
clf

subplot 231
plot(dCaf(:,idx-1)./rho_sw.*1e3,-depths.*100)
title('Ca')
ylabel('depth (cm)')
xlabel('concentration (mmol/kg)')
hold on 
scatter(data_st3_Ca(:,2),-data_st3_Ca(:,1),'filled')

subplot 232
plot(dMnf(:,idx-1)./rho_sw.*1e6,-depths.*100)
title('Mn')
hold on 
scatter(data_st3_Mn(:,2),-data_st3_Mn(:,1),'filled')
ylabel('depth (cm)')
xlabel('concentration (umol/kg)')

subplot 233
plot(OmegaC,-depths.*100)
title('Omega')
hold on 
scatter(data_st3_OmegaC(:,2),-data_st3_OmegaC(:,1),'filled')
ylabel('depth (cm)')
xlabel('with respect to calcite')

subplot 234
plot(-log10(H),-depths.*100)
title('pH')
hold on 
scatter(data_st3_pH(:,2),-data_st3_pH(:,1),'filled')
ylabel('depth (cm)')
xlabel('pH total scale')

subplot 235
plot(dFef(:,idx-1)./rho_sw.*1e6,-depths.*100)
title('Fe')
ylabel('depth (cm)')
xlabel('concentration (umol/kg)')

subplot 236
plot(fdO2,-depths.*100,'blue-','LineWidth',4)
hold on
plot(fdtNO3,-depths.*100,'green-','LineWidth',4)
plot(fpMnO2,-depths.*100,'red-','LineWidth',4)
plot(fpFeOH3,-depths.*100,'black-','LineWidth',4)
plot(fdtSO4,-depths.*100,'yellow-','LineWidth',4)
plot(fdCH4,-depths.*100)
title('degradation pathways')



%% d13C

delta13C_DIC=(((dt13CO2f./dt12CO2f)./0.0111)-1).*1000;
delta13C_PIC=(((p13calcitef./p12calcitef)./0.0111)-1).*1000;
delta13C_POC=((((pro13cf+pso13cf+pfo13cf)./(pro12cf+pso12cf+pfo12cf))./0.0111)-1).*1000;

figure(4)

subplot 311
pcolor([1:1:idx-1],-depths,delta13C_DIC(:,1:idx-1))
title('delta13C DIC')
colorbar
shading flat
ylabel('depth (m)')
xlabel('time (a)')

subplot 312
pcolor([1:1:idx-1],-depths,delta13C_PIC(:,1:idx-1))
title('delta13C PIC')
colorbar
shading flat
ylabel('depth (m)')
xlabel('time (a)')

subplot 313
pcolor([1:1:idx-1],-depths,delta13C_POC(:,1:idx-1))
title('delta13C POC')
colorbar
shading flat
ylabel('depth (m)')
xlabel('time (a)')

%% d13C plots along with data 

figure(5)
clf

subplot 311
plot(delta13C_DIC(:,idx-1),-depths.*100)
title('d13C DIC')
hold on 
scatter(data_st3_DIC_d13C(:,2),-data_st3_DIC_d13C(:,1),'filled')
ylabel('permil')
xlabel('depth (cm)')

subplot 312
plot(delta13C_PIC(:,idx-1),-depths.*100)
title('d13C PIC')
hold on 
scatter(data_st3_PIC_d13C(:,2),-data_st3_PIC_d13C(:,1),'filled')
ylabel('permil')
xlabel('depth (cm)')

subplot 313
plot(delta13C_POC(:,idx-1),-depths.*100)
title('d13C POC')
hold on 
scatter(data_st3_POC_d13C(:,2),-data_st3_POC_d13C(:,1),'filled')
ylabel('permil')
xlabel('depth (cm)')

%% rates 

figure(6)
clf

subplot 311
plot(Rd_calcite,-depths.*100)
title('Calcite dissolution rate')
hold on 
ylabel('depth (cm)')
xlabel('rate (mol/m3/a)')

subplot 312
plot(Rp_calcite,-depths.*100)
title('Calcite precipitation rate')
hold on 
ylabel('depth (cm)')
xlabel('rate (mol/m3/a)')

subplot 313
plot(log10(Rf_tot),-depths.*100)
hold on 
plot(log10(Rs_tot),-depths.*100)
title('Fast and slow-decay POC degradation rate')
ylabel('depth (cm)')
xlabel('rate (mol/m3/a)')







