%% this script saves the last resolved-year variables required to run the RADI model, to set them as initial conditions.
   
    
    dtalkic=dtalkf(:,idx-1);            %[mol/m3]
    dtCO2ic=dtCO2f(:,idx-1);            %[mol/m3]
    dO2ic=dO2f(:,idx-1);                %[mol/m3]
    dCaic=dCaf(:,idx-1);                %[mol/m3]
    dtNO3ic=dtNO3f(:,idx-1);            %[mol/m3]
    dtSO4ic=dtSO4f(:,idx-1);            %[mol/m3]
    dtPO4ic=dtPO4f(:,idx-1);            %[mol/m3]
    dtNH4ic=dtNH4f(:,idx-1);            %[mol/m3]
    dtH2Sic=dtH2Sf(:,idx-1);            %[mol/m3]
    dMnic=dMnf(:,idx-1);                %[mol/m3]
    dFeic=dFef(:,idx-1);                %[mol/m3]
    dt12CO2ic=dt12CO2f(:,idx-1);        %[mol/m3]
    dt13CO2ic=dt13CO2f(:,idx-1);        %[mol/m3]
    
    % initial condition for solids
    pcalciteic=pcalcitef(:,idx-1);      %[mol/m3]
    p12calciteic=p12calcitef(:,idx-1);  %[mol/m3]
    p13calciteic=p13calcitef(:,idx-1);  %[mol/m3]
    paragoniteic=paragonitef(:,idx-1);  %[mol/m3]
    pfocic=pfocf(:,idx-1);              %[mol/m3]
    pfo12cic=pfo12cf(:,idx-1);          %[mol/m3]
    pfo13cic=pfo13cf(:,idx-1);          %[mol/m3]
    psocic=psocf(:,idx-1);              %[mol/m3]
    pso12cic=pso12cf(:,idx-1);          %[mol/m3]
    pso13cic=pso13cf(:,idx-1);          %[mol/m3]
    pMnO2ic=pMnO2f(:,idx-1);            %[mol/m3]
    pFeOH3ic=pFeOH3f(:,idx-1);          %[mol/m3]
    pclayic=pclayf(:,idx-1);            %[mol/m3] 
    procic=procf(:,idx-1);              %[mol/m3]
    pro12cic=pro12cf(:,idx-1);              %[mol/m3]
    pro13cic=pro13cf(:,idx-1);              %[mol/m3]

    clearvars -except dtalkic dtCO2ic dO2ic dCaic dtNO3ic dtSO4ic dtPO4ic...
        dtNH4ic dtH2Sic dMnic dFeic dt12CO2ic dt13CO2ic pcalciteic p12calciteic...
        p13calciteic paragoniteic pfocic pfo12cic pfo13cic psocic pso12cic...
        pso13cic pMnO2ic pFeOH3ic pclayic procic pro12cic pro13cic
    
    save('ini.mat')