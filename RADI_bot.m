TA(z_length) = TA(z_length) + t_res * (D_TA / tort2(z_length) * 2 * ((TA(z_length-1) - TA(z_length)) / z_res(z_length).^2) ...   %diffusion
    +alpha(z_length) * (TAw - TA(z_length)) ... %irrigation
    +TotR(1,z_length));

DIC(z_length) = DIC(z_length) + t_res * (D_DIC / tort2(z_length) * 2 * ((DIC(z_length-1) - DIC(z_length)) / z_res(z_length).^2) ...   %diffusion
    +alpha(z_length) * (DICw - DIC(z_length)) ...  %irrigation
    +TotR(2,z_length));

Calcite(z_length) = Calcite(z_length) + t_res * (D_bio(z_length) * 2 * ( (Calcite(z_length-1) - Calcite(z_length))  / z_res(z_length).^2)...%diffusion
    - APPW(z_length) * (-sigma(z_length)*Calcite(z_length-1) + sigma(z_length)*Calcite(z_length))/z_res(z_length)...  %advection
    +TotR(3,z_length));

Aragonite(z_length) = Aragonite(z_length) + t_res * (D_bio(z_length) * 2 * ( (Aragonite(z_length-1) - Aragonite(z_length))  / z_res(z_length).^2)...%diffusion
    - APPW(z_length) * (-sigma(z_length)*Aragonite(z_length-1) + sigma(z_length)*Aragonite(z_length))/z_res(z_length)... %advection
    +TotR(4,z_length));

O2(z_length) = O2(z_length) + t_res * (D_O2 / tort2(z_length) * 2 * ((O2(z_length-1) - O2(z_length)) / z_res(z_length).^2) ...  %diffusion
    + alpha(z_length) * (O2w - O2(z_length)) ... %irrigation
    +TotR(5,z_length));

OC_labile(z_length) = OC_labile(z_length) + t_res * (D_bio(z_length) * 2 * ( (OC_labile(z_length-1) - OC_labile(z_length)) / z_res(z_length).^2)... %diffusion
    - APPW(z_length) * (-sigma(z_length)*OC_labile(z_length-1) + sigma(z_length)*OC_labile(z_length))/z_res(z_length)... %advection
    +TotR(6,z_length));

OC_refractory(z_length) = OC_refractory(z_length) + t_res * (D_bio(z_length) * 2 * ( (OC_refractory(z_length-1) - OC_refractory(z_length)) / z_res(z_length).^2)... %diffusion
    - APPW(z_length) * (-sigma(z_length)*OC_refractory(z_length-1) + sigma(z_length)*OC_refractory(z_length))/z_res(z_length)... %advection
    +TotR(7,z_length));

Ca(z_length) = Ca(z_length) + t_res * (D_Ca / tort2(z_length) * 2 * ((Ca(z_length-1) - Ca(z_length)) / z_res(z_length).^2) ... %diffusion
    +alpha(z_length) * (Caw - Ca(z_length)) ... %irrigation
    +TotR(8,z_length));

MnO2(z_length) = MnO2(z_length) + t_res * (D_bio(z_length) * 2 * ( (MnO2(z_length-1) - MnO2(z_length))  / z_res(z_length).^2)... %diffusion
    - APPW(z_length) * (-sigma(z_length) * MnO2(z_length-1) + sigma(z_length)*MnO2(z_length))/z_res(z_length)...  %advection
    +TotR(9,z_length));

FeOH3(z_length) = FeOH3(z_length) + t_res * (D_bio(z_length) * 2 * ( (FeOH3(z_length-1) - FeOH3(z_length))  / z_res(z_length).^2)...%diffusion
    - APPW(z_length) * (-sigma(z_length)*FeOH3(z_length-1) + sigma(z_length)*FeOH3(z_length))/z_res(z_length)...%advection
    +TotR(10,z_length));

NO3(z_length) = NO3(z_length) + t_res * (D_NO3 / tort2(z_length) * 2 * ((NO3(z_length-1) - NO3(z_length)) / z_res(z_length).^2) ...  %diffusion
    + alpha(z_length) * (NO3w - NO3(z_length)) ...  %irrigation
    +TotR(11,z_length));

SO4(z_length) = SO4(z_length) + t_res * (D_SO4 / tort2(z_length) * 2 * ((SO4(z_length-1) - SO4(z_length)) / z_res(z_length).^2) ...  %diffusion
    + alpha(z_length) * (SO4w - SO4(z_length)) ... %irrigation
    +TotR(12,z_length));

PO4(z_length) = PO4(z_length) + t_res * ( (phi(z_length) * D_PO4 / tort2(z_length) + phiS(z_length) * D_bio(z_length) * KPA) ... %diffusion
    * 2 * ((PO4(z_length-1) - PO4(z_length)) / z_res(z_length).^2) ...  %diffusion
    + phi(z_length) * alpha(z_length) * (PO4w - PO4(z_length)) ... %irrigation
    +TotR(13,z_length)) / (phi(z_length)+phiS(z_length)*KPA); %reaction

NH4(z_length) = NH4(z_length) + t_res * ( (phi(z_length) * D_NH4 / tort2(z_length) + phiS(z_length) * D_bio(z_length) * KNA) ... %diffusion
    * 2 * ((NH4(z_length-1) - NH4(z_length)) / z_res(z_length).^2) ...     %diffusion
    + phi(z_length) * alpha(z_length) * (NH4w - NH4(z_length)) ...  %irrigation
    +TotR(14,z_length)) / (phi(z_length)+phiS(z_length)*KNA);  %reaction

H2S(z_length) = H2S(z_length) + t_res * (D_H2S / tort2(z_length) * 2 * ((H2S(z_length-1) - H2S(z_length)) / z_res(z_length).^2) ...  %diffusion
    + alpha(z_length) * (H2Sw - H2S(z_length)) ... %irrigation
    +TotR(15,z_length));

Mn2p(z_length) = Mn2p(z_length) + t_res * (D_Mn2p / tort2(z_length) * 2 * ((Mn2p(z_length-1) - Mn2p(z_length)) / z_res(z_length).^2) ...  %diffusion
    + alpha(z_length) * (Mn2pw - Mn2p(z_length)) ... %irrigation
    +TotR(16,z_length));

Fe2p(z_length) = Fe2p(z_length) + t_res * (D_Fe2p / tort2(z_length) * 2 * ((Fe2p(z_length-1) - Fe2p(z_length)) / z_res(z_length).^2) ...  %diffusion
    + alpha(z_length) * (Fe2pw - Fe2p(z_length)) ...   %irrigation
    +TotR(17,z_length));

Clay(z_length) = Clay(z_length) + t_res * (D_bio(z_length) * 2 * ((Clay(z_length-1) - Clay(z_length))  / z_res(z_length).^2) ...%diffusion
    - APPW(z_length) * (-sigma(z_length)*Clay(z_length-1) + sigma(z_length)*Clay(z_length))/z_res(z_length));... %advection    