%% top boundary condition: prescribed solid fluxes and diffusive boundary layer control on solutes
TA(1) = TA(1) + t_res * ( D_TA(1) / tort2(1) * (2*TA(2) - 2*TA(1) + TR(1) * (TAw - TA(1))) / (z_res(1).^2) ... %diffusion
    - u(1) *  -1 * TR(1) * (TAw - TA(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * (TAw - TA(1)) ... %irrigation
    + TotR(1,1)); %reaction

DIC(1) = DIC(1) + t_res * ( D_DIC(1) / tort2(1) * (2*DIC(2) - 2*DIC(1) + TR(1) * (DICw - DIC(1))) / (z_res(1).^2) ... %diffusion
    - u(1) *  -1 * TR(1) * (DICw - DIC(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * (DICw - DIC(1)) ... %irrigation
    + TotR(2,1)); %reaction

Calcite(1) = Calcite(1) + t_res * (D_bio(1) * ( 2 * Calcite(2) - 2 * Calcite(1) +... %diffusion
    2 * z_res(1) * (Fcalcite - phiS(1) * w(1) * Calcite(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
    - w(1) * -1 * (Fcalcite - phiS(1) * w(1) * Calcite(1)) / (D_bio(1) * phiS(1))... %advection
    + TotR(3,1)); %reaction

Aragonite(1) = Aragonite(1) + t_res * (D_bio(1) * ( 2 * Aragonite(2) - 2 * Aragonite(1) +... %diffusion
    2 * z_res(1) * (Faragonite - phiS(1) * w(1) * Aragonite(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
    - w(1) * -1 * (Faragonite - phiS(1) * w(1) * Aragonite(1)) / (D_bio(1) * phiS(1))... %advection
    + TotR(4,1)); %reaction

O2(1) = O2(1) + t_res * ( D_O2 / tort2(1) * (2*O2(2) - 2*O2(1) + TR(1) * (O2w - O2(1))) / (z_res(1)^2) ... %diffusion
    - u(1) * -1 * TR(1) * ( O2w - O2(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * ( O2w - O2(1) ) ... %irrigation
    + TotR(5,1)); %reaction

OC_labile(1) = OC_labile(1) + t_res * (D_bio(1) * ( 2 * OC_labile(2) - 2 * OC_labile(1) +... %diffusion
    2 * z_res(1) * (Flabile - phiS(1) * w(1) * OC_labile(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
    -w(1) * -1 * (Flabile - phiS(1) * w(1) * OC_labile(1)) / (D_bio(1) * phiS(1))... %advection
    + TotR(6,1)); %reaction

OC_refractory(1) = OC_refractory(1) + t_res * (D_bio(1) * ( 2 * OC_refractory(2) - 2 * OC_refractory(1) +... %diffusion
    2 * z_res(1) * (2*Frefractory - phiS(1) * w(1) * OC_refractory(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
    -w(1) * -1 * (2*Frefractory - phiS(1) * w(1) * OC_refractory(1)) / (D_bio(1) * phiS(1))... %advection
    +TotR(7,1)); %reaction

Ca(1) = Ca(1) + t_res * ( D_Ca(1) / tort2(1) * (2*Ca(2) - 2*Ca(1) + TR(1) * (Caw - Ca(1))) / (z_res(1).^2) ... %diffusion
    - u(1) *  -1 * TR(1) * (Caw - Ca(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * (Caw - Ca(1)) ... %irrigation
    + TotR(8,1)); %reaction

MnO2(1) = MnO2(1) + t_res * (D_bio(1) * ( 2 * MnO2(2) - 2 * MnO2(1) +... %diffusion
    2 * z_res(1) * (Fmno2 - phiS(1) * w(1) * MnO2(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
    - w(1) * -1 * (Fmno2 - phiS(1) * w(1) * MnO2(1)) / (D_bio(1) * phiS(1))... %advection
    + TotR(9,1)); %reaction                                                                                  

FeOH3(1) = FeOH3(1) + t_res * (D_bio(1) * ( 2 * FeOH3(2) - 2 * FeOH3(1) +... %diffusion
    2 * z_res(1) * (Ffeoh3 - phiS(1) * w(1) * FeOH3(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
    - w(1) * -1 * (Ffeoh3 - phiS(1) * w(1) * FeOH3(1)) / (D_bio(1) * phiS(1))... %advection
    + TotR(10,1)); %reaction

NO3(1) = NO3(1) + t_res * ( D_NO3(1) / tort2(1) * (2*NO3(2) - 2*NO3(1) + TR(1) * (NO3w - NO3(1))) / (z_res(1).^2) ... %diffusion
    - u(1) *  -1 * TR(1) * (NO3w - NO3(1)) / (2*z_res(1)) ...  %advection
    + alpha(1) * (NO3w - NO3(1))... %irrigation
    + TotR(11,1));  %reaction

SO4(1) = SO4(1) + t_res * ( D_SO4(1) / tort2(1) * (2*SO4(2) - 2*SO4(1) + TR(1) * (SO4w - SO4(1))) / (z_res(1).^2) ...  %diffusion
    - u(1) *  -1 * TR(1) * (SO4w - SO4(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * (SO4w - SO4(1))... %irrigation
    + TotR(12,1)); %reaction

PO4(1) = PO4(1) + t_res * ( (phi(1) * D_PO4 / tort2(1) + phiS(1) * D_bio(1) * KPA) * ( 2 * PO4(2) - 2 * PO4(1)... %diffusion
    + TR(1) * (PO4w - PO4(1) )) / (z_res(1).^2) ... %diffusion
    - (phi(1) * u(1) + phiS(1) * w(1) * KPA - phi(1) * D_PO4 * DFF(1) - KPA * DBF(1)) * -1 * TR(1)... %advection
    * (PO4w - PO4(1)) / (2*z_res(1)) ... %advection
    + phi(1) * alpha(1) * (PO4w - PO4(1)) ... %irrigation
    +TotR(13,1))./ (phi(z_length)+phiS(z_length)*KPA); %reaction

NH4(1) = NH4(1) + t_res * ( (phi(1) * D_NH4 / tort2(1) + phiS(1) * D_bio(1) * KNA) * ( 2 * NH4(2) - 2 * NH4(1)... %diffusion
    + TR(1) * (NH4w - NH4(1) )) / (z_res(1).^2) ...  %diffusion
    - (phi(1) * u(1) + phiS(1) * w(1) * KNA - phi(1) * D_NH4 * DFF(1) - KNA * DBF(1)) * -1 * TR(1)... %advection
    * (NH4w - NH4(1)) / (2*z_res(1)) ... %advection
    + phi(1) * alpha(1) * (NH4w - NH4(1)) ...  %irrigation
    + TotR(14,1)) / (phi(z_length)+phiS(z_length)*KNA); %reaction

H2S(1) = H2S(1) + t_res * ( D_H2S(1) / tort2(1) * (2*H2S(2) - 2*H2S(1) + TR(1) * (H2Sw - H2S(1))) / (z_res(1).^2) ...  %diffusion
    - u(1) *  -1 * TR(1) * (H2Sw - H2S(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * (H2Sw - H2S(1)) ... %irrigation
    + ( phiS(1) / phi(1) ) * (OC_labile(1) * klabile(1) + OC_refractory(1) * krefractory(1) ) * fSO4(1) / 2 ...  %reaction
    +TotR(15,1)); %reaction

Mn2p(1) = Mn2p(1) + t_res * ( D_Mn2p(1) / tort2(1) * (2*Mn2p(2) - 2*Mn2p(1) + TR(1) * (Mn2pw - Mn2p(1))) / (z_res(1).^2) ...  %diffusion
    - u(1) *  -1 * TR(1) * (Mn2pw - Mn2p(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * (Mn2pw - Mn2p(1)) ...   %irrigation
    + ( phiS(1) / phi(1) ) * 2 * (OC_labile(1) * klabile(1) + OC_refractory(1) * krefractory(1) ) * fMnO2(1)... %reaction
    +TotR(16,1)); %reaction

Fe2p(1) = Fe2p(1) + t_res * ( D_Fe2p(1) / tort2(1) * (2*Fe2p(2) - 2*Fe2p(1) + TR(1) * (Fe2pw - Fe2p(1))) / (z_res(1).^2) ... %diffusion
    - u(1) *  -1 * TR(1) * (Fe2pw - Fe2p(1)) / (2*z_res(1)) ... %advection
    + alpha(1) * (Fe2pw - Fe2p(1)) ... %irrigation
    + TotR(17,1)); %reaction

Clay(1) = Clay(1) + t_res * (D_bio(1) * ( 2 * Clay(2) - 2 * Clay(1) +... %diffusion
    2 * z_res(1) * (Fclay - phiS(1) * w(1) * Clay(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ... %diffusion
    - w(1) * -1 * (Fclay - phiS(1) * w(1) * Clay(1)) / (D_bio(1) * phiS(1))); %advection
