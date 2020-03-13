function [H] =  calc_pCO2(dic,pt,sit,bt,ta,ff,k1,k2,k1p,k2p,k3p,kb,kw,ksi,fg,H)

%Example FORTRAN subroutine: solve carbonate system for pC02
%M. Follows, T. Ito, S. Dutkiewicz
%D. Carroll, 2019

%routine variables
%real pt,sit,ta,pCO2,dic,H,ff,bt,k1,k2,k1p,k2p,k3p,kb,kw,ksi

%local variables
%real gamm,co2s,hg,cag,bohg,h3po4g,h2po4g,hpo4g,po4g, & siooh3g,denom,dummy,fg

%dic = dissolved inorganic carbon
%pt = dissolved inorganic phosphorus
%sit = dissolved inorganic silica
%bt = dissolved inorganic boron
%ta = total alkalinity 
%ca = carbonate alkalinity
%H = [H+]
%pCO2 = partial pressure CO2
%ff = fugacity of CO2

%k1, k2 = carbonate equilibrium coeffs; 
%kw = dissociation of water
%klp, k2p, k3p = phosphate equilibrium coefficients
%ksi, kb = silicate and borate equilibrium coefficients

%equilibrium relationships from DOE handbook (DOE, 1994), coefficients evaluated elsewhere and passed in.
%first guess of [H+]: from last timestep .*OR.* fixed for cold start hg = H
%estimate contributions to total alk from borate, silicate, phosphate

hg = H;

bohg = bt.*kb./(hg + kb);
siooh3g = sit.*ksi./(ksi + hg);
denom = hg.*hg.*hg + (k1p.*hg.*hg) + (k1p.*k2p.*hg) + (k1p.*k2p.*k3p);
h3po4g = (pt.*hg.*hg.*hg)./denom;
h2po4g = (pt.*k1p.*hg.*hg)./denom;
hpo4g = (pt.*k1p.*k2p.*hg)./denom;
po4g = (pt.*k1p.*k2p.*k3p)./denom;

%estimate carbonate alkalinity
fg = (-bohg - (kw./hg)) + hg - hpo4g - 2.0.*po4g + h3po4g - siooh3g;
cag = ta + fg;

%improved estimate of hydrogen ion conc
gamm = dic./cag;
%gamm = 1;
dummy = (1.0-gamm).*(1.0-gamm).*k1.*k1 - 4.0.*k1.*k2.*(1.0 - 2.0.*gamm);
H = 0.5.*((gamm-1.0).*k1 + sqrt(dummy));

%evaluate [CO2.*]
%co2s = dic / (1.0 + (k1/H) + (k1.*k2/(H.*H)));

%c evaluate surface pCO2
%pCO2 = co2s/ff; %mol kg^-1

end