function total = RADI_top(VAR, D_VAR, VARw, TotR_slice, TR, alpha, u, ...
    tort2, t_res, z_res)
diffusion = t_res*(D_VAR(1)/tort2(1)*(2*VAR(2) - 2*VAR(1) + ...
    TR(1)*(VARw - VAR(1)))/z_res(1)^2);
advection = u(1)*TR(1)*(VARw - VAR(1))/(2*z_res(1));
irrigation = alpha(1)*(VARw - VAR(1));
reaction = TotR_slice(1);
total = VAR(1) + diffusion + advection + irrigation + reaction;
