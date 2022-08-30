function mu_0 = VicFc_Get_elastoviscousNum(L, u_0, a)
%% Function to calculate the mu_0
% 
% This is to calculate the elasto-viscous number
% a: pillar diameter

B = 6.9e-26;  % Bending rigidity
mu = 6.1e-3;  % Dynalic viscosity
d = 8e-9; % Diameter of the actin filament

mu_0 = 8 * pi * mu * L^4 * u_0 / (B * a * -log((d/L)^2 * exp(1)));

end