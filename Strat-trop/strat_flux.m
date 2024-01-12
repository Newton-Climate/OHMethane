% CH4 concentrations (you might adjust these based on data or starting assumptions)

load('case5_H2O_strat_real_inversion.mat')
strat_conc = out_anal.nh_ch4_strat;  % ppb (for illustration)
trop_conc = out_anal.nh_ch4;  % ppb (for illustration)

% Exchange timescale (you'll have to make an informed estimate or find data for this)
tau_exchange = ems_anal(:,15); % year (just a placeholder value)

% Lifetime of methane in the stratosphere
tau_lifetime_strat = 7; % years

% Simulation time step
dt = 0.01; % year (for example, representing about 3.65 days)

% Stratosphere values
P_strat = 100e2; % converting hPa to Pa
T_strat = 253; % K
% Volume calculation for simplicity - this is an illustrative value
radius_earth = 6371e3; % m
height_strat = 8e3; % considering stratosphere from 17km to 25km
V_strat = pi * height_strat * (radius_earth^2 + (radius_earth + height_strat)^2 + radius_earth * (radius_earth + height_strat));
V_strat = V_strat/2;

% Troposphere values
P_trop = 500e2; % converting hPa to Pa
T_trop = 248; % K
% Volume calculation for simplicity - this is an illustrative value
height_trop = 17e3; % 0 to 17km
V_trop = pi * height_trop * (radius_earth^2 + (radius_earth + height_trop)^2 + radius_earth * (radius_earth + height_trop));
V_trop = V_trop/2;
% Call the function
nh_flux = methane_flux(strat_conc, trop_conc, tau_exchange, P_strat, T_strat, V_strat, P_trop, T_trop, V_trop);

sh_flux = methane_flux(out_anal.sh_ch4_strat, out_anal.sh_ch4, tau_exchange, P_strat, T_strat, V_strat, P_trop, T_trop, V_trop);


% Display the results
%disp(['Flux: ', num2str(flux), ' Tg/yr']);
%disp(['New Stratosphere CH4 Concentration: ', num2str(new_strat_conc), ' ppb']);



function flux = methane_flux(strat_conc, trop_conc, tau_exchange, P_strat, T_strat, V_strat, P_trop, T_trop, V_trop)
    % Constants
    R = 8.314;                 % J/(mol·K)
    molar_mass_CH4 = 16.04e-12;  % Tg/mol

    % Convert ppb to mole fraction
    strat_mole_frac = strat_conc * 1e-9;
    trop_mole_frac = trop_conc * 1e-9;

    % Calculate moles of CH4 in each box
    moles_CH4_strat = (P_strat .* V_strat) ./ (R * T_strat) .* strat_mole_frac;
    moles_CH4_trop = (P_trop .* V_trop) ./ (R * T_trop) .* trop_mole_frac;

    % Calculate mass of CH4 in each box
    mass_CH4_strat = moles_CH4_strat * molar_mass_CH4;
    mass_CH4_trop = moles_CH4_trop * molar_mass_CH4;

    % Compute flux in Tg/yr
    flux = (mass_CH4_strat - mass_CH4_trop) ./ tau_exchange;

    % Update stratosphere CH4 concentration based on first-order decay
%    decay_strat = -strat_conc / tau_lifetime_strat * dt;  % Using the first-order rate equation
%    new_strat_conc = strat_conc + decay_strat;  % Update the concentration

end
