%function [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip)
%
%     Calculate tip diameter for a given pipette_resistance
%
%       Input -         params.rho_pip -> ion concentrations in pipette (mM)
%                               e.g. [135 15 150] -> 150mM, 15mM, 150mM
%                       params.D       - diffusion coefficient
%                               e.g. [2e-9, 1.3e-9, 2e-9] in m^2/s
%                       params.z      unit charge
%                                   e.g. [1 1 -1]
%                       params.theta - tip angle in radians
%                       params.T     -  temperature in kelvin
%                       parmas.r_pip - inner radius of pipette in meters
%
%                       R_pip -> Resistance in Ohms
%
%       Outputs
%                       tip_diameter -> required tip diameter in meters
%                       sigma -> solution conductivity in S/m



function [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip)

    if ~isfield(params,'R'); params.R = 8.314; end; % Gas Constant in SI
	if ~isfield(params,'F'); params.F = 96485; end; % Faraday's constant
	if ~isfield(params,'T'); params.T = 293 ;  end; % Temperature in Kelvin.
    if ~isfield(params,'r_pip'); di = 0; else; di = 1 / (2 * params.r_pip); end; % Inverse of pipette inner diameter
    
    sigma = sum(params.rho_pip .* params.D .* params.z .* params.z) * params.F^2 / (params.R * params.T) ; % Conductivity in S/m
    alpha = R_pip * sigma * pi * tan(params.theta/4) ;
    tip_diameter =  1 / (alpha + di); % Tip Diameter in m
    