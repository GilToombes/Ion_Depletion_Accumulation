%function P_mem = calculate_PGHK_from_I(I_mem, rho_int, rho_ext, V_mem, z, A_mem, T)
%
%     Calculate GHK permeability required to give a membrane current, I_mem
%
%       Input -         I_mem -> membrane current carried by each species in Amps 
%                                   e.g. [1e-9, 0, 0]  
%                                   1nA Current carried by species 1
%                       rho_int -> intra-cellular concentration of each species in mM 
%                                   e.g. [150, 0, 150]
%                       rho_ext -> extra-cellular concentration of each species in mM 
%                                   e.g. [0, 150, 150]
%                       V_mem  ->  membrane voltage in Volts
%                                   e.g. V_mem = 1e-1; 100mV
%                       z      ->  unit charge of each species
%                                   e.g. [1 1 -1]
%                       A_mem  ->  membrane area in m^2
%                                   
%                       T      -> optional - temperature in Kelvin
%                                   default T = 293 ; 20 C
%
%       Outputs
%                       P_mem -> GHK membrane permeability in m/s
%                                is row vector

function P_mem = calculate_PGHK_from_I(I_mem, rho_int, rho_ext, V_mem, z, A_mem, T)

	if (nargin<7)
		T = 293; % Temperature in Kelvin
    end
	R = 8.314 ; % Gas constant in J/mol.K
    F =	96485 ; % Faraday's constant in Coulomb/mol

    u = z*F*V_mem/(R*T) ; 
	uz = (abs(u)<1e-6); % Dealing with Uncharged particles and zero membrane voltge
    Ius = A_mem * F  .* z .* (rho_int - rho_ext.*exp(-u)) .* (u + uz) ./ (1 - exp(-u) + uz) ; % Current in permeabilty was 1.
    
    P_mem = I_mem ./ Ius; % Membrane Permeability
    
    