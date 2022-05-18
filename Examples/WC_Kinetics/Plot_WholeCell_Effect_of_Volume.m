% function plot_WholeCell_Effect_of_Volume
%
%   150mM KCl in pipette and bath, 100mV holding.
%   step increase in membrane conductance to g_K+ = 20nS at t = 0 to t = 20s
%   3 difference cell volumes, 4, 1, 0.25 pL.
%   Plot I, rho, versus t.

function plot_WholeCell_Effect_of_Volume

	% Configuration parameters
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	% Solution Parameters
        params.D =  [1,  1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ and Cl-.
		params.z =  [1, -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [150,  150  ]';  	% 150mM NaCl in bath
		params.rho_pip =  [150,  150  ]';  	% 150mM KCl in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
        params.theta = 10 * pi/180;  % Tip Angle
        [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
		params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
        params.Nseg  =  11; 		% Number of segments between the pipette electrode and tip
		params.V_hold = 100e-3; 	% Initial holding voltage
        params.Rseries = 4.5e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Cell and Membrane parameters
        %params.membrane_x = 1e-9;   			% Position of membrane patch back from the tip (in m)
		params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0]'; 		% Membrane permability to each ion in m/s
        params.membrane_Pf = 0e-5   ; 			% Membrane permeability to water (m/s)
		params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell
        
	% Initialize all remaining parameters
    params = initialize_system(params);
    
    % Calculate membrane permeability for an initial current of 2nA K+
    I_mem = [2e-9, 0 ]'; % Initial current of 2nA by K+
    P_GHK = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, params.V_hold, params.z, params.A_cell, params.T);
    
    % Report Access Resistance and Series Resistance Correction
	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params)/1e6), params.Rseries/1e6)
    fprintf('\n Membrane Permeability = [%e, %e] m/s', P_GHK);
    
    % Create a protocol
    % 20 second activation of 2nA current
    %perm = 3.5035e-8 * [1 0]'; % Membrane permeability
    perm = P_GHK; % Membrane permeability
    intervals.duration =        [5,     1e-3,  0.999,  4,   15,   1e-3   0.999, 4] ; % Duration of each of intervals
    intervals.N_eval =          [2,     10,     10,    10,  10,   10,    10 , 10] ; % Number of timepoints
    intervals.variables.membrane_perm = { 0*perm, perm,   perm, perm,   perm,   0*perm, 0*perm, 0*perm}; % Membrane permeability
    protocol = create_protocol( intervals );
       
    
    % Create volumes
    lambda = [0.25e-15, 1e-15, 4e-15] ;  
    
    figure(1); clf;
    
    for j = 1:length(lambda)
    
        % Initialize state and volume
        paramst = params;
        paramst.Lambda(paramst.Nseg-1) = lambda(j); % set cell volume
        [paramst, results] = evolve_backward_euler(paramst, inf);
        [paramst, results] = evolve_backward_euler(paramst, protocol);
        
        % Current versus time
        subplot(2,1,1); plot(protocol.t, results.I_pip * 1e9,'k'); hold on;
        subplot(2,1,2); plot(protocol.t, results.rho_intra(1,:),'k' ); hold on;
    end
        
    % Now label plots
        subplot(2,1,1); axis off;
        xscale = [0 30];
        axis([xscale 0 2]); 
        plot([ 10 20], [0.5 0.5],'k', 'Linewidth',1); % 10s time bar
        plot( [2.5, 2.5], [0.5 1.5], 'k', 'Linewidth',1); % 1nA current scale bar
        
        subplot(2,1,2);
        box off
        axis([xscale 85, 155]);
        ylabel('\rho (mM)');
        set(gca,'YTick',[110, 130 150],'TickDir','out');
        
        save_figure('Fig_WholeCell_Effect_of_Volume',[ 6 10]);   
     