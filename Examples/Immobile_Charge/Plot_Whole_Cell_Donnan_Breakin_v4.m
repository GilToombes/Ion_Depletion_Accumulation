% function Fig_Whole_Cell_Donnan
%
%   Plot concentration and voltage-chang at the end of pipette containing
%   150mM KCl as a function of IK+*Raccess, and |INa+| * Raccess


function Fig_Whole_Cell_Donnan_Breakin

	% Configuration parameters
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	% Solution Parameters
        params.D =  [1, 0.682, 1.0382,  0.33,  0]' * 1.957e-9; 	% Diffusion coefficients of K+ , Na+, Cl-, gluconate and immobile ions
		params.z =  [1, 1, -1, -1, -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [15, 135, 150,  0, 0  ]';  	% 150mM KCl in bath
		params.rho_pip =  [135, 15, 10, 140, 0  ]';  	% 150mM KCl in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
        params.theta = 10 * pi/180;  % Tip Angle
        [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
        params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
    	params.Nseg  =  21; 		% Number of segments between the pipette electrode and tip
		params.V_hold = -80e-3; 	% Initial holding voltage
        params.Rseries = 0e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Cell and Membrane parameters
        params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0 0 0 0 ]'; 		% Membrane permability to each ion in m/s
        params.membrane_Pf = 0e-5   ; 			% Membrane permeability to water (m/s)
		params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell

    % Initialize all remaining parameters
    params = initialize_system(params);
    
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
   	fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params)/1e6), sum(params.Rseries/1e6))

    % Exclude the immobile anions from concentration changes and fluxes.  
    params.Constraints.components = {[1 4], [1 4], [1 1]};
    params = initialize_system(params);
    
    % Calculate rho/V versus immobile anions 
    % Set up immobile anions
    rho_im = 100; % Concentration of immobile anions in cell (mM)
    d_patch = 3e-6; % Concentration ramps linearly starting from d_pip back in pipette tip
    [tmp, N_patch ] = min(abs(params.x + d_patch)); 
    
    for j = N_patch:params.membrane_N
            r = (j-N_patch)/(params.membrane_N -N_patch) * rho_im;
            params.rho(:,j) = [params.rho(1,j), params.rho(2,j), params.rho(3,j), params.rho(4,j)-r,  r]' ;
    end
    
    % Now need to initialize before breaking patch.
    params.patch_func = @patch_current_GHK;
    params.patch_perm = [ 10, 0 , 0, 0, 0 ]' * 1e-5; 		% Patch permability to each ion in m/s
    params.patch_N = N_patch;
    params = evolve_backward_euler(params,0.1);
    params = evolve_backward_euler(params,10);
    params.t = -8;
    
    t_val = [-8 0, [0:1e-4:0.9e-3, 0.001:0.001:0.009, 0.01:0.01:0.09, 0.1:0.1:0.9, 1:30]];
    
    for k=1:length(t_val)
        
        [params, results] = evolve_backward_euler(params, t_val(k));
        if (t_val(k) >=0)&isfield(params,'patch_func') 
            % Remove patch
            params = rmfield(params,'patch_func');
        end
    
        rho(:,k) = results.rho_intra;
        V(k) = (results.V_mem - params.V_hold) ;
        J(:,k) = params.J(:,params.membrane_N-1) ;
        
    end
    
    
    left = 0.15;
    width = 0.75;
    
    
    plotpos = [left, 0.8, width, 0.18;
               left, 0.64,width, 0.12 ;
               left, 0.43, width, 0.18;
               left, 0.25,width, 0.18 ;
               left, 0.05,width, 0.18;];
    
    figure(1); clf;
        x_range = [-8 30];
    
        % Molecular Currents
        subplot('Position', plotpos(1,:));
             plot(t_val, J(1,:) * params.F * 1e9 ,'b'); hold on ; 
             plot(t_val, J(2,:) * params.F * 1e9,'y');
             plot(t_val, -1 * J(3,:) * params.F * 1e9,'c');
             plot(t_val, -1 * J(4,:) * params.F * 1e9,'r');
             I_val = (J(1,:) + J(2,:) - J(3,:) - J(4,:)) * params.F;
             plot(t_val, I_val * 1e9 ,'k'); hold on;
             
            axis([x_range -1.3 1.3]);
            axis off;
           plot( [10 20], [-1 -1],'k', 'Linewidth',1); % 10s scale bar
           plot( [-3 -3 ], 0.1 + [0 1], 'k', 'Linewidth', 1); % 1nA current bar
             
        
        % Net Current
        %subplot('Position', plotpos(2,:));
        %     axis([x_range -1.3 1.3]);
        %     axis off;
        %    plot( [10 20], [-1 -1],'k', 'Linewidth',1); % 10s scale bar
        %    plot( [-3 -3 ], 0.1 + [0 1], 'k', 'Linewidth', 1); % 1nA current bar
   
        % Voltage
        subplot('Position', plotpos(2,:));
            plot(t_val, V * 1e3,'k');
                axis([x_range -9 0]);
                box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
                ylabel('\DeltaV (mV)'); 
            
        % Concentration
        subplot('Position', plotpos(3,:));
            plot(t_val, rho(1,:),'b'); hold on ;
            plot(t_val, rho(4,:),'r'); hold on ;
            plot(t_val, rho(5,:),'k'); hold on ;
            axis([x_range, 30 220]);
             box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
             ylabel('\rho_{cell} (mM)');
             
        subplot('Position', plotpos(4,:));
            plot(t_val, rho(2,:),'y'); hold on ;
            plot(t_val, rho(3,:),'c'); hold on ;
                
            axis([x_range, 0 30]);
             box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
                
             
             
        % Chemical Potentials
        subplot('Position', plotpos(5,:));
            mu_K =  log(rho(1,:)/params.rho_pip(1))       + V / (params.R * params.T / params.F);
            mu_Na = log(rho(2,:)/params.rho_pip(2))       + V / (params.R * params.T / params.F);
            mu_Cl = log(rho(3,:)/params.rho_pip(3))       - V / (params.R * params.T / params.F);
            mu_Gluc=log(rho(4,:)/params.rho_pip(4))        - V / (params.R * params.T / params.F);
    
            plot(t_val, mu_K,'b'); hold on ;
            plot(t_val, mu_Na,'y'); hold on ;
            plot(t_val, mu_Cl,'c'); 
            plot(t_val, mu_Gluc,'r'); 
            axis([x_range, -1.3 0.5]);
             box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
            ylabel('\Delta\mu (k_BT)');
    
            
        save_figure('Fig_WC_Donnan_Breakin',[ 10 14]);         
        
        
        