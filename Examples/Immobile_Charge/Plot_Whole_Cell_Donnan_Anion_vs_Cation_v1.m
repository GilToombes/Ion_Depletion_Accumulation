% function Fig_Whole_Cell_Donnan_Anion_vs_Cation_v1


function Fig_Whole_Cell_Donnan_Anion_vs_Cation_v1

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
        %params.membrane_x = 1e-9;   			% Position of membrane patch back from the tip (in m)
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
    params.V_hold = -100e-3;
    
    params_im0 = params;
    for j = N_patch:params.membrane_N
            r = (j-N_patch)/(params.membrane_N -N_patch) * rho_im;
            params.rho(:,j) = [params.rho(1,j), params.rho(2,j), params.rho(3,j), params.rho(4,j)-r,  r]' ;
    end
    params_im100 = params;
    
    % Initialiaze
    params_im0 = evolve_backward_euler(params_im0,inf);
    params_im100 = evolve_backward_euler(params_im100,inf);
    
    % I versus t for chloride an K+ conductance
    figure(1); clf;
    fnameA = 'Fig_Donnan_Anion_vs_Cation_partA';
    sizeA = [ 8 11];
    save_figure(fnameA,sizeA); 
    
    % Chloride conductance
    idx = 3; % Cl-
    I_mem = [0, 0, -1e-10, 0, 0 ]'; % Initial current of -100pA by K+
    V_mem = -100e-3; % Voltage in mV
    P_GHK = calculate_PGHK_from_I(I_mem, params_im0.rho_pip, params_im0.rho_bath, V_mem-I_mem(idx) * R_pip, params_im0.z, params_im0.A_cell, params_im0.T);
    perm  = [ 0, 0, P_GHK(3), 0, 0 ]'; 		% Membrane permability to each ion in m/s
    
    V_rev = params_im0.R * params_im0.T / params.F * log( params_im0.rho_pip(idx) / params_im0.rho_bath(idx));
    fprintf('\n\n Chloride Permeability = %e m/s. V_Nernst_Cl- = %.2f mV. g_Cl-(chord) = %.2f nS.', P_GHK(idx), V_rev * 1e3, I_mem(idx)/(V_mem-V_rev)*1e9);
    
    intervals.duration =                [2,     0.1,    19.9,    0.1,    1.9 ] ; % Duration of each of intervals
    intervals.N_eval =                  [2,     5,      20,     5,      2 ] ;     % Number of timepoints
    intervals.variables.membrane_perm = {0*perm,perm,   perm,   0*perm, 0*perm};  % Membrane permeability
    protocol = create_protocol( intervals );
    [paramst, results1 ] = evolve_backward_euler(params_im0, protocol);
    [paramst, results100 ] = evolve_backward_euler(params_im100, protocol);
    
    subplot(4,2,1);   
        plot(protocol.t, results1.I_pip*1e9,'k'); hold on;   axis([0 24 -0.1 0]);
        axis off; 
        plot([10 20],-0.02 + [0 0 ], 'k', 'Linewidth',1); % 10s time bar
        plot([1 1], [-0.075 -0.025], 'k', 'Linewidth',1); % 50pA scale bar
    subplot(4,2,3);   plot(protocol.t, results1.rho_intra(3,:), 'k'); axis([0 24 5 10]);
        box off; ylabel('\rho_{Cl-}'); set(gca,'YTick',[6 8 10],'TickDir','out');
        h = gca; h.XAxis.Visible = 'off';   
    
    subplot(4,2,2);   
        plot(protocol.t, results100.I_pip*1e9,'k'); hold on; axis([0 24 -0.1 0]);
        axis off; 
    subplot(4,2,4);   plot(protocol.t, results100.rho_intra(3,:),'k');axis([0 24 5 10]);
        box off; h = gca; h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off';
    
    % Now calculate I_init/ss as function of anion concentration.
    rho_im = [0:10:100]; 
    V_hold = -100e-3;
    for k = 1:length(rho_im)
        
        % Set immobile charge density
        paramst = params;
        for j = N_patch:params.membrane_N
                r = (j-N_patch)/(paramst.membrane_N -N_patch) * rho_im(k);
                  paramst.rho(:,j) = [paramst.rho(1,j), paramst.rho(2,j), paramst.rho(3,j), paramst.rho(4,j)-r,  r]' ;
        end
        
        % Initialiaze state
        paramst.V_hold = V_hold; paramst.perm = 0 * perm;  paramst = evolve_backward_euler(paramst, inf);
        
        % Measure initial and ss current
        paramst.membrane_perm = perm;
        [paramst, result1] = evolve_backward_euler(paramst, 0.01);
        [paramst, result1] = evolve_backward_euler(paramst, 0.01);
        I_init(k) = result1.I_pip;
        [paramst, result1] = evolve_backward_euler(paramst, inf);
        I_ss(k) = result1.I_pip;
    end
    
    figure(2); clf;
    fnameB = 'Fig_Donnan_Anion_vs_Cation_partB';
    sizeB = [ 4 11];
    save_figure( fnameB, sizeB); 
    subplot(2,1,1);
    plot(rho_im, -1 * I_init*1e12,'k'); hold on
    plot(rho_im, -1 * I_ss*1e12); hold on
    axis([0 100 50 100]);
    ylabel('|I_{Cl-}| (pA)');
    xlabel('\rho_{immobile} (mM)'); set(gca,'TickDir','out');
    
        
    % K conductance
    V_hold = 50e-3;
    idx = 1; % Cl-
    I_mem = [2e-9, 0, 0, 0, 0 ]'; % Initial current of 2nA by K+
    P_GHK = calculate_PGHK_from_I(I_mem, params_im0.rho_pip, params_im0.rho_bath, V_hold-I_mem(idx) * R_pip, params_im0.z, params_im0.A_cell, params_im0.T);
    V_rev = -1 * params_im0.R * params_im0.T / params.F * log( params_im0.rho_pip(idx) / params_im0.rho_bath(idx));
    fprintf('\n\n K_ Permeability = %e m/s. V_Nernst_K+ = %.2f mV. g_K+(chord) = %.2f nS.', P_GHK(idx), V_rev * 1e3, I_mem(idx)/(V_hold-V_rev)*1e9);
    perm = [P_GHK(idx), 0, 0, 0, 0]';
    
    intervals.variables.membrane_perm = {0*perm,perm,   perm,   0*perm, 0*perm};  % Membrane permeability
    protocol = create_protocol( intervals );
    params_im0.V_hold = V_hold;   params_im0 = evolve_backward_euler(params_im0,inf);
    params_im100.V_hold = V_hold; params_im100 = evolve_backward_euler(params_im100,inf);
     
    [paramst, results1 ] = evolve_backward_euler(params_im0, protocol);
    [paramst, results100 ] = evolve_backward_euler(params_im100, protocol);
    
    figure(1);
    subplot(4,2,5);   
        plot(protocol.t, results1.I_pip*1e9,'k'); hold on;   axis([0 24 0 2.4]);
        axis off; 
        plot([10 20],0.5 + [0 0 ], 'k', 'Linewidth',1); % 10s time bar
        plot([1 1], [0 1], 'k', 'Linewidth',1); % 1nA scale bar
    subplot(4,2,7);   
        plot(protocol.t, results1.rho_intra(1,:),'k'); axis([0 24 100 200]);
        box off; ylabel('\rho_{K+}'); set(gca,'YTick',[100 140 180],'TickDir','out');
        h = gca; h.XAxis.Visible = 'off';   
    
    subplot(4,2,6);   
        plot(protocol.t, results100.I_pip*1e9,'k'); hold on; axis([0 24 0 2.4]);
            axis off
    subplot(4,2,8);   plot(protocol.t, results100.rho_intra(1,:),'k');axis([0 24 100 200]);
            axis off
    
    save_figure(fnameA,sizeA); 
            
    % Plot effect of immobile anions versus concentraiton
     % Now calculate I_init/ss as function of anion concentration.
    rho_im = [0:10:100]; 
    V_hold = 50e-3;
    for k = 1:length(rho_im)
        
        % Set immobile charge density
        paramst = params;
        for j = N_patch:params.membrane_N
                r = (j-N_patch)/(paramst.membrane_N -N_patch) * rho_im(k);
                  paramst.rho(:,j) = [paramst.rho(1,j), paramst.rho(2,j), paramst.rho(3,j), paramst.rho(4,j)-r,  r]' ;
        end
        
        % Initialiaze state
        paramst.V_hold = V_hold; paramst.perm = 0 * perm;  paramst = evolve_backward_euler(paramst, inf);
        
        % Measure initial and ss current
        paramst.membrane_perm = perm;
        [paramst, result1] = evolve_backward_euler(paramst, 0.01);
        [paramst, result1] = evolve_backward_euler(paramst, 0.01);
        I_init(k) = result1.I_pip;
        [paramst, result1] = evolve_backward_euler(paramst, inf);
        I_ss(k) = result1.I_pip;
    end
    
    figure(2); 
    subplot(2,1,2);
    plot(rho_im, I_init*1e9,'k'); hold on
    plot(rho_im, I_ss*1e9); hold on
    axis([0 100 1.5 2.4]);
    ylabel('|I_{K+}| (nA)');
    xlabel('\rho_{immobile} (mM)');set(gca,'TickDir','out');
    save_figure(fnameB,sizeB); 
    
    