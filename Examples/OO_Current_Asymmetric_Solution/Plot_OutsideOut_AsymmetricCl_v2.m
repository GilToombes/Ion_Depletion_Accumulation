% function Fig12_OutsideOut_AsymmetricCl_v2
%
%   Show effect of asymmetric Cl concentration on anion current
%   Outside-out pipette containing 150mM K+, 10mM Cl-, and 140mM gluconate

%   Plot concentration and voltage-change at the end of pipette containing
%   150mM KCl as a function of IK+*Raccess, and |INa+| * Raccess


function Fig_OutsideOut_Cl_Asymmetric

	% Configuration parameters
		params.config = 'outside-out';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	% Solution Parameters
        params.D =  [1, 0.33, 1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ , gluconate and Cl-.
		params.z =  [1, -1 -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [150, 10, 140  ]';  	% 150mM K+, 10mM  Gluconate-, 140mM Cl- in bath
		params.rho_pip =  [150, 140, 10  ]';  	% 150mM K+, 140mM Gluconate-, 10mM  Cl- in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
        params.theta = 10 * pi/180;  % Tip Angle
        [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
        params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
        params.Nseg  =  11; 		% Number of segments between the pipette electrode and tip
		params.V_hold = 0e-3; 	% Initial holding voltage
        params.Rseries = 4.5e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Cell and Membrane parameters
        params.membrane_x = 1e-9;   			% Position of membrane patch back from the tip (in m)
		params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0 0e-4]'; 		% Membrane permability to each ion in m/s
        %params.membrane_Pf = 5e-5   ; 			% Membrane permeability to water (m/s)
		%params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		%params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell

	% Initialize all remaining parameters
    params = initialize_system(params);
    
    % Report Access Resistance and Series Resistance Correction
    fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
   	fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', resistance(params)/1e6, params.Rseries/1e6)

    perm =   1.0e-4 * [0 0 1]' ;
    intervals.duration =                [50e-3, 1e-5,   9.9e-4, 9e-3,   9e-2,   5e-3,   45e-3 ] ; % Duration of each of intervals
    intervals.N_eval =                  [2,     40,     10,     20,     20,     10,     20] ;     % Number of timepoints
    intervals.variables.membrane_perm = {0*perm,perm,   perm,   perm,   perm,   0*perm, 0*perm};  % Membrane permeability
    protocol = create_protocol( intervals );
    V = [25e-3, -100e-3];
    pf = [0.25:0.25:1];
    
    figure(1); clf;
    
    params.Ct = 0 * params.Cl;
    params.Cl = 0 * params.Cl;

    paramst = params;
    Vrev = (params.R * params.T/params.F) * log( params.rho_pip(3) / params.rho_bath(3)  );
    
    for k = 1:length(V)
    for j=1:length(pf)
        
        % Equilibrate voltage
        paramst         = params;
        paramst.V_hold  = V(k);
        [paramst, results] = evolve_backward_euler(paramst, inf);
    
        p = perm * pf(j);
        intervals.variables.membrane_perm = {0*p, p,   p, p, p, 0*p, 0*p};  % Membrane permeability
        protocol = create_protocol( intervals );
        
        [ptemp, results] = evolve_backward_euler(paramst, protocol);
        
        paramst.membrane_perm = p; 
        [paramst, results_f] = evolve_backward_euler(paramst, inf);
        
        subplot(3,2,k);
        plot(protocol.t, results.I_pip*1e9,'k'); hold on;
        
        subplot(3,2,k+2);
        plot(protocol.t, results.rho_intra(3,:),'k'); hold on;
        
        I_max(k,j) = max(abs(results.I_pip)) * sign(results_f.I_pip);
        I_fin(k,j) = results_f.I_pip;
        rho_fin(j) = results_f.rho_intra(3);
        
    end
    end
    
        subplot(3,2,1)
            title('V_{hold} = 25mV')
            axis([0 0.2 -0.01 1.5]);
            hold on
            plot([0.05 0.15], [0.2, 0.2],'k','LineWidth',1); % 100ms time bar
            plot([0.02 0.02], [0.25,1.25],'k', 'LineWidth',1); % 1 nA scale bar
            axis off
            
        subplot(3,2,2)
            title('V_{hold} = -100mV')
            axis([0 0.2 -0.21 0.005 ]);
            plot( [0.05, 0.15], [-0.175 -0.175], 'k', 'Linewidth',1); % 100ms time bar
            plot([0.02 0.02], [-0.05,-0.15],'k', 'LineWidth',1); % 100 pA scale bar    
            axis off;
            
            
        subplot(3,2,3)
            axis([ 0 0.2 5 55]);
            ylabel('\rho_{Cl} (mM)')
            box off;
            set(gca,'XTick', []);
            set(gca,'YTick',[10 30 50],'TickDir','out');
            
            
        subplot(3,2,4)
            axis([ 0 0.2 5.5 10.5]);
            ylabel('\rho_{Cl} (mM)')
            
            set(gca,'XTick', []); box off;
            set(gca,'YTick',[6,8, 10],'TickDir','out');
            
        subplot(3,2,5)
            Ii = [0, I_max(1,:)*1e9];
            If = [0, I_fin(1,:)*1e9];
            plot( Ii, If,'k'); hold on
            plot( Ii, If,'ko'); hold on
            plot([0 1.5], [0 1.5],'k:');
            axis image; box off
            axis([0 1.6 0 1.6]);
            set(gca,'TickDir','out','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5]);
            xlabel('|I_{init}| (nA)');
            ylabel('|I_{final}| (nA)');
            
        subplot(3,2,6)
            Ii = [0, -I_max(2,:)*1e12];
            If = [0, -I_fin(2,:)*1e12];
            plot( Ii, If,'k'); hold on
            plot( Ii, If,'ko'); hold on
            axis image; box off;
            axis([0 215 0 215]); 
            plot([0 200], [0 200],'k:');
            set(gca,'TickDir','out','XTick',[0 100 200],'YTick',[0 100 200]);
            xlabel('|I_{init}| (pA)');
            ylabel('|I_{final}| (pA)');
            
        
	save_figure('Fig_OO_Cl_Asymmetric',[ 10 12.5]);       
            
    function result = I_GHK(params, V, rho_in, rho_out)
        
        idx = 1;
        Am = params.A(params.membrane_N);
        u = params.z(idx)*params.F*V/(params.R*params.T) ; 
        uz = (abs(u)<1e-6); % Uncharged particles and zero membrane voltge
        result = params.F * Am * params.membrane_perm(idx) .* (u + uz) ./ (1 - exp(-u) + uz) .* (rho_in - rho_out.*exp(-u));