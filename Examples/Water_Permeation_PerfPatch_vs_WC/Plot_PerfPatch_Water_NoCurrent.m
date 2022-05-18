% function Fig_WaterFlow_NoCurrent
%
%   Plot water flow, volume, concentration, and molar currents when
%   Cell containing 150mM KCl is exposed to solution containing 150mM KCl +
%   60mM sucrose
%
%   Consider 2 cases -> perforated patch
%                    -> whole-cell


function Fig_WaterFlow_NoCurrent

	% Configuration parameters
        % For now we will jury rig the perforated patch config
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
                       
	% Solution Parameters
        params.D =  [1, 1.0382, 0.25]' * 1.957e-9; 	% Diffusion coefficients of K+ , Cl-, sucrose ions
                                                % Sucrose ~ 5e-10 m^/2        English, A.C., and M. Dole. 1950. Diffusion of Sucrose in Supersaturated Solutions. J. Am. Chem. Soc. 72:3261–3267. doi:10.1021/ja01163a132.
		params.z =  [1, -1, 0 ]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [150, 150, 0  ]';  	% 150mM KCl in bath
		params.rho_pip =  [150, 150, 0  ]';  	% 150mM KCl in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
        params.theta = 10 * pi/180;  % Tip Angle
        [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
        params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
    	params.Nseg  =  21; 		% Number of segments between the pipette electrode and tip
		params.V_hold = 0e-3; 	% Initial holding voltage
        params.Rseries = 0e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Cell and Membrane parameters
        params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0 0  ]'; 		% Membrane permability to each ion in m/s
        params.membrane_Pf = 5e-5   ; 			% Membrane permeability to water (m/s)
        params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
        params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell

    % Now set up 2 configs 
        % Perforated Patch Config
        params_perf = params;
        params_perf.Constraints.type = 'solvent_cell_volume';  % Allow solvent flow into cell to change cell volume.
        params_perf = initialize_system(params_perf);
    
        % Report Access Resistance and Series Resistance Correction
       	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params_perf.r_tip * 2e6);
       	fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params_perf)/1e6), sum(params_perf.Rseries/1e6))

        % Set up perforated patch
        d_patch = 3e-6; % Concentration ramps linearly starting from d_pip back in pipette tip
        [tmp, N_patch ] = min(abs(params_perf.x + d_patch)); 
        params_perf.patch_func = @patch_current_GHK;
        params_perf.patch_perm = [ 2, 2 , 0 ]' * 1e-3; 		% Patch permability to each ion in m/s
        params_perf.patch_N = N_patch;
    
        % Initialize all remaining parameters
        % Now need to initialize before breaking patch.
        params_perf = evolve_backward_euler(params_perf,0.1);
        params_perf = evolve_backward_euler(params_perf,10);
        params_perf.t = -2;
    
        % Conventional Whole-Cell Config
        params_wc = params;
        params_wc.Constraints.type = 'solvent';  % Allow solvent flow but keep all compartment volumes constant.
        params_wc = initialize_system(params_wc);
    
        % Report Access Resistance and Series Resistance Correction
       	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params_wc.r_tip * 2e6);
       	fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params_wc)/1e6), sum(params_wc.Rseries/1e6))
    
        % And initialize state
        params_wc = initialize_system(params_wc);
        params_wc = evolve_backward_euler(params_wc,inf);
        params_wc.t = -2;

        
        % Protocol for exposing cells to 60mM sucrose
        t_val = [-2 ,[0:0.01:0.1, 0.2:0.1:1.9, 2:0.5:10]];
        rho_bath_osm = [150, 150, 60]'; % Add 60mM sucrose to external solution
        params_set = {params_perf, params_wc}; % Two configs to evolve
        fnames = {'fig_Perf_water','fig_WC_water'};
        
        % Plot positions
        left = [0.1, 0.6];
        bottom = [0.925, 0.75, 0.55, 0.35, 0.15, 0.2];
        width = [0.35, 0.35];
        height = [ 0.02, 0.15, 0.15, 0.15, 0.15, 0.3];
      
        for j=1:length(params_set) % Do each config in turn
        
            params = params_set{j};
        
            for k=1:length(t_val)
                [params, results] = evolve_backward_euler(params, t_val(k));
        
                % Step increase in osmotic pressure
                if (t_val(k) >= 0)
                    params.rho(:, params.Nseg) = rho_bath_osm;
                end
        
                rho(:,k) = results.rho_intra;
                V(k) = (results.V_mem - params.V_hold) ;
                Lambda(k) = params.Lambda(params.membrane_N);
                fsol(k) = params.Jsol(params.membrane_N);
                fpip(k) = params.Jsol(N_patch);
                J_pip_K(k)  = params.J(1,N_patch);
                J_pip_Cl(k) = params.J(2,N_patch);
            end
    
            figure(1); clf;
            t_range = [-2 10];
            
            % Bar for application of solution
            subplot('position',[left(1), bottom(1), width(1), height(1)]);
            plot( [-2 0], [0 0],'k','Linewidth',1); hold on
            plot( [ 0 10], [0 0],'b','Linewidth',1); hold on
            axis([t_range, -1, 1]); axis off;
            
            % Flow
            subplot('position',[left(1), bottom(2), width(1), height(2)]);
            plot(t_val, fsol*1e18,'k'); hold on
            plot(t_val, fpip*1e18,'m:'); hold on
            axis([t_range, 0, 60]);
            box off; h = gca;  h.XAxis.Visible = 'off';    
            ylabel('f (fL/s)');
            
            % Cell Volume
            subplot('position',[left(1), bottom(3), width(1), height(3)]);
            plot(t_val, Lambda*1e15,'k'); hold on
            axis([t_range, 0.70, 1.02]);
            box off; h = gca;  h.XAxis.Visible = 'off';    
            ylabel('\Lambda_{cell} (pL)');

            % Concentration
            subplot('position',[left(1), bottom(4), width(1), height(4)]);
            plot(t_val, rho(1,:),'b'); hold on;
            plot(t_val, rho(2,:),'g:');
            axis([t_range, 148, 172]);
            box off; h = gca;  h.XAxis.Visible = 'off';    
            ylabel('\rho_{cell} (mM)');
        
            % Solvent flux
            subplot('position',[left(1), bottom(5), width(1), height(5)]);

            plot(t_val,  J_pip_K *1e15,'b'); hold on
            plot(t_val, J_pip_Cl *1e15,'g:'); hold on
            plot(4+ [0 2], [0 0],'k','Linewidth',1); % Add time bar
            switch j
                case 1
                    axis([t_range, -6.4, 2.2]);
                case 2
                    axis([t_range, -0.2, 8.4]);
            end
            box off; h = gca;  h.XAxis.Visible = 'off';    
            ylabel('J_{pip} (fmol/s)');

            % Plot concentration versus position 
            subplot('position',[left(2), bottom(6), width(2), height(6)]);
            
            
            % Rescale distance from pipette tip
            xval = 1e6 * params.x; % Xposition in microns

            
            % Rescale pipette distances
            idx_pip = [1:(params.membrane_N - 1)]; % Pipette indices
            xpip = xval(idx_pip); % Pipette x values
            lpip = 18; % Pick scaling length in microns
            xpiprs = lpip * ( lpip./(lpip - xpip) - 1);
            xval(idx_pip) = xpiprs;

            % Rescale cell distances
            xcell = 8; % Cell wall at 8
            xval(params.membrane_N) = xcell;
            idx_cell = params.membrane_N + [-1, 0]; 
            
            
            plot(xval(idx_pip), params.rho(1,idx_pip),'b'); hold on
            plot(xval(idx_pip), params.rho(2,idx_pip),'g:'); hold on
            plot(xval(idx_cell), params.rho(1,idx_cell),'b'); hold on
            plot(xval(idx_cell), params.rho(2,idx_cell),'g:'); hold on
            plot( [xcell, 20 ], params.rho(1,params.Nseg) + [0, 0],'b'); hold on
            plot( [xcell, 20], params.rho(2,params.Nseg)+ [0, 0],'g:'); hold on
            axis([-18 20 145 175]) ; 
            ylabel('\rho (mM)')
            box off
            
            
            save_figure(fnames{j},[ 10 10]);       
            
            
       end
            
            
            
       
        
        
        