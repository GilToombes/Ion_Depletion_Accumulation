% function Fig
%
%	Example of a calculation for the whole-cell configuration.
%       K-channel time dependent pulse.
% 

function Fig5_WholeCell_IVrhomu_vs_time_v2
    
 	% Configuration parameters
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	% Solution Parameters
		params.D =  [1, 1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ and Cl-.
		params.z =  [1, -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [150, 150  ]';  	% 150mM KCl in bath
		params.rho_pip =  [150, 150  ]';  	% 150mM NaCl in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
        params.theta = 10 * pi/180;  % Tip Angle
        [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
		params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
        params.Nseg  =  3; 		    % Number of segments between the pipette electrode and bath electrode
		params.V_hold = 100e-3; 	% Initial holding voltage
        params.Rseries = 0e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
        
	% Cell and Membrane parameters
        %params.membrane_x = 1e-6;   			% Position of membrane 
		params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0 ]'; 		% Membrane permability to each ion in m/s
        % params.membrane_Pf = 0   ; 			% Membrane permeability to water (m/s)
		params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell
        
	% Initialize all remaining parameters
    params = initialize_system(params);
    
    % Calculate membrane permeability for an initial current of 2nA K+
    I_mem = [2e-9, 0 ]'; % Initial current of 2nA by K+
    P_GHK = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, params.V_hold, params.z, params.A_cell, params.T);
    
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', resistance(params)/1e6, params.Rseries/1e6)
    fprintf('\n Membrane Permeability = [%e, %e] m/s', P_GHK);
    
	% Equilibrate initial state
	[params, results] = evolve_backward_euler(params, inf);
    params.t = -20e-6; % Start time at -20 microseconds to show initial step.
    
    % Need to catch each time step to plot relevant quantities
    t_val_q = 1e-6 * [-20:230]; % Timescale for electric steady-state
    t_val_rho = 2e-2 * [1:1000]; % Timescale for chemical steady-state
    t_val = [t_val_q, t_val_rho];
    Nq = length(t_val_q);
    idx_q = [1:Nq];
    Ne = length(t_val);
    idx_rho = [(Nq+1):Ne];
    

    % Store results
    idx_p = 1; % index for pipette current
    idx_m = params.Nseg-1; % index for membrane current
    idx_cell =params.Nseg-1; % Index for cell compartment
    JK = zeros(Ne, 1);
    JCl =zeros(Ne, 1);
    rho_K = zeros(Ne, 1);
    rho_Cl = zeros(Ne, 1);
    V = zeros(Ne, 1);
    Imem = zeros(Ne,1);
    Ielec = zeros(Ne,1);
    
    % Evolve State
    for j=1:Ne
        
        if (t_val(j)>=0)
           params.membrane_perm = P_GHK;  % Make the membrane K-selective after t = 0
        end
    
        [params, results] = evolve_backward_euler(params, t_val(j));
        JK(j) = params.J(1,idx_p);
        JCl(j) = params.J(2,idx_p);
        rho_K(j) = params.rho(1,idx_cell);
        rho_Cl(j)= params.rho(2,idx_cell);
        V(j) = params.V(idx_cell);
        Imem(j) = (params.J(1,idx_m) - params.J(2,idx_m)) *params.F*1e9;
        Ielec(j) = (params.J(1,idx_p) - params.J(2,idx_p)) *params.F*1e9;
    end
    
    
    figure(1); clf
    fig_dim = [12 12];
    fname = 'Fig_Example_WholeCell_TimeDependence_v2';
    save_figure(fname,fig_dim)
       
        yp = [0.08, 0.32, 0.57, 0.75]; % Position of rows
        yp = [0.0667    0.2667    0.4750    0.6250    0.8333]; % Position of rows
        yh = [0.1667    0.1667    0.1250    0.1667    0.1000]; % Height of rows
        xp = [0.1 0.5]; % Position of columns
        xw = [0.375 0.375];
        t_range_q = [-20 230];
        
        % Conductance on top
        subplot('position', [xp(1) yp(5) xw(1), yh(5)]);
        y_range = [0 21]; % nS
        tbase = t_val_q*1e6;
        g = Imem ./ (V + (params.R * params.T / params.F) * log( rho_K/ params.rho_bath(1))) ;
        tbase = t_val_q*1e6;
        plot(tbase, g(idx_q) ,'b'); hold on
        plot(tbase, g(idx_q) * 0 ,'Color',[0.5 0.5 0.5]); hold on
        axis([t_range_q, y_range]);
        box off; h = gca;  h.XAxis.Visible = 'off';    
        ylabel('g_{mem} (nS)'); 
        
        subplot('position', [xp(2) yp(5) xw(2), yh(5)]);
        tbase = t_val_rho;
        plot(tbase, g(idx_rho) ,'b'); hold on;
        plot(tbase, g(idx_rho)*0 ,'Color',[0.5 0.5 0.5]); hold on;
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        axis([0 20 y_range]);
        
        % Currents on 4
        subplot('position', [xp(1) yp(4) xw(1), yh(4)]);
        y_range = [0 1];y_range = [0 2];
        tbase = t_val_q*1e6;
        plot(tbase, Ielec(idx_q) ,'k--'); hold on; plot(tbase,Imem(idx_q),'k'); plot(tbase, 1e9 * JK(idx_q) * params.F,'b');plot(tbase, -1e9 * JCl(idx_q) * params.F,'g');
        axis([t_range_q, y_range]);
        box off; h = gca;  h.XAxis.Visible = 'off';    
        ylabel('I (nA)'); 
        
        subplot('position', [xp(2) yp(4) xw(2), yh(4)]);
        tbase = t_val_rho;
        plot(tbase, Imem(idx_rho) ,'k'); hold on; plot(tbase,Ielec(idx_rho),'k--'); plot(tbase, 1e9 * JK(idx_rho) * params.F,'b');plot(tbase, -1e9 * JCl(idx_rho) * params.F,'g');
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        axis([0 20 y_range]);
             
        
        % Voltage Next
        idx_v = 3;
        subplot('position', [xp(1) yp(idx_v) xw(1), yh(idx_v)]);
        y_range = [95 100]; y_range = [90 100];
        
        tbase = t_val_q*1e6;
        plot(tbase,V(idx_q)*1e3,'k'); 
        axis([t_range_q, y_range ]);
        box off; h = gca;  h.XAxis.Visible = 'off';    
        ylabel('V_{cell} (mV)');
        
        subplot('position', [xp(2) yp(idx_v) xw(2), yh(idx_v)]);
        tbase = t_val_rho;
        plot(tbase,V(idx_rho)*1e3,'k');
        axis([0 20 y_range ]);
        box off ; h = gca; h.YAxis.Visible = 'off'; h.XAxis.Visible = 'off';           
        
               
        % Concentrations Third
        idx_v = 2;
        subplot('position', [xp(1) yp(idx_v) xw(1), yh(idx_v)]);
        y_range = [125 150.1]; y_range = [110 150.1];
        
        tbase = t_val_q*1e6;
        plot(tbase,rho_K(idx_q),'b'); hold on;
        plot(tbase,rho_Cl(idx_q),'g--'); 
        axis([t_range_q, y_range ]);
        box off; h = gca;  h.XAxis.Visible = 'off';  
        ylabel('\rho_{cell} (mV)');set(gca,'YTick',[110 130 150]);
        
        subplot('position', [xp(2) yp(idx_v) xw(2), yh(idx_v)]);
        tbase = t_val_rho;
        plot(tbase,rho_K(idx_rho),'b'); hold on;
        plot(tbase,rho_Cl(idx_rho),'g--'); 
        axis([0 20 y_range ]);
        box off; h = gca; h.YAxis.Visible = 'off';      h.XAxis.Visible = 'off';   
        
        % Chemical Potential Fourth
        idx_v = 1;
        subplot('position', [xp(1) yp(idx_v) xw(1), yh(idx_v)]);
        y_range = [-0.2 0.4]; y_range = [-0.8 0.4];
        delta_mu_K  = log(rho_K/params.rho_pip(1))  -  params.F * (params.V_hold - V) / (params.R * params.T);
        delta_mu_Cl = log(rho_Cl/params.rho_pip(2)) +  params.F * (params.V_hold - V) / (params.R * params.T);
        
        tbase = t_val_q*1e6;
        plot(tbase,delta_mu_K(idx_q),'b'); hold on;
        plot(tbase,delta_mu_Cl(idx_q),'g'); 
        axis([t_range_q, y_range ]);
        box off
        ylabel('\Delta\mu_{cell}/RT '); set(gca,'YTick',[-0.4, 0, 0.4, 0.8]);
        xlabel('(\mus)'); set(gca,'XTick',[0 50 100 150 200]); 
        
        
        subplot('position', [xp(2) yp(idx_v) xw(2), yh(idx_v)]);
        tbase = t_val_rho;
        plot(tbase,delta_mu_K(idx_rho), 'b'); hold on;
        plot(tbase,delta_mu_Cl(idx_rho),'g'); 
        axis([0 20 y_range ]);
        box off ; h = gca; h.YAxis.Visible = 'off';        
        
        xlabel('(s)'); set(gca,'XTick',[5 10 15 20]); 
    
   save_figure(fname,fig_dim)
        
        