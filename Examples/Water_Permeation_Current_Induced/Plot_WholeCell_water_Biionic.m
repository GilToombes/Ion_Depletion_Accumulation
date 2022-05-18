% function Fig13_Ionic_Selectivity_v1
%
%   Illusrate effect of ion depletion/accumulation on measurement of ionic
%   selectivity during prolonged activation.
%
%   Simulation Details 
%               (based on fig1/2 of  Nat. Neurosci. volume 18, pages 1577–1583 
%
%              Whole Cell
%
%              Extracellular Solution = 150mM NMDG-Cl or 150mM Na-Cl
%
%              Intracellular Solution = 150mM NaCl
%              Pipette Resistance = 5.00 MOhms
%              Pipette Solution Conductance = 1.93 S/m. 
%              Pipette Tip Diameter = 0.75 microns.
%   
%              Membrane Permeability = [ 1 1/20  0]' * 3e-6 m/s
%                       in presence of agonist
%                       i.e. P_Na+/P_NMDG+ = 20 / 1
%              Cell Area =  1e-9 m^2  = 1000 square microns
%              Cell Volume = 1e-15 m^3 = 1pL
%
%              Protocol 1 - Ext Solution is constantly 150mM NMDG-Cl
%                   
%                   V_hold = -90mV
%                   Voltage ramps - -90mV to -20mV over 0.5 seconds
%                   V_stim = -60mV  voltage during prolonged activation
%                   T_stim = 15s ; activation duration
%
%                   Initialize for 0.5 seconds at V_hold
%                   Apply agonist and perform Ramp
%                   Wait 0.5 seconds at V_hold in absence of agonist
%                   Apply agonist for T_stim seconds holding at V_stim
%                   Wait 0.5 seconds at V_hold in absence of agonist
%                   Apply agonist and perform Ramp
%                   Wait 0.5 seconds at V_hold in absence of agonist
%               
%                    Plot :  
%                          Holding Voltage vs t
%                          Current vs t
%                          Individual Ionic currents vs t
%                          Cell concentrations vs t
%                          Delta_V = Vpip-Vcell vs t
%                          I versus V for voltage ramps
%                   For 2 ramps
%                           Slope Conductance = 106.10 nS. Reversal voltage = -74.73 mv
%                           Slope Conductance = 87.64 nS. Reversal voltage = -44.57 mv
%
%               Then repeat protocol 1 with different evels of membrane
%               permeability and plot
%                           Vrev of final voltage ramp 
%                               versus
%                           Slope conductance of initial voltage ramp
%
%               Compare Pf = 0 m/s and Pf = 5e-5 m/s for water permeation    

function Fig_Ionic_Selectivity_Water_v1

	% Configuration parameters
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        %params.Constraints.type = 'basic_perfusion';  % Option specifying which parameters are variable
        params.Constraints.type = 'solvent';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
                                            
	% Solution Parameters
		params.D =  [0.682, 0.33, 1.0382]' * 1.957e-9; 	% Diffusion coefficients of Na+, NMDG+, Cl-
		params.z =  [1, 1, -1]';                    % Relative charge of        Na+, NMDG+, Cl-
		params.rho_bath = [0, 150, 150  ]';           % 150mM NMDG-Cl in bath
		params.rho_pip =  [150, 0, 150]';           % 150mM NaCl in pipette
        
	% Pipette parameters
        R_pip = 5e6;                % Pipette resistance in Ohms
        params.theta = 10 * pi/180;  % Tip Angle
        params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
        params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
		[tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
        params.Nseg  =  11; 		% Number of segments between the pipette electrode and tip
		params.V_hold = -90e-3; 	% Initial holding voltage Vhold = -60mV
%       params.Rseries = 4e6;       % Series resistance correction in Ohms
%       params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Cell and Membrane parameters
        %params.membrane_x = 1e-9;   			% Position of membrane patch back from the tip (in m)
		params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
        mem_p = [ 1 1/20  0]' * 3e-6;           % Membrane permability to each ion in m/s
		params.membrane_perm = 0 * mem_p;                   % Membrane permability to each ion in m/s
        membrane_Pf = 5e-5;	% Membrane permeability to water (m/s)
        params.membrane_Pf = membrane_Pf   ; 		
		params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell

	% Initialize all remaining parameters
    params = initialize_system(params);
    %params.Constraints.segments{2} = [2,params.Nseg-1]; % Switch to current-clamp mode
    %params.Constraints.segments{1}  = [1, params.Nseg]; % Allow perfusion of solution in bath and pipette
    %params = initialize_system(params);
    %params.Constraints.segments{2}

    params.Cl = params.Cl * 0; % Compensate capacitance.
    params.Ct = params.Ct * 0; % Compensate capacitance.
        
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. ', sum(resistance(params)/1e6));
    fprintf('\n Membrane Permeability = [%e, %e, %e] m/s', mem_p');
    
    
    % Now define protocols
    V_hold1 = -90e-3; % Holding voltage
    V_hold2 = -30e-3; % Holding voltage
    
    T_rest = 0.5; % Rest periods of 1 second
    N_rest = 5; 
    T_stim = 15; % 15 seconds stimulation
    N_stim = 30; % Number of timepoints during stimulation
    V_stim = -60e-3; % Stimulation voltage
    rho_bstim = [0, 150, 150  ]';           % 150mM NMDG-Cl in bath
	rho_b = [150 0 150]'; % Use 150mM NMDG-CL during ramp
   
    % Ramps
    T_ramp = 0.5; % 500ms voltage ramps
    N_ramp = 60; %  Number of steps during ramp
    V_ramp1 = [-90e-3, -20e-3]; % -90mV to -20mV voltage ramp
    V_ramp2 = [-30e-3, 30e-3]; % -90mV to -20mV voltage ramp
    
    intervals1.duration=                 [T_rest,    T_ramp,     T_rest,     T_stim,     T_rest,     T_ramp,     T_rest]; % Duration of each of intervals
    intervals1.N_eval =                  [2,         N_ramp,     N_rest,     N_stim,     N_rest,     N_ramp,     N_rest]; % Number of timepoints during each interval
    intervals1.variables.V_hold =        {V_hold1,   V_ramp1,    V_hold1,    V_stim,     V_hold1,    V_ramp1,    V_hold1}; % Holding voltage
    intervals1.variables.membrane_perm=  {0*mem_p,   mem_p,      0*mem_p,    mem_p,      0*mem_p,    mem_p,      0*mem_p};% Membrane Permeability.
    
    protocol1 = create_protocol(intervals1); % Create Protocol 
    
    % Add Vrev versus Slope conductance
    sf = [0.001, 0.01, 0.02, 0.03, 0.05, 0.065, 0.08, 0.1,0.12, 0.145, 0.17, 0.195, 0.225, 0.26 0.3, 0.35, 0.42, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.55];
    % Calculate for different membrane permeabilities
    
    for j=1:length(sf)
    
       % Initialize state
       
       % Change permeability in protocol
       mp = mem_p * sf(j);
       intervals1.variables.membrane_perm=  {0*mp,   mp,      0*mp,    mp,      0*mp,    mp,      0*mp};% Membrane Permeability.
       protocol1 = create_protocol(intervals1); % Create Protocol 
    
       % Simulate Protocol 1    
       
       % No water first
       paramst = params; paramst.membrane_Pf = 0;
       paramst.V_hold = V_hold1 ;  [paramst, results] = evolve_backward_euler(paramst, inf); 
       [paramst, results_nw, paramst_full_nw] = evolve_backward_euler(paramst, protocol1);

       % Now with water permeating
       paramst = params; paramst.membrane_Pf = membrane_Pf;
       paramst.V_hold = V_hold1 ;  [paramst, results] = evolve_backward_euler(paramst, inf); 
       [paramst, results_wat, paramst_full_wat] = evolve_backward_euler(paramst, protocol1);
       
       
       % Calculate Conductance/ V_rev
       idx_ramp = { [protocol1.index(2,1):protocol1.index(2,2)], [protocol1.index(6,1):protocol1.index(6,2)]};
       for k = 1:length(idx_ramp)
           
            results = results_nw;
            Vr = results.V_command(idx_ramp{k}) * 1e3; % Ramp voltages in mV
            Vm = results.V_mem(idx_ramp{k}) * 1e3; % Membrane voltage in mV
            Ir = results.I_pip(idx_ramp{k})*1e9; % Ramp currents in nA
            slope_nw(j,k) = (max(Ir) - min(Ir)) / (max(Vr) - min(Vr)) * 1e3; % Slope conductance in nS
            Vrev_nw(j,k) = interp1(Ir, Vr, 0); % Reversal potential in mV
            Vmrev_nw(j,k) = interp1(Ir, Vm, 0); % Reversal potential in mV at membrane
        
            results = results_wat;
            Vr = results.V_command(idx_ramp{k}) * 1e3; % Ramp voltages in mV
            Vm = results.V_mem(idx_ramp{k}) * 1e3; % Membrane voltage in mV
            Ir = results.I_pip(idx_ramp{k})*1e9; % Ramp currents in nA
            slope_wat(j,k) = (max(Ir) - min(Ir)) / (max(Vr) - min(Vr)) * 1e3; % Slope conductance in nS
            Vrev_wat(j,k) = interp1(Ir, Vr, 0); % Reversal potential in mV
            Vmrev_wat(j,k) = interp1(Ir, Vm, 0); % Reversal potential in mV at membrane
        
            %    fprintf('\n Slope Conductance = %.2f nS. Reversal voltage = %.2f mV. Vrev at membrane = %.2f mV.', slope(j,k) , Vrev(j,k), Vmrev(j,k) );
       end
    end
    
    % Add in origin
    slope_nw = [0, 0; slope_nw];
    slope_wat = [0, 0; slope_wat];
    Vt =  -1 * params.R * params.T / params.F *  log(  mem_p(1)  / mem_p(2) ) * 1e3; % membrane reversal potential
    Vrev_nw = [Vt, Vt; Vrev_nw];
    Vrev_wat = [Vt, Vt; Vrev_wat];
    Vmrev_nw = [Vt, Vt; Vmrev_nw];
    Vmrev_wat = [Vt, Vt; Vmrev_wat];
    
    figure(1); clf;
    
    plot(slope_nw(:,1), Vrev_nw(:,2),'k--'); hold on
    plot(slope_wat(:,1), Vrev_wat(:,2),'k'); hold on
    
    plot(slope_nw(:,1), Vmrev_nw(:,2),'g--'); hold on
    plot(slope_wat(:,1), Vmrev_wat(:,2),'g'); hold on
        
    axis([0 115 -77.5 -25])
        set(gca,'YTick',[-75 -60 -45 -30]','TickDir','out');
        set(gca,'XTick',[0 50 100],'TickDir','out');
        xlabel('g_{slope} (nS)');
        ylabel('V_{rev} (mV)'); 
        box off;
    
    save_figure('Fig_Ionic_Selectivity_Water_v1',[ 9 9]);       
            
    function result = I_GHK(params, V, rho_in, rho_out)
        
        idx = 1;
        Am = params.A(params.membrane_N);
        u = params.z(idx)*params.F*V/(params.R*params.T) ; 
        uz = (abs(u)<1e-6); % Uncharged particles and zero membrane voltge
        result = params.F * Am * params.membrane_perm(idx) .* (u + uz) ./ (1 - exp(-u) + uz) .* (rho_in - rho_out.*exp(-u)); 