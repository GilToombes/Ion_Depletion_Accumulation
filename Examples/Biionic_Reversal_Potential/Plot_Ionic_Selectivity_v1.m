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
%              Protocol 2 - Ext Solution 150mM NMDG-Cl only during prolonged activation
%
%                   V_hold = -30mV
%                   Voltage ramps = -30mV to 30mV over 0.5 seconds
%
%                   Initialize for 0.5 seconds at V_hold in 150mM NaCl
%                   Apply agonist and perform Ramp in 150mM NaCl
%                   Wait 0.5 seconds at V_hold in absence of agonist in 150mM NaCl
%                   Apply agonist for T_stim seconds holding at V_stim in 150mM NMDG
%                   Wait 0.5 seconds at V_hold in absence of agonist in 150mM NaCl
%                   Apply agonist and perform Ramp in 150mM NaCl
%                   Wait 0.5 seconds at V_hold in absence of agonist in 150mM NaCl
%
%                   Plot :
%                           Holding voltage versus T
%                           Current vs t
%                           I versus V for ramps
%
%                   For Ramps
%                           Slope Conductance = 178.60 nS. Reversal voltage = -0.76 mv
%                           Slope Conductance = 252.85 nS. Reversal voltage = 10.80 mv

function Fig13_Ionic_Selectivity_v1

	% Configuration parameters
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic_perfusion';  % Option specifying which parameters are variable
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
        %params.membrane_Pf = 5e-5   ; 			% Membrane permeability to water (m/s)
		params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell

	% Initialize all remaining parameters
        params = initialize_system(params);
        params.Cl = params.Cl * 0; % Compensate capacitance.
        params.Ct = params.Ct * 0; % Compensate capacitance.
        
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. ', sum(resistance(params)/1e6));
    fprintf('\n Membrane Permeability = [%e, %e, %e] m/s', mem_p');
    
    paramst = params;
    
    % Now define protocols
    
    % Protocol 1 - Use just NMDG external
    % Protocol 2 - Only have NMDG external during prolonged stimulation
    
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
    
    intervals2.duration=                 [T_rest,    T_ramp,     T_rest,     T_stim,     T_rest,     T_ramp,     T_rest]; % Duration of each of intervals
    intervals2.N_eval =                  [2,         N_ramp,     N_rest,     N_stim,     N_rest,     N_ramp,     N_rest]; % Number of timepoints during each interval
    intervals2.variables.V_hold =        {V_hold2,   V_ramp2,    V_hold2,    V_stim,     V_hold2,    V_ramp2,    V_hold2}; % Holding voltage
    intervals2.variables.membrane_perm=  {0*mem_p,   mem_p,      0*mem_p,    mem_p,      0*mem_p,    mem_p,      0*mem_p};% Membrane Permeability.
    intervals2.variables.rho_bath =      {rho_b,     rho_b,       rho_b,     rho_bstim,   rho_b,     rho_b,       rho_b}; % Bath solution
    
    protocol1 = create_protocol(intervals1); % Create Protocol 
    protocol2 = create_protocol(intervals2); % Create Protocol 
    
    
    [paramst, results] = evolve_backward_euler(params, inf); % Initialize state
    
    % Simulate Protocol 1    
    [paramst, results, paramst_full] = evolve_backward_euler(paramst, protocol1);
        
    %Initiliaze Plot
    figure(1); clf;
            
            t_range = [0, 18]; % Timebased
            
            % Plot Positions
            %     left, bottom,     width,      height
            l1 = 0.12; w1 = 0.4;
            hv = 0.05;
            hc = 0.15;
            hic = 0.15;
            hconc = 0.15;
            hverr = 0.08;
                 
            pos =   [l1,  0.94,     w1,  hv; % Voltage
                     l1,  0.77,     w1,  hc; % Current
                     l1,   0.6,     w1,  hic ; % Individual Currents                 
                     l1   0.42,     w1,  hconc; % Concentration
                     l1,  0.32,     w1,  hverr;  % Voltage Error
                     0.55,  0.70, 0.35, 0.28;  % IV
                     0.65, 0.39, 0.325,   0.25; % Vrev versus slope
                     l1,   0.25 ,   w1,   hv; % Voltage Symmetric Solutions
                     l1,   0.01  ,  w1,  0.235  ; % Current Symmetric Solutions
                     0.55  0.03,   0.35,  0.28;];   % IV Voltage Symmetric
                 
        % Voltage and stimuli
        subplot('Position', pos(1,:));
        plot(protocol1.t, results.V_command*1e3,'k'); hold on
        axis([t_range, -100 -10]); axis off
    
        % ATP Bar
        plot(protocol1.t(protocol1.index(2,:)), -15 + [0,0],'k','Linewidth',1);
        plot(protocol1.t(protocol1.index(4,:)), -15 + [0,0],'k','Linewidth',1);
        plot(protocol1.t(protocol1.index(6,:)), -15 + [0,0],'k','Linewidth',1);
        
        % Net current    
        subplot('Position',pos(2,:));
        plot(protocol1.t, results.I_pip*1e9,'k'); hold on;
        plot([0,18],[0 0], 'k:');
        axis([t_range, -4.1, 6.1]); axis off;
        plot([2 7], [-4 -4],'k','Linewidth',1); % Time bar
        plot( [0.2 0.2], [2 4],'k','Linewidth',1); % Current bar
 
        % Individual membrane currents
        subplot('Position',pos(3,:));
        % Extract membrane currents
        I_Na_mem= zeros(1,length(paramst_full));
        I_NMDG_mem = zeros(1,length(paramst_full));
        I_Na_pip= zeros(1,length(paramst_full));
        I_NMDG_pip = zeros(1,length(paramst_full));
        I_Cl_pip = zeros(1,length(paramst_full));
   
            for k=1:length(paramst_full)
                I_Na_mem(k) = paramst_full(k).J(1,params.membrane_N) * params.F;
                I_Na_pip(k) = paramst_full(k).J(1,params.membrane_N-1) * params.F;
                I_NMDG_mem(k) = paramst_full(k).J(2,params.membrane_N) * params.F;
                I_NMDG_pip(k) = paramst_full(k).J(2,params.membrane_N-1) * params.F;
                I_Cl_pip(k)= -1 * paramst_full(k).J(3,params.membrane_N-1) * params.F;
            end
        plot(protocol1.t, I_Na_mem*1e9,'b'); hold on
        plot(protocol1.t, I_Na_pip*1e9,'b--'); hold on
        plot(protocol1.t, I_NMDG_mem*1e9,'g'); hold on
        plot(protocol1.t, I_NMDG_pip*1e9,'g--'); hold on
        plot(protocol1.t, I_Cl_pip*1e9,'y--'); hold on
        axis([t_range, -7.75 ,10.75]);
        ylabel('I (nA)');
        set(gca,'XTick',[],'YTick',[-5, 0, 5, 10]);
        
    % Concentrations
        subplot('Position',pos(4,:));
        plot(protocol1.t, results.rho_intra(1,:),'b'); hold on
        plot(protocol1.t, results.rho_intra(2,:),'g'); hold on
        plot(protocol1.t, results.rho_intra(3,:)),'k'; hold on
        axis([t_range, 0 , 300]);
        ylabel('\rho_{cell} (mM)');       set(gca,'XTick',[]);
    
    % Voltage difference
        subplot('Position',pos(5,:));
        plot(protocol1.t, (results.V_mem - results.V_command)* 1e3,'k'); hold on
        axis([t_range, -30 , 30]);
        ylabel('\Delta V (mV)');    set(gca,'XTick',[]);
        set(gca,'TickDir','out');
        
    % IV
        subplot('Position', pos(6,:));
        idx_ramp = { [protocol1.index(2,1):protocol1.index(2,2)], [protocol1.index(6,1):protocol1.index(6,2)]};
        col = 'kr';

        % Evaluate each ramp
        for j = 1:length(idx_ramp)
            Vr = results.V_command(idx_ramp{j}) * 1e3; % Ramp voltages in mV
            Ir = results.I_pip(idx_ramp{j})*1e9; % Ramp currents in nA
            slope = (max(Ir) - min(Ir)) / (max(Vr) - min(Vr)) * 1e3; % Slope conductance in nS
            Vrev = interp1(Ir, Vr, 0); % Reversal potential in mV
            fprintf('\n Slope Conductance = %.2f nS. Reversal voltage = %.2f mv', slope , Vrev );
            plot(Vr, Ir,col(j)); hold on
        end
        ax = gca;   ax.XAxisLocation = 'origin';    ax.YAxisLocation = 'origin';
  
        axis([-90 -20 -4 7]);
        ylabel('\Delta V (mV)');
        xlabel('V (mV)');
        ylabel('I (nA)');
        set(gca,'TickDir','out'); 
        box off;
        
    % Add Vrev versus Slope conductance
    sf = [0.001, 0.01, 0.03, 0.05, 0.08, 0.1,0.14, 0.22, 0.3, 0.35, 0.5, 0.65, 0.8, 1, 1.55];
    
    % Calculate for different membrane permeabilities
    for j=1:length(sf)
    
       % Initialize state
       paramst = params; paramst.V_hold = V_hold1 ;  [paramst, results] = evolve_backward_euler(paramst, inf); 
    
       % Change permeability in protocol
       mp = mem_p * sf(j);
       intervals1.variables.membrane_perm=  {0*mp,   mp,      0*mp,    mp,      0*mp,    mp,      0*mp};% Membrane Permeability.
       protocol1 = create_protocol(intervals1); % Create Protocol 
    
       % Simulate Protocol 1    
       [paramst, results, paramst_full] = evolve_backward_euler(paramst, protocol1);

       % Calculate Conductance/ V_rev
       idx_ramp = { [protocol1.index(2,1):protocol1.index(2,2)], [protocol1.index(6,1):protocol1.index(6,2)]};
       for k = 1:length(idx_ramp)
            Vr = results.V_command(idx_ramp{k}) * 1e3; % Ramp voltages in mV
            Vm = results.V_mem(idx_ramp{k}) * 1e3; % Membrane voltage in mV
            Ir = results.I_pip(idx_ramp{k})*1e9; % Ramp currents in nA
            slope(j,k) = (max(Ir) - min(Ir)) / (max(Vr) - min(Vr)) * 1e3; % Slope conductance in nS
            Vrev(j,k) = interp1(Ir, Vr, 0); % Reversal potential in mV
            Vmrev(j,k) = interp1(Ir, Vm, 0); % Reversal potential in mV at membrane
        %    fprintf('\n Slope Conductance = %.2f nS. Reversal voltage = %.2f mV. Vrev at membrane = %.2f mV.', slope(j,k) , Vrev(j,k), Vmrev(j,k) );
       end
    end
    
    % Add in origin
    slope = [0, 0; slope];
    Vt =  -1 * params.R * params.T / params.F *  log(  mem_p(1)  / mem_p(2) ) * 1e3; % membrane reversal potential
    Vrev = [Vt, Vt; Vrev];
    Vmrev = [Vt, Vt; Vmrev];
    
    subplot('Position', pos(7,:));
        plot(slope(:,1), Vrev(:,2),'k'); hold on
        plot(slope(:,1), Vmrev(:,2),'k--');
        axis([0 115 -77.5 -28])
        set(gca,'YTick',[-75 -60 -45 -30]','TickDir','out');
        xlabel('\sigma_{slope} (nS)');
        ylabel('V_{rev} (mV)'); 
        box off;
        
    % Now add protocol using symmetric solutions to detect concentration change.
    [paramst, results] = evolve_backward_euler(params, inf); % Initialize state
    [paramst, results, paramst_full] = evolve_backward_euler(paramst, protocol2);
        
    % Voltage and stimuli
    subplot('Position', pos(8,:));
       plot(protocol1.t, results.V_command*1e3,'k'); hold on
       axis([t_range, -61 30]); axis off
    
        % ATP Bar
        plot(protocol1.t(protocol1.index(2,:)), -15 + [0,0],'k','Linewidth',1);
        plot(protocol1.t(protocol1.index(4,:)), -15 + [0,0],'k','Linewidth',1);
        plot(protocol1.t(protocol1.index(6,:)), -15 + [0,0],'k','Linewidth',1);
    
    % Sodium - NMDG
        plot(protocol1.t(protocol1.index(4,:)), 15 + [0,0],'b','Linewidth',1);
     
    % Net current    
    subplot('Position',pos(9,:));
        plot(protocol1.t, results.I_pip*1e9,'k'); hold on;
        plot([0,18],[0 0], 'k:');
        axis([t_range, -9.9, 6.1]); axis off;
        plot([2 7], [-4 -4],'k','Linewidth',1); % Time bar
        plot( [0.2 0.2], [2 4],'k','Linewidth',1); % Current bar
    
    % IV
        subplot('Position', pos(10,:));
        idx_ramp = { [protocol1.index(2,1):protocol1.index(2,2)], [protocol1.index(6,1):protocol1.index(6,2)]};
        col = 'kr';
        fprintf('\n\n Symmetric solutions\n');
        
        % Evaluate each ramp
        for j = 1:length(idx_ramp)
            Vr = results.V_command(idx_ramp{j}) * 1e3; % Ramp voltages in mV
            Ir = results.I_pip(idx_ramp{j})*1e9; % Ramp currents in nA
            slope = (max(Ir) - min(Ir)) / (max(Vr) - min(Vr)) * 1e3; % Slope conductance in nS
            Vrev = interp1(Ir, Vr, 0); % Reversal potential in mV
            fprintf('\n Slope Conductance = %.2f nS. Reversal voltage = %.2f mV', slope , Vrev );
            plot(Vr, Ir,col(j)); hold on
        end
        
        ax = gca;   ax.XAxisLocation = 'origin';    ax.YAxisLocation = 'origin';
        axis([-30 30 -6 6]);
        ylabel('\Delta V (mV)');
        xlabel('V (mV)');
        ylabel('I (nA)');
        set(gca,'TickDir','out'); box off;
    
        
    
	save_figure('Fig_WholeCell_Biionic_Vrev_Shift1',[ 9 14]);       
            
    function result = I_GHK(params, V, rho_in, rho_out)
        
        idx = 1;
        Am = params.A(params.membrane_N);
        u = params.z(idx)*params.F*V/(params.R*params.T) ; 
        uz = (abs(u)<1e-6); % Uncharged particles and zero membrane voltge
        result = params.F * Am * params.membrane_perm(idx) .* (u + uz) ./ (1 - exp(-u) + uz) .* (rho_in - rho_out.*exp(-u)); 