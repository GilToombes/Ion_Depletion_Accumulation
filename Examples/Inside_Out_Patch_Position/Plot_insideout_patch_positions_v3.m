% function plot_insideOut_patch_positions
%
%   Show effect of membrnae patch position in inside out configuration
%
%   Pipette   - R = 5 MOhms with 150mM KCl
%   
%   Membrane patch located at 
%   Plot concentration and voltage-chang at both sides of membrg
%    as a function of IK+*Raccess, and |INa+| * Raccess


function plot_insideout_patch_positions

	% Configuration parameters
		params.config = 'inside-out';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	
    % Solution Parameters
        params.D =  [1,  1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ and Cl-.
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
		params.Nseg  =  21; 		% Number of segments between the pipette electrode and tip
		params.V_hold = 0e-3; 	% Initial holding voltage
        params.Rseries = 4.5e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
	% Cell and Membrane parameters
        ltip = params.r_tip / sin(params.theta/2);
        params.membrane_x = ltip;   			% Position of membrane patch back from the tip (in m)
		params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0 ]'; 		% Membrane permability to each ion in m/s

    % Calculate membrane permeability for an initial current of 2nA K+
    paramst = initialize_system(params);
    I_mem = [2e-9, 0 ]'; % Initial current of 2nA by K+
    V_mem = 100e-3; % 100mV 
    P_GHK = calculate_PGHK_from_I(I_mem, paramst.rho_pip, paramst.rho_bath, V_mem, paramst.z, paramst.A(paramst.membrane_N), paramst.T);
    
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(paramst)/1e6), paramst.Rseries/1e6)
    
    % Report Access Resistance and Series Resistance Correction
   	%fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', resistance(params)/1e6, params.Rseries/1e6)

    
     
     intervals.duration =                [50e-3,    1e-5,   9.9e-4, 9e-3,   19e-2,  1e-5, 1e-4, 1e-3,   4e-3,   45e-3 ] ; % Duration of each of intervals
     intervals.N_eval =                  [2,        40,     10,     10,     20,     10,   10,   10,     10,     20] ; % Number of timepoints
    %intervals.variables.V_hold =        {0e-3,     V(j),   V(j),   V(j),   V(j) ,  0e-3, 0e-3, 0e-3,   0e-3,   0e-3};

    
    % Plot SS concentration in Patch Pipette
    mem_pos = 3e-6;  
    paramst = params;
    paramst.Nseg = 60;
    paramst.membrane_x = mem_pos;  			% Position of membrane patch back from the tip (in m)
    paramst.membrane_perm = P_GHK / (1/2 + mem_pos/(2*ltip))^2 /  1.7153;  % 1 nA current
    paramst.V_hold = V_mem;
    paramst.Rseries = 0;       % Series resistance correction in Ohms
    paramst = initialize_system(paramst);
    [paramst, results] = evolve_backward_euler(paramst, inf);
    
    figure(1); clf;
    panel_bottom = [ 0.7, 0.4, 0.1];
    panel_heights = [0.25, 0.25, 0.25];
    pl = 0.1; % Left postion of panels
    pw = 0.8 ; % Panel width
    xrange = [-1 8];
    
    x = paramst.x*1e6;
    Nm = paramst.membrane_N;
    Ns = paramst.Nseg;
    x_intra = [-2, x([1:Nm])];
    rho_intra = [paramst.rho_bath(1), paramst.rho(1,[1:Nm]) ];
    Vf = paramst.V(Ns);
    
    V_intra = ([0, paramst.V([1:Nm])] -  Vf) * 1e3;
        
    idx_e = [(Nm+1):Ns];
    x_extra =  x(idx_e);
    rho_extra = paramst.rho(1,idx_e);
    V_extra = (paramst.V(idx_e) -  Vf)* 1e3;

    % Plot Pipette and Membrane Shape
    idx = 1;
    subplot('position',[ pl, panel_bottom(idx), pw, panel_heights(idx)]);
        x = [0 max(xrange)];
        y = params.r_tip * 1e6 + x * sin(params.theta/2);
        xmem = paramst.membrane_x * 1e6;
        ymem = params.r_tip * 1e6 + xmem * sin(params.theta/2);
        
        plot(x,y,'k'); hold on
        plot(x,-1 * y, 'k');
        plot( xmem + [0 0], ymem * [-1 1],'b');         
        axis equal
        axis([xrange, -1.3, 1.3]);
        axis off
    
    % Plot Voltage
    idx = 2;
    subplot('position',[ pl, panel_bottom(idx), pw, panel_heights(idx)]);
    Vsub = 90;
    plot(x_intra, V_intra-Vsub,'k'); hold on
    plot(x_extra, V_extra,'k'); hold on
    YTickLabel = [0 1 2  98 99 100];
    YTick =      [0 1 2  98-Vsub, 99-Vsub, 100-Vsub];
    set(gca,'YTick',YTick,'YTickLabel',YTickLabel);
    axis([xrange, 0 10.5]);
    box off
    set(gca,'XTick',[]);
    ylabel('V (mV)')
    
    % Plot Concentration
    idx = 3;
    subplot('position',[ pl, panel_bottom(idx), pw, panel_heights(idx)]);
    plot(x_intra, rho_intra,'b','LineWidth',1); hold on
    plot(x_intra, rho_intra,'LineWidth',1,'Color',[0.5 0.5 0.5],'LineStyle','--'); hold on
    plot(x_extra, rho_extra,'b','LineWidth',1); hold on
    plot(x_extra, rho_extra,'LineWidth',1,'Color',[0.5 0.5 0.5],'LineStyle','--'); hold on
    axis([xrange, 132 172]); box off
    ylabel('\rho (mM)');
    xlabel('x (\mum)');
    
 	save_figure('Fig_InsideOut_patch_static_v3',[ 15 15]);  
    
        
    figure(1); clf;

    V = [-100:25:100] * 1e-3;
    mem_pos = [1e-6, 3e-6, 10e-6];

    for j=1:length(mem_pos)
        for k =1:length(V)
        
            paramst = params;
            paramst.membrane_x = mem_pos(j);  			% Position of membrane patch back from the tip (in m)
    
            %mem_perm =  1.09e-4 * [ 1 0 ]'; % Membrane permeability
            %p = mem_perm / (1 + mem_pos(j)/ltip)^2;    
            p = P_GHK / (1/2 + mem_pos(j)/(2*ltip))^2;    
         
            intervals.variables.membrane_perm = {p,         p,      p,      p,      p,      p,      p,  p,      p,      p };
            intervals.variables.V_hold =        {0e-3,     V(k),   V(k),   V(k),   V(k) ,  0e-3, 0e-3, 0e-3,   0e-3,   0e-3};

            protocol = create_protocol( intervals );
    
            paramst.V_hold = V(k);
            paramst = initialize_system(paramst);
            paramst.Ct = 0 * paramst.Cl;
            paramst.Cl = 0 * paramst.Cl;
            [paramst, results] = evolve_backward_euler(paramst, inf);
            [paramst, results] = evolve_backward_euler(paramst, protocol);

            subplot(3,3,j);
                hold on; plot(protocol.t, results.I_pip*1e9,'k'); hold on;

            subplot(3,3,j+3);
                hold on; plot(protocol.t, results.rho_intra(1,:),'k'); hold on;
                
            subplot(3,3,j+6);
                hold on; plot(protocol.t, results.rho_extra(1,:),'k'); hold on;
        end
           fprintf(1,'\n\n Membrane position = %.2f microns, P_K+ = %e m/s', mem_pos(j) * 1e6 , p(1));
    end
    
    subplot(3,3,1)
            axis([0 0.3 -2.1 2.1])
            ylabel('I (nA)');
            title('x = 1\mu m')
            
        subplot(3,3,2)
            axis([0 0.3 -2.1 2.1])
            title('x = 3 \mu m')
        
        subplot(3,3,3)
            axis([0 0.3 -2.1 2.1])
            title('x = 10 \mu m')
    
        subplot(3,3,4)
            ylabel('\rho_{K+,intra} (mM)');
            axis([0 0.3 105 195])
    
        subplot(3,3,5)
            axis([0 0.3 105 195])
    
        subplot(3,3,6)
            axis([0 0.3 105 195])
    
        subplot(3,3,7)
            axis([0 0.3 105 195])
            xlabel('t (s)');
            ylabel('\rho_{K+,extra} (mM)');
            
        subplot(3,3,8)
            axis([0 0.3 105 195])
            xlabel('t (s)');
        
        subplot(3,3,9)
            axis([0 0.3 105 195])
            xlabel('t (s)');
            
	save_figure('Fig_InsideOut_patch_position_v3',[ 15 15]);       
            
    