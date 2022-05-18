% function [params_f, results, params_full] = evolve_backward_euler(params, protocol)
%
% 	Evolve system state using the backward Euler method
%
% 	Inputs
%		params -> system state
%       Optional protocol
%           no protocol -> evaluate steady-state
%           protocol = vector -> set of time points
%           protocol = structure
%               protocol.t ->   vector of M time points
%               protocol.variables -> structure containing variables that
%                                     change during the protocol
%                                     e.g. V_hold, rho_bath, etc.
% 	Ouput
%		params_f    -> final state
%		results.I_pip -> electrical current in pipette electrode at each time point
%		results.rho_intra -> concentrations at intra-cellular side of membrane for each time point
%		results.rho_extra -> concentrations at extra-cellular side of membrane for each time point
%		results.V_mem     -> membrane voltage at each time point
%		results.V_command -> command voltage at each time point (useful if doing series resistance compenstation
%
%       params_full -> vector of state at each point.

function [params_f, results, params_full ] = evolve_backward_euler(params, protocol)

   % Information about parameters
   Np = params.Nseg;    % Number of elements
   Ns = length(params.z); % Number of species
   Nm = params.membrane_N; % Position of membrane
   Nbath =  params.n_bath ; % Bath electode
   Npip = params.n_pip; % Pipette electrode
   Npip_current = max(Npip-1,1); % Index of current fluxes at pipette end     
   if isfield(params,'membrane_func')
            Nm = params.membrane_N; % Index of membrane        
   else
           Nm = 1;
   end
   params_i = params; % Initial parameters
   params_f = params; % Final parameters

   % Protocol Info
   Nvar = 0; % Number of variables that change during protocol is 0 by default
   Nsteps = 1; % Number of time steps to evaluate
   if (nargin<2) % Determine steady state if no protocol specified
        t_val = inf; % Determine steady state
   else
        if isstruct(protocol) % Full Protocol
            t_val = protocol.t; % Time points
            Nsteps = length(t_val); % Number of time points
        
            % Variables which change during timesteps
            if isfield(protocol,'variables')
                v_names = fieldnames(protocol.variables); % Name of variables
                Nvar = length(v_names); % Number of variables
            end
        else
            t_val = protocol; % Protocol is just set of time points
            Nsteps = length(t_val);
        end
   end
    
	% Initialize Output structures
    results.I_pip = zeros(1,Nsteps);  % Electrical current flowing at pipette electrode
	results.rho_intra = zeros(Ns,Nsteps); % Density at each time point
	results.rho_extra = zeros(Ns,Nsteps); % Density at each time point
	results.V_mem = zeros(1,Nsteps);    % Membrane Voltage at each point
    results.V_command = zeros(1, Nsteps); % Command Voltage at each point
	if (nargout>2) ; params_full(Nsteps) = params; end % Full reporting structure
        
    
    % Take a series of time steps.
    for k=1:Nsteps  

		params_i = params;% Initial State
		    
        % Advance time and set new boundary conditions
        if (t_val ~= inf); 
            
            % Standard time step
            dt = t_val(k) - params.t; % Time step
            
            % Adjust all variables that changed during timestep
            params.t = t_val(k); 
            for m=1:Nvar
                params.(v_names{m}) = protocol.variables.(v_names{m})(:,k);
            end
        else
            % Determine equilibrium state
            dt = t_val;
        end

		% Refine final state until all parameters are within cut-off error or exceed Niter_max
        % Use simple Newton-Raphson method on each step
        j = 0; Niter_max = 10; 
        err = All_constraints(params, params_i, dt); % determine initial error
        err_worst = max(abs(err));

        while (j<Niter_max)&(err_worst>1)
                       
                  % Calculate the current Hessian matrix
                  H = Hessian_General(params, params_i, dt);
                  
                  % Precondition the Hessian to improve condition number
                  % prior to matrix division
                  lc = max(abs(H),[],2);   Hl = (diag (1./lc) * H);
                  rc = max(abs(Hl),[],1);     Hlr = Hl  * diag(1./rc);
                  %fprintf('\n cond(H) = %e, cond(Hl) = %e, cond(Hlr) = %e', cond(H), cond(Hl), cond(Hlr));
                  
                  % Calculate the change in each parameter
                  dp = ( (err * diag(1./rc)) / Hlr ) * diag(1./lc);
                  %fprintf('\n%d %e %e',[[1:params.Constraints.N]; err; dp]);
                  %fprintf('%e %e \n',max(abs(err)), max(abs(dp)) );
        
                  % Update all the parameters
                  for m=1:length(params.Constraints.param)
                        pname = params.Constraints.param{m};
                        idx_v = params.Constraints.idx_v{m};
                        idx_c = params.Constraints.idx_c{m};
                        for n=1:length(idx_v)
                            params.(pname)(idx_v(n)) = params.(pname)(idx_v(n))- dp(idx_c(n));
                        end
                  end

                  % Determine the new errors
                  err = All_constraints(params, params_i, dt);
                  err_worst = max( abs(err) ); j = j+1;		
        end

        % Report if did not converge
        if (err_worst > 1)
            fprintf('\nStep size did not converge.\n');
            [err_worst, idx] = max( abs(err));
            fprintf('\nWorst Error Index = %d, Error = %f.', idx, err_worst);
        end
        
		% Calculate current and store parameters for this step
        results.I_pip(k) = params.F * sum(params.J(:,Npip_current) .* params.z); % Current flowing at pipette electrode
		results.V_mem(k) = params.V(Nm) - params.V(Nm + 1); % Membrane voltage
        results.V_command(k) = params.V(Npip) * (1 - 2*(Npip>1.5));
        results.rho_intra(:,k) = params.rho(:,Nm); % Intra-cellular Concentrations
        results.rho_extra(:,k) = params.rho(:,Nm+1); % Extra-cellular Concentrations
        
        % If desired, give full state at each timepoint
        if (nargout>2)
            params_full(k) = params;
        end
        
    end

   params_f = params;

function [result, full_constraints] = All_constraints(params, params_i, delta_t)
    % Return a vector reporting all constraints
    % Inputs
    %       params ->   values of parameters at end of timestep
    %       params_i -> values of parameters at the start of the timestep.
    %       delta_t ->  duration of timestep in seconds.
    %
    % Outputs
    %       result -> Vector reporting value of all constraints
    %       full_constraints -> Cell array reporting value of each
    %                           constraint function 
    %                           (primarily of use for debugging constraint
    %                           functions)
       
    Nc = params.Constraints.N; % Number of constraint functions
    result = zeros(1, Nc); % Vector containing contraint values.  Correct final state gives result(j) = 0.
    
    % Process each contraint function in turn
    for j = 1 :length(params.Constraints.param)
        
        p_name =  params.Constraints.param{j}; % Parameter name
        c_func =  params.Constraints.fhandle{j}; % Constraint function for that parameter
        idx_v = params.Constraints.idx_v{j}; % Indices of parameters that are constrained
        idx_c = params.Constraints.idx_c{j}; % Corresponding indices in constraint vector
        err_sf = params.Constraints.err_tol{j}; % Scale error by error tolerance so error < 1 is acceptable.
        
        Cv = c_func(params, params_i, delta_t); % Evaluate constraint function
        result(idx_c) = Cv(idx_v) / err_sf;       
        full_constraints{j} = Cv;
    end
    
function result = Hessian_General(params, params_i, delta_t)
    % Numerically evaluate the Hessian Matrix for all free parameters
    % Inputs
    %       params ->   values of parameters at end of timestep
    %       params_i -> values of parameters at the start of the timestep.
    %       delta_t ->  duration of timestep in seconds.
    %
    % Outputs
    %       result ->  Nc * Nc matrix with Hessian
    
    
    Nc = params.Constraints.N; % Number of constraints
    result = zeros(Nc, Nc);
    
    base_error = All_constraints(params, params_i, delta_t);
    params_test = params;

    % Process each contraint in turn
    for j = 1 :length(params.Constraints.param)
        
        p_name =  params.Constraints.param{j}; % Parameter name
        %c_func =  params.Constraints.fhandle{j}; % Constraint function for that parameter
        idx_v = params.Constraints.idx_v{j}; % Indices of parameters that are constrained
        idx_c = params.Constraints.idx_c{j}; % Corresponding indices in constraint vector
        sv = params.Constraints.stepsize{j}; % Step size for evaluating derivatives.

        % Process each element of parameter in turn
        for k=1:length(idx_v)
            params_test.(p_name)(idx_v(k)) = params_test.(p_name)(idx_v(k)) + sv ; 
            result(idx_c(k),:) = (All_constraints(params_test, params_i, delta_t) - base_error)/sv;
            params_test.(p_name) = params.(p_name);
        end
    end


function Test_Hessian(H, params, params_i, delta_t)
    % Test function used to verify if a Hessian routine is working
    % correclty
    %
    % Inputs
    %       H -> Nc * Nc array with Hessian
    %       params ->   values of parameters at end of timestep
    %       params_i -> values of parameters at the start of the timestep.
    %       delta_t ->  duration of timestep in seconds.
    

    % Check each parameter which is contrained
    for j=1:length(params.Constraints.param)
   
        % Parameter and indices
        pname = params.Constraints.param{j};
        idx_v = params.Constraints.idx_v{j};
        idx_c = params.Constraints.idx_c{j};
 
        % Vary each index a little
        params_t = params;
        params_t.(pname)(idx_v) = params.(pname)(idx_v) + rand(size(idx_v)) * params.Constraints.stepsize{j};
   
        % Calculate the effect on the error
        err = All_constraints(params, params_i, delta_t);
        err_n = All_constraints(params_t, params_i, delta_t);
        de = err_n - err;

        % Now compute Hessian version
        u = zeros(1, params.Constraints.N);
        u(idx_c) = params_t.(pname)(idx_v) - params.(pname)(idx_v);
        He = u * H;
    
        % Report result
        fprintf('\nParameter - %s\n',pname);
        fprintf('Exact, Hessian, Diff \n',pname);
        fprintf('%e , % e,  %e\n', [de; He; He-de])
   
    end