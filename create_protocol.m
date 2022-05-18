% function protocol = create_protocol( intervals )
% 
%   Create an experimental protocol
%	
%   Input -    
%		intervals.duration = [t1, t2, t3, ..., tp]  duration of each of P intervals in seconds
%		intervals.N_eval =   [N1, N2, N3, ..., Np]  number of time points during each interval
%		
%       Optional 
%       intervals.variables 
%            set of fields which change between or during intervals
%            such as V_hold, rho_bath, or membrane.perm
%            e.g. intervals.variables.V_hold = {-60e-3, [-100e-3,100e-3],-60e-3}
%                               1st interval - hold voltage at -60mV
%                               2nd interval - ramp from -100mV to 100mV
%                               3rd interval - hold voltage at -60mV
%
%   Output -
%       protocol.t -> Time value for each step
%       protocol.index -> indices for each interval 
%                         e.g. j-th interval corresponds to time values from 
%                                   protocol.index(j,1) to
%                                   protocol.index(j,2)
%       protocol.variables. -> Value for each variable during each timestep
%

function protocol = create_protocol( intervals )

    % If no input, create a simple 3-interval protocol 
    if (nargin<1)|~isfield(intervals,'duration') 
            intervals.duration = [0.1, 0.5, 0.1]; % Interval duration in seconds 
            intervals.N_eval =   [10 , 50,  10]; % Number of time points during each interval
    end
    
    N_int = length(intervals.duration); % Number of intervals
    Ne = sum(intervals.N_eval) + N_int; % Number of time points
    
    protocol.t = zeros(1,Ne); % Time points 
    protocol.index = zeros(N_int,2); % Index for each interval
    
    % Initialize output variables
    if isfield(intervals, 'variables')
        % v_name = fields(intervals.variables);
        v_name = fieldnames(intervals.variables);
         Nv = length(v_name); % Number of variables
        for k=1:Nv
            vs = size(intervals.variables.(v_name{k}){1},1);
            protocol.variables.(v_name{k}) = zeros(vs, Ne);
        end
    else
         Nv = 0; % No varaibles which change during intervals
    end
    
    
    % Process each interval in turn
    t_cur = 0; % Current time.
    idx = 1;
    for j=1:N_int
        
        Ns = intervals.N_eval(j); % Number of steps in interval;
        isi = [0:Ns]; % Index within step
        isi_frac = isi/Ns; isi_frac(1) = 0.01/Ns; % Fractional values with first time point offset by 1% 
        
        protocol.t(idx + isi) = t_cur + isi_frac * intervals.duration(j); % Time values during interval
        protocol.index(j,:) = [idx, idx+Ns]; % Starting and finishing indices of interval
        
        % Calculate values for each variable during step
        for k=1:Nv
            
            val= intervals.variables.(v_name{k}){j}; % Value of variable on step
            if (size(val,2)==1) % Change of variable during step
                    val(:,2) = val(:,1); 
            end
            protocol.variables.(v_name{k})(:,idx+isi) = val(:,1) * (1-isi_frac) + val(:,2)*isi_frac;
        end
        
        t_cur = t_cur + intervals.duration(j);
        idx = idx + Ns + 1;
     end

