function param = loadParameters_WilsonCowan_func
%% loadParameters_WilsonCowan_func.m     
%
% Contains all the parameters of the Wilson-Cowan model and the 
% computational parameters for the simulations/calculations. It is
% necessary to make an instance of the function before you can use the 
% parameters.
%
% Example:
% >> param = utils.loadParameters_WilsonCowan_func;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%   (example: param.G = 2 if you want to change the global coupling strength).
%
% 2. The dependent parameters need to be manually updated when other independent 
%    parameters are changed.
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%
    % =====================================================================
    %               DEFAULT INDEPENDENT MODEL PARAMETERS
    % ===================================================================== 
        
    param.A         = [];      % connectivity matrix
        
    % coupling parameters
    param.w_EE      = 16;      % excitatory to excitatory coupling [unitless]
    param.w_EI      = 12;      % inhibitory to excitatory coupling [unitless]
    param.w_IE      = 15;      % excitatory to inhibitory coupling [unitless]
    param.w_II      = 3;       % inhibitory to inhibitory coupling [unitless]
    param.G         = 1;       % global coupling strength [unitless]

    % input
    param.G_E       = 0.5;       % scaling constant of excitatory drive [unitless]
    param.G_I       = 0;       % scaling constant of inhibitory drive [unitless]

    % firing rate parameters
    param.a_E       = 1.5;     % excitatory gain
    param.a_I       = 1.5;     % inhibitory gain
    param.mu_E      = 3;       % excitatory firing threshold
    param.mu_I      = 3;       % inhibitory firing threshold

    % synaptic parameters
    param.tau_E     = 2.5e-3;  % excitatory time constant [s]
    param.tau_I     = 3.75e-3; % inhibitory time constant [s]
    param.v         = 10;      % propagation speed [m/s]

    % noise
    param.sigma_E   = 5e-5;    % excitatory noise strength
    param.sigma_I   = 5e-5;    % inhibitory noise strength
        
    % =====================================================================
    %                    COMPUTATIONAL PARAMETERS
    % ===================================================================== 

    param.tstep     = 0.001;   % time step
    param.tmax      = 100;     % maximum time
        
    % =====================================================================
    %                     DEPENDENT PARAMETERS
    % ===================================================================== 
        
    param.N = size(param.A, 2);
    param.tspan = [0, param.tmax];
    param.T = 0:param.tstep:param.tmax;
    
end
