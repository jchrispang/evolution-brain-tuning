function param = loadParameters_reducedWongWang_func
%% loadParameters_reducedWongWang_func.m     
%
% Contains all the parameters of the reduced Wong-Wang model and the 
% computational parameters for the simulations/calculations. It is
% necessary to make an instance of the function before you can use the 
% parameters.
%
% Example:
% >> param = utils.loadParameters_reducedWongWang_func;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%    (example: param.w = 0.5 if you want to change the recurrent connection strength).
%
% 2. The depdendent parameters need to manually updated when other independent 
%    parameters are changed.
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, University, 2021

%%
    % =====================================================================
    %               DEFAULT INDEPENDENT MODEL PARAMETERS
    % ===================================================================== 
        
    param.A       = [];       % connectivity matrix
        
    % coupling parameters
    param.w       = 0.2;      % recurrent connection strength [unitless]
    param.J       = 0.2609;   % synaptic coupling [nA]
    param.G       = 0.2;      % global scaling constant [unitless]
        
    % input
    param.I       = 0.33;     % excitatory subcortical input [nA]  % previous value 0.3
        
    % firing rate parameters
    param.a       = 270;      % [(V nC)^(-1)]
    param.b       = 108;      % [Hz]
    param.d       = 0.154;    % [s]
        
    % synaptic parameters
    param.tau_s   = 0.1;      % [s]
    param.gamma_s = 0.641;    % [unitless]
        
    % noise
    param.sigma   = 0.003; 
        
    % =====================================================================
    %                    COMPUTATIONAL PARAMETERS
    % ===================================================================== 

    param.tstep   = 0.01;     % time step
    param.tmax    = 1000;     % maximum time
        
    % =====================================================================
    %                     DEPENDENT PARAMETERS
    % ===================================================================== 
        
    param.N = size(param.A, 2);
    param.tspan = [0, param.tmax];
    param.T = 0:param.tstep:param.tmax;
    
end
