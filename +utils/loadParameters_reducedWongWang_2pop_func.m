function param = loadParameters_reducedWongWang_2pop_func
%% loadParameters_reducedWongWang_2pop_func.m     
%
% Contains all the parameters of the reduced Wong-Wang model with one 
% excitatory/inhibitory population and computational parameters for the 
% simulations/calculations. It is necessary to make an instance of the 
% function before you can use the 
% parameters.
%
% Example:
% >> param = utils.loadParameters_reducedWongWang_2pop_func;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%   (example: param.w_E = 0.2 if you want to change the excitatory recurrent connection strength).
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
        
    param.A       = [];     % connectivity matrix
        
    % coupling parameters
    param.w_E     = 1;      % excitatory scaling of external input [unitless]
    param.w_I     = 0.7;    % inhibitory scaling of external input [unitless]
    param.w_p     = 1.4;    % recurrent excitation weight [unitless]
    param.J_E     = 0.15;   % excitatory synaptic coupling [nA]
    param.J_I     = 1;      % inhibitory synaptic coupling [nA]
    param.G       = 0.2;    % global scaling constant [unitless]
    param.lambda  = 0;      % long-range feedforward inhibition [unitless]

    % input
    param.I0      = 0.42;   % effective external input [nA]  % previous value 0.382
    param.Iext    = 0;      % external stimulation [nA]
                            % 0 for resting state

    % excitatory firing rate parameters
    param.a_E     = 310;    % [(V nC)^(-1)]
    param.b_E     = 125;    % [Hz]
    param.d_E     = 0.16;   % [s]

    % inhibitory firing rate parameters
    param.a_I     = 615;    % [(V nC)^(-1)]
    param.b_I     = 177;    % [Hz]
    param.d_I     = 0.087;  % [s]

    % synaptic parameters
    param.tau_E   = 0.1;    % [s] (also equal to tau_NMDA)
    param.tau_I   = 0.01;   % [s] (also equal to tau_GABA)
    param.gamma   = 0.641;  % [unitless]

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
