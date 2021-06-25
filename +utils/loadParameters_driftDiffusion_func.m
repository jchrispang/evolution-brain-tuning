function param = loadParameters_driftDiffusion_func
%% loadParameters_driftDiffusion_func.m     
%
% Contains all the parameters of the drift diffusion model and the 
% computational parameters for the simulations/calculations. It is
% necessary to make an instance of the function before you can use the 
% parameters.
%
% Example:
% >> param = utils.loadParameters_driftDiffusion_func;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%    (example: param.thres = 0.1 if you want to change the decision threshold).
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
    
    param.beta    = 1;      % drift rate (can be a constant or a column vector)
    param.thres   = 1;      % decision threshold [unitless]
    param.sigma   = 1;      % noise [unitless]
        
    % =====================================================================
    %                    COMPUTATIONAL PARAMETERS
    % ===================================================================== 

    param.tstep   = 0.01;   % time step
    param.tmax    = 100;    % maximum time
        
    % =====================================================================
    %                     DEPENDENT PARAMETERS
    % ===================================================================== 
        
    param.N = size(param.A, 2);
    param.tspan = [0, param.tmax];
    param.T = 0:param.tstep:param.tmax;
    
end
