function [H_E, H_I] = calc_activationFunction_WilsonCowan(param, S_E, S_I, P_E, P_I)
% calc_activationFunction_WilsonCowan.m
%
% Calculates the activation_function of the Wilson Cowan model
%
% Inputs: param : parameters of the model
%         S_E   : excitatory firing rate variable [NxT]
%         S_I   : inhibitory firing rate variable [NxT]
%                 T = time length of simulation
%
% Outputs: H_E  : excitatory activation_function [NxT]
%          H_I  : inhibitory activation_function [NxT]
%                 T = time length of simulation
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%

x_E = param.w_EE*S_E - param.w_EI*S_I + param.G*param.A*S_E + param.G_E*P_E;
x_I = param.w_IE*S_E - param.w_II*S_I + param.G_I*P_I;

H_E = activation_function(x_E, param.a_E, param.mu_E);
H_I = activation_function(x_I, param.a_I, param.mu_I);

end

function H = activation_function(x, a, mu)

    H = 1./(1 + exp(-a*(x - mu)));
end