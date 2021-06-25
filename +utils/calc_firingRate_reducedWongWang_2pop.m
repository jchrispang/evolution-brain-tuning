function [H_E, H_I] = calc_firingRate_reducedWongWang_2pop(param, S_E, S_I)
% calc_firingRate_reducedWongWang_2pop.m
%
% Calculates the firing rate of the reduced Wong-Wang model with one 
% excitatory/inhibitory population
%
% Inputs: param : parameters of the model
%         S_E   : excitatory synaptic gating variable [NxT]
%         S_I   : inhibitory synaptic gating variable [NxT]
%                 T = time length of simulation
%
% Outputs: H_E  : excitatory firing rate [NxT]
%          H_I  : inhibitory firing rate [NxT]
%                 T = time length of simulation
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%

x_E = param.w_E*param.I0 + param.w_p*param.J_E.*S_E + ...
      param.G*param.J_E*param.A*S_E - param.J_I*S_I + param.Iext;
x_I = param.w_I*param.I0 + param.J_E.*S_E + ...
      param.lambda*param.G*param.J_E*param.A*S_E - S_I;  
    
H_E = firing_rate(x_E, param.a_E, param.b_E, param.d_E);
H_I = firing_rate(x_I, param.a_I, param.b_I, param.d_I);
  
end

function H = firing_rate(x, a, b, d)

    H = (a*x - b)./(1 - exp(-d*(a*x - b)));
end