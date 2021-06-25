function H = calc_firingRate_reducedWongWang(param, S)
% calc_firingRate_reducedWongWang.m
%
% Calculates the firing rate of the reduced Wong-Wang model with one excitatory population
%
% Inputs: param : parameters of the model
%         S     : synaptic gating variable [NxT]
%                 T = time length of simulation
%
% Output: H     : firing rate [NxT]
%                 T = time length of simulation
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, University, 2021

%%

x = param.w*param.J*S + param.G*param.J*param.A*S + param.I;
H = (param.a*x - param.b)./(1 - exp(-param.d*(param.a*x - param.b)));