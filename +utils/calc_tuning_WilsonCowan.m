function [SEmean, SImean] = calc_tuning_WilsonCowan(param, wEE_vec, tpre, num_trials)
% calc_tuning_WilsonCowan.m
%
% Calculate tuning curve from Wilson-Cowan model with respect to
% excitatory to excitatory connection strength wEE
%
% Inputs: param      : parameters of the model
%         wEE_vec    : excitatory to excitatory recurrent connection strengths (vector)
%         tpre       : burn time (float)
%         num_trials : number of trials (integer)
%
% Outputs: Smean     : mean synaptic gating
%          Hmean     : mean firing rate
%          Hmax      : maximum firing rate
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%

% calculating index of start of steady state
time_steady_ind = dsearchn(param.T', tpre)+1;

% initialization
SEmean = zeros(param.N, length(wEE_vec), num_trials);
SImean = zeros(param.N, length(wEE_vec), num_trials);

% simulation
for w_ind = 1:length(wEE_vec)
    param.w_EE = wEE_vec(w_ind);
    
    for trial = 1:num_trials
        sol = models.WilsonCowan(param);

        SEmean(:,w_ind,trial) = mean(sol.y_E(:,time_steady_ind:end),2);
        SImean(:,w_ind,trial) = mean(sol.y_I(:,time_steady_ind:end),2);
    end
end

