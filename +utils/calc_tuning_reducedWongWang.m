function [Smean, Hmean, Hmax] = calc_tuning_reducedWongWang(param, w_vec, tpre, num_trials)
% calc_tuning_reducedWongWang.m
%
% Calculate tuning curve from reduced Wong-Wang model with one excitatory
% population with respect to recurrent connection strength w
%
% Inputs: param      : parameters of the model
%         w_vec      : recurrent connection strengths (vector)
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
Smean = zeros(param.N, length(w_vec), num_trials);
Hmean = zeros(param.N, length(w_vec), num_trials);
Hmax = zeros(param.N, length(w_vec), num_trials);

% simulation
for w_ind = 1:length(w_vec)
    param.w = w_vec(w_ind);
    
    for trial = 1:num_trials
        sol = models.reducedWongWang(param);
        S = sol.y;                                             % synaptic gating
        H = utils.calc_firingRate_reducedWongWang(param, S);   % firing rate

        Smean(:,w_ind,trial) = mean(S(:,time_steady_ind:end),2);
        Hmean(:,w_ind,trial) = mean(H(:,time_steady_ind:end),2);
        Hmax(:,w_ind,trial) = max(H(:,time_steady_ind:end),[],2);
    end
end

