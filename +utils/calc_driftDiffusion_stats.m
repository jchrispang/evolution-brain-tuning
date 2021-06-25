function stats = calc_driftDiffusion_stats(y, T, thres)
% calc_driftDiffusion_stats.m
%
% Calculate various statistics of the drift diffusion model.
%
% Inputs: y      : matrix of evidence [N x T x num_trials]
%         T      : time [vector]
%         thres  : decision threshold (float)
%
% Output: stats  : statistics (struct)
%
% Original: James Pang, QIMR Berghofer, 2020

%%

N = size(y,1);
num_trials = size(y,3);

% time to make decision
[~,decision_ind] = max(abs(y)==thres,[],2);
stats.t_decision = squeeze(T(decision_ind));

% decisions across trials
stats.Pp = sum(y>=thres,3)/num_trials;
stats.Pn = sum(y<=-thres,3)/num_trials;
stats.Pu = sum(abs(y)<thres,3)/num_trials;

% decision across nodes
% stats.Np = squeeze(sum(y>=thres,1)/N)';
% stats.Nn = squeeze(sum(y<=-thres,1)/N)';
% stats.Nu = squeeze(sum(abs(y)<thres,1)/N)';

