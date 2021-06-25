function stats = calc_response_stats(x, F, F_transition)
% calc_response_stats.m
%
% Calculate various statistics of F
%
% Inputs: x             : x vector (vector)
%         F             : function value at x [N by length(x)]
%         F_transition  : transition value of F (float)
%
% Output: stats         : statistics (struct)
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, University, 2021

%%

if nargin<3
    F_transition = 0.3;
end

if size(F,2)~=length(x)
    error('Error. Length of x must be the same as the number of columns of F.')
end

for j=1:size(F,1)
    stats.max(j,1) = max(F(j,:));
    stats.min(j,1) = min(F(j,:));
    stats.Fval_10(j,1) = stats.min(j) + 0.1*(stats.max(j)-stats.min(j));
    stats.Fval_90(j,1) = stats.min(j) + 0.9*(stats.max(j)-stats.min(j));
    if isnan(stats.Fval_10(j,1)) || isnan(stats.Fval_90(j,1))
        stats.xval_10(j,1) = nan;
        stats.xval_90(j,1) = nan;
        stats.dynamic_range(j,1) = nan;
        stats.xval_transition(j,1) = nan;
    else
        stats.xval_10(j,1) = x(dsearchn(F(j,:)', stats.Fval_10(j)));
        stats.xval_90(j,1) = x(dsearchn(F(j,:)', stats.Fval_90(j)));
        stats.dynamic_range(j,1) = 10*log10(stats.xval_90(j)/stats.xval_10(j));
        stats.xval_transition(j,1) = x(dsearchn(F(j,:)', F_transition));
    end
end

stats.dx_transition = max(stats.xval_transition) - min(stats.xval_transition);

