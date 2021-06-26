function W_norm = norm_matrix(W, normalization)
% norm_matrix.m
%
% Normalize connectome according to method
%
% Inputs: W              : connectivity matrix [NxN]
%         normalization  : normalization method to use [string]
%                          'maximum', 'in_strength', 'out_strength'
%
% Output: W_norm         : normalized matrix [NxN]
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%

if nargin<2
    normalization = 'maximum';
end

if strcmpi(normalization, 'maximum')
    W_norm = W/max(W,[],'all');
elseif strcmpi(normalization, 'in_strength')
    node_strength = sum(W,2);
    W_norm = bsxfun(@rdivide, W, node_strength);
elseif strcmpi(normalization, 'out_strength')
    node_strength = sum(W,1);
    W_norm = bsxfun(@rdivide, W, node_strength);
end

W_norm(isnan(W_norm)) = 0;
