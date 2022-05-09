function [Iinj, heterogeneity_text] = calc_Iinj(W, heterogeneity_type, varargin)
% calc_Iinj.m
%
% Calculate vector of Iinj based on heterogeneity_type
%
% Inputs: W                   : connectivity matrix [NxN]
%         heterogeneity_type  : type of heterogeneity (string)
%                               'homogeneous', 'wang', 'hierarchical2',
%                               'hierarchical1', 'degree'
%
% Outputs: Iinj               : input current per node [Nx1]
%          heterogeneity_text : name of heterogeneity (string)
%
% Original: James Pang, Monash University, 2022

%%

temp = load('data/Wang2018_MFMem/Estimated_Parameter_Lausanne.mat', 'external_input');
Wang = struct();
Wang.I = temp.external_input;

normalization = 'maximum';
N = size(W, 1);
A = norm_matrix(W, normalization);
s = sum(A, 2);


switch heterogeneity_type
    case 'homogeneous'
        I0 = varargin{1};
        Iinj = I0*ones(N, 1);
        
        heterogeneity_text = sprintf('homogeneous_%.2f', I0);
        
    case 'wang'
        Iinj = Wang.I;
        
        heterogeneity_text  = 'wang';
        
    case 'hierarchical2'
        Imin = varargin{1};
        Imax = varargin{2};
        Iinj = Imax - (Imax - Imin)*((s - min(s))/(max(s) - min(s))).^2;

        heterogeneity_text = sprintf('hierarchical2_%.2f-%.2f', Imin, Imax);
        
    case 'hierarchical1'
        Imin = varargin{1};
        Imax = varargin{2};
        Iinj = Imax - (Imax - Imin)*((s - min(s))/(max(s) - min(s))).^1;

        heterogeneity_text = sprintf('hierarchical1_%.2f-%.2f', Imin, Imax);
        
    case 'degree'
        Imin = varargin{1};
        Imax = varargin{2};
        [~,s_sort_ind] = sort(s, 'descend');
        Iinj = zeros(N, 1);
        Iinj(s_sort_ind) = linspace(Imin, Imax, N);

        heterogeneity_text = sprintf('degree_%.2f-%.2f', Imin, Imax);    
end
