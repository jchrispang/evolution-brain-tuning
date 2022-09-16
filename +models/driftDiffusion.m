function sol = driftDiffusion(param, lambda, y0)
% driftDiffusion.m
%
% Solves the dynamic MFM model with one excitatory/inhibitory population
%
% Inputs: param  : parameters of the model
%         lambda : self-coupling (float)
%         y0     : initial conditions [Nx1]
%
% Output: sol   : output (struct)
%                 fields: x, y
%                 x = time [1xT]
%                 y = solutions [NxT]
%
% Equation:
% \dot{y} = beta_i + lambda*y - Ly + sigma*v_i(t)
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2022

%%

if nargin<3
    y0 = zeros(param.N,1);
end
if nargin<2
    lambda = 0;
end

sol.x = param.T;
sol.y = zeros(param.N, length(param.T));

% calculate laplacian matrix
W = param.A;
W_A = double(W>0);
W_D = diag(sum(W_A,2));
W_L = W_D - W_A;

% wiener process
w_coef = param.sigma; 
dW = sqrt(param.tstep)*randn(param.N,length(param.T)); 

% setting initial condition
sol.y(:,1) = y0;

for k = 2:length(param.T)
    dy = param.beta.*ones(param.N,1) + lambda*sol.y(:,k-1) - W_L*sol.y(:,k-1);
    sol.y(:,k) = sol.y(:,k-1) + dy*param.tstep + w_coef*dW(:,k);
    
    if sum(sol.y(:,k-1)>=param.thres) > 0
        sol.y(sol.y(:,k-1)>=param.thres, k-1) = param.thres;
        sol.y(sol.y(:,k-1)>=param.thres, k) = param.thres;
    end
    if sum(sol.y(:,k-1)<=-param.thres) > 0
        sol.y(sol.y(:,k-1)<=-param.thres, k-1) = -param.thres;
        sol.y(sol.y(:,k-1)<=-param.thres, k) = -param.thres;
    end
end