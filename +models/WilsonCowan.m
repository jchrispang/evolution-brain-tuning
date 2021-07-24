function sol = WilsonCowan(param, P_E, P_I, y0)
% WilsonCowan.m
%
% Solves the Wilson-Cowan model without delay using the Euler-Maruyama method
%
% Inputs: param : parameters of the model
%         y0    : initial conditions [Nx1]
%         P_E   : excitatory drive per node [Nx1]
%         P_I   : inhibitory drive per node [Nx1]
%
% Output: sol   : solutions [NxT]
%
% Equation:
% \dot{S^E_i} = (1/tau_E)*{-S^E_i + [1 - S^E_i]*H(x^E_i) + sigma_E*xi_i}
% \dot{S^I_i} = (1/tau_I)*{-S^I_i + [1 - S^I_i]*H(x^I_i) + sigma_I*xi_i}
%    H(x^E_i) = 1/[1 + exp(-a_E*(x^E_i - mu_E))]
%    H(x^I_i) = 1/[1 + exp(-a_I*(x^E_i - mu_I))]
%       x^E_i = w_EE*S^E_i - w_EI*S^I_i + G*sum_j(A_ij*S^E_j) + G^E*P^E_i
%       x^I_i = w_IE*S^E_i - w_II*S^I_i + G^I*P^I_i
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%

sol.x = param.T;
sol.y_E = zeros(param.N, length(param.T));
sol.y_I = zeros(param.N, length(param.T));

if nargin<4
    y_E = 0.001*ones(param.N,1);
    y_I = 0.001*ones(param.N,1);
    y0 = [y_E, y_I];
end
if nargin<3
    P_I = ones(param.N,1);
end
if nargin<2
    P_E = ones(param.N,1);
end

% wiener process
w_coef_E = param.sigma_E; 
w_coef_I = param.sigma_I; 
dW_E = sqrt(param.tstep)*randn(param.N,length(param.T));  
dW_I = sqrt(param.tstep)*randn(param.N,length(param.T)); 

% setting initial condition
y = y0;
sol.y_E(:,1) = y(:,1);
sol.y_I(:,1) = y(:,2);

for k = 2:length(param.T)
    dy = ode_fun(param, y, P_E, P_I);
    y = y + dy*param.tstep + [w_coef_E*dW_E(:,k), w_coef_I*dW_I(:,k)];
    sol.y_E(:,k) = y(:,1);
    sol.y_I(:,k) = y(:,2);
end

end

% ode part of the equations
function dy = ode_fun(param, y, P_E, P_I)
    y_E = y(:,1);
    y_I = y(:,2);
    
    x_E = param.w_EE*y_E - param.w_EI*y_I + param.G*param.A*y_E + param.G_E*P_E;
    x_I = param.w_IE*y_E - param.w_II*y_I + param.G_I*P_I;

    H_E = activation_function(x_E, param.a_E, param.mu_E);
    H_I = activation_function(x_I, param.a_I, param.mu_I);
    
    dy(:,1) = (1/param.tau_E)*(-y_E + (1 - y_E).*H_E);
    dy(:,2) = (1/param.tau_I)*(-y_I + (1 - y_I).*H_I);
end

function H = activation_function(x, a, mu)

    H = 1./(1 + exp(-a*(x - mu)));
end

    
    