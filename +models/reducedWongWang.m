function sol = reducedWongWang(param, y0)
% reducedWongWang.m
%
% Solves the reduced Wong-Wang model with one excitatory population
% using the Euler-Maruyama method
%
% Inputs: param : parameters of the model
%         y0    : initial conditions [Nx1]
%
% Output: sol   : output (struct)
%                 fields: x, y
%                 x = time [1xT]
%                 y = solutions [NxT]
%
% Equation:
% \dot{S_i} = -S_i/tau_s + gamma_s*(1-S_i)*H(x_i) + sigma*v_i(t)
%    H(x_i) = (a*x_i - b)/[1 - exp(-d*(a*x_i - b))]
%       x_i = w*J*S_i + G*J*sum_j(A_ij*S_j) + I
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%

if nargin<2
    y0 = 0.001*ones(param.N,1);
end

sol.x = param.T;
sol.y = zeros(param.N, length(param.T));

% wiener process
w_coef = param.sigma; 
dW = sqrt(param.tstep)*randn(param.N,length(param.T)); 

% setting initial condition
sol.y(:,1) = y0;

for k = 2:length(param.T)
    dy = ode_fun(param, sol.y(:,k-1));
    sol.y(:,k) = sol.y(:,k-1) + dy*param.tstep + w_coef*dW(:,k);
end

end

% ode part of the equations
function dy = ode_fun(param, y)
    x = param.w*param.J.*y + param.G*param.J*param.A*y + param.I;
    H = (param.a*x - param.b)./(1 - exp(-param.d*(param.a*x - param.b)));
    dy = -y/param.tau_s + param.gamma_s*(1 - y).*H;
end
