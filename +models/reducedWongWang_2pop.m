function sol = reducedWongWang_2pop_fast(param, y0)
% reducedWongWang_2pop_fast.m
%
% Solves the reduced Wong-Wang model with one excitatory/inhibitory population
% using the Euler-Maruyama method
%
% Inputs: param : parameters of the model
%         y0    : initial conditions [Nx1]
%
% Output: sol   : solutions [NxT]
%
% Equation:
% \dot{S^E_i} = -S^E_i/tau_E + gamma*(1-S^E_i)*H(x^E_i) + sigma*v_i(t)
% \dot{S^I_i} = -S^I_i/tau_I + H(x^I_i) + sigma*v_i(t)
%    H(x^E_i) = (a_E*x^E_i - b_E)/[1 - exp(-d_E*(a_E*x^E_i - b_E))]
%    H(x^I_i) = (a_I*x^I_i - b_I)/[1 - exp(-d_I*(a_I*x^I_i - b_I))]
%       x^E_i = w_E*I0 + w_+*J_E*S^E_i + G*J_E*sum_j(A_ij*S^E_j) - J*S^I_i + Iext
%       x^I_i = w_I*I0 + J_E*S^E_i + lambda*G*J_E*sum_j(A_ij*S^E_j) - S^I_i
%
% Original: James Pang, QIMR Berghofer, 2020
% Revised:  James Pang, Monash University, 2021

%%

sol.x = param.T;
sol.y_E = zeros(param.N, length(param.T));
sol.y_I = zeros(param.N, length(param.T));

if nargin<2
    y_E = 0.001*ones(param.N,1);
    y_I = 0.001*ones(param.N,1);
    y0 = [y_E, y_I];
end

%wiener process
w_coef = param.sigma; 
dW_E = sqrt(param.tstep)*randn(param.N,length(param.T));  
dW_I = sqrt(param.tstep)*randn(param.N,length(param.T)); 

% setting initial condition
y = y0;
sol.y_E(:,1) = y(:,1);
sol.y_I(:,1) = y(:,2);

for k = 2:length(param.T)
    dy = ode_fun(param, y);
    y = y + dy*param.tstep + w_coef*[dW_E(:,k), dW_I(:,k)];
    sol.y_E(:,k) = y(:,1);
    sol.y_I(:,k) = y(:,2);
end

end

% ode part of the equations
function dy = ode_fun(param, y)
    y_E = y(:,1);
    y_I = y(:,2);
    
    x_E = param.w_E*param.I0 + param.w_p*param.J_E.*y_E + ...
          param.G*param.J_E*param.A*y_E - param.J_I*y_I + param.Iext;
    x_I = param.w_I*param.I0 + param.J_E.*y_E + ...
          param.lambda*param.G*param.J_E*param.A*y_E - y_I;  
    
    H_E = firing_rate(x_E, param.a_E, param.b_E, param.d_E);
    H_I = firing_rate(x_I, param.a_I, param.b_I, param.d_I);
    
    dy(:,1) = -y_E/param.tau_E + param.gamma*(1 - y_E).*H_E;
    dy(:,2) = -y_I/param.tau_I + H_I;
end

function H = firing_rate(x, a, b, d)

    H = (a*x - b)./(1 - exp(-d*(a*x - b)));
end
