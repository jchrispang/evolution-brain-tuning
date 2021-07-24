function [acf, lags, tau] = calc_timescales(F, maxlag, dt, tau_0)
% calc_timescales.m
%
% Calculate timescales of F
%
% Inputs: F      : function vale [NxT]
%         maxlag : maximum time lag (float)
%         dt     : time step (float)
%
% Outputs: acf   : autocorrelation function [Nxlags]
%          lags  : time lags [vector]
%          tau   : timescales [vector]
%
% Original: James Pang, QIMR Berghofer, 2021
% Revised:  James Pang, Monash University, 2021

%%

if nargin<4
    tau_0 = 0.01;
end

N = size(F,1);
tau = zeros(size(F,1), 1);

for node = 1:N
%     [acf(node,:), temp_lags, ~] = autocorr(F(node,:), 'numlags', floor(maxlag/dt));
%     lags = temp_lags*dt;
    [acf_temp, temp_lags] = xcorr(F(node,:)-mean(F(node,:)), floor(maxlag/dt));
    acf_temp = acf_temp((floor(length(temp_lags)/2)+1):end);
    acf(node,:) = acf_temp/acf_temp(1);
    lags = temp_lags((floor(length(temp_lags)/2)+1):end)*dt;
    
    fun = @(beta,x) beta(1)*exp(-x/beta(2))+beta(3);
    beta0 = [1, tau_0, 0];
    
    % mdl = fitnlm(lags, acf(node,:), fun, beta0);
%     beta = nlinfit(lags, acf(node,:), fun, beta0);
    
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    beta = lsqcurvefit(fun, beta0, lags, acf(node,:), [], [], options);  
    
    tau(node) = beta(2);
end