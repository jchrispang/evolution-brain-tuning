% demo_simulation.m

%% load connectome data
load('data/connectome_human.mat')
load('data/connectome_chimp.mat')
load('data/connectome_macaque.mat')
load('data/connectome_marmoset.mat')

%% generating and analyzing time series from the reduced Wong-Wang model

% load predefined reduced Wong-Wang model parameters
param = utils.loadParameters_reducedWongWang_func;

% define connectome matrix A
type = 'human';
param.A = eval(sprintf('connectome_%s', type));  % replace with your own connectivity matrix
param.N = size(param.A, 2);

% normalize connectivity matrix with respect to maximum weight
normalization = 'maximum';
param.A = utils.norm_matrix(param.A, normalization);  

% define simulation time 
tpre =  100;                                     % burn time to remove transient
tpost = 50;                                      % steady-state time (it is advisable to increase this to reach ergodicity)
param.tmax = tpre + param.tstep + tpost;         
param.tspan = [0, param.tmax];
param.T = 0:param.tstep:param.tmax;

% simulate model using predefined initial conditions (y0 = 0.01)
sol = models.reducedWongWang_fast(param);

% obtain steady-state synaptic gating and firing rate time series per region
time_steady_ind = dsearchn(param.T', tpre)+1;                       % index of start of steady state
S_steady = sol.y(:,time_steady_ind:end);                            % synaptic gating
H_steady = utils.calc_firingRate_reducedWongWang(param, S_steady);  % firing rate

% redefine time vector to remove burn time
T_steady = (0:size(S_steady,2)-1)*param.tstep;

% calculate some stats of S_steady per region
maxlag = 5;
[acf, lags, tau] = utils.calc_timescales(S_steady, maxlag, param.tstep); % timescale
Smean = mean(S_steady,2);                                                % mean synaptic gating

% plot steady-state synaptic gating and firing rate time series
figure;
subplot(1,2,1)
plot(T_steady, S_steady)
xlabel('time')
ylabel('regional synaptic gating, S')

subplot(1,2,2)
plot(T_steady, H_steady)
xlabel('time')
ylabel('regional firing rate, H')

% plot some regional stats
figure;
yyaxis left
plot(1:param.N, Smean, '.-')
xlabel('region')
ylabel('mean synaptic gating')

yyaxis right
plot(1:param.N, tau, '.-')
xlabel('region')
ylabel('timescale, \tau')


%% calculating and analyzing synaptic gating tuning curve from the reduced Wong-Wang model   
   
% load predefined reduced Wong-Wang model parameters
param = utils.loadParameters_reducedWongWang_func;

% define connectome matrix A
type = 'human';
param.A = eval(sprintf('connectome_%s', type));  % replace with your own connectivity matrix
param.N = size(param.A, 2);

% normalize connectivity matrix with respect to maximum weight
normalization = 'maximum';
param.A = utils.norm_matrix(param.A, normalization);  

% define simulation time 
tpre =  100;                                     % burn time to remove transient
tpost = 50;                                      % steady-state time (it is advisable to increase this to reach ergodicity)
param.tmax = tpre + param.tstep + tpost;         
param.tspan = [0, param.tmax];
param.T = 0:param.tstep:param.tmax;

% define vector of recurrent connection strengths and number of trials
w_vec = 0.05:0.05:1;                             % make sure to increase resolution for a proper analysis
num_trials = 1;

% calculate tuning curve per region
[Smean, Hmean, Hmax] = utils.calc_tuning_reducedWongWang(param, w_vec, tpre, num_trials);

% calculate some stats of Smean
S_transition = 0.3;                              % transition value of Smean
stats = utils.calc_response_stats(w_vec, Smean, S_transition);

% plot synaptic gating tuning curve and dynamic range per region
figure;
subplot(1,2,1)
plot(w_vec, Smean)
xlabel('recurrent strength, w')
ylabel('regional synaptic gating, S')

subplot(1,2,2)
plot(1:param.N, stats.dynamic_range, 'k.-')
xlabel('region')
ylabel('dynamic range')
