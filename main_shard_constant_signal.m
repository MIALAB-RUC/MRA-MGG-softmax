%% 文章名称，作者
%代码的网址
clc;clear all;close all;
addpath('../utilities')

%% Check that Manopt is available
if ~exist('manopt_version', 'file')
    error(sprintf(['Please get Manopt 3.0 or later,\n' ...
                   'available from http://www.manopt.org.\n' ...
                   'Then, run importmanopt.m to add it to the path.'])); %#ok<SPERR>
end
fprintf('This could take a few minutes.\n');
N = 101;
x_true = zeros(N, 1);
x_true(31:60) = 1;

%% Start the parallel pool
parallel_nodes = 2;
if isempty(gcp('nocreate'))
    parpool(parallel_nodes, 'IdleTimeout', 240);
end
M = 25000;
alpha_true = 0.2;
sigma1_true = 10;
sigma2_true = 0.1;

%% Generate M observations with noise variance sigma^2
[X_data,T,shifts] = generate_observations_mix(x_true, M,alpha_true, sigma1_true,sigma2_true);

%% Estimate x from the data through an adaptive variational model
alpha=0.5;
sigma1=50;
sigma2=0.15;
lam=5;
r=20*lam;
tic_EM_proposed = tic();
u=MRA_MGG_softmax(X_data,alpha,sigma1,sigma2,lam,r);
time_proposed=toc(tic_EM_proposed);
relerr_proposed = relative_error(x_true, u); % 2-norm, up to integer shifts
x_proposed= align_to_reference(u, x_true);

%% It is informative to compare against EM
sigma=4.8;
tic_original_em = tic();
x_original_em = MRA_EM_priori_bispectrum(X_data, sigma);
time_em = toc(tic_original_em);
relerr_em = relative_error(x_true, x_original_em); % 2-norm, up to integer shifts
x_original_em = align_to_reference(x_original_em, x_true);

%%  It is informative to compare against invariant 
% Estimate the invariants once (this is common to many methods)
tic_invariants = tic();
Emean = mean(mean(X_data));
mean_est = mean(X_data);
mean_est = repmat(mean_est,N,1);
X_data = X_data - mean_est;
Y_hat = fft(X_data,[],1);   
Y_power = mean(abs(Y_hat).^2,2)-N*sigma^2;
Y_power = max(0, Y_power);
meanB = get_bispectrum(Y_hat);
B_phase = get_B_phase(meanB);
time_invariants=toc(tic_invariants);

%% gap
tic_bispectrum_inversion_1= tic();
Est_phase_1 = get_phase_from_bispectrum_gap(B_phase,N); 
x_est_1  = get_recon(Emean,Y_power,Est_phase_1);
time_bispectrum_inversion_1=toc(tic_bispectrum_inversion_1);
time_invfeatmra_1=time_invariants+time_bispectrum_inversion_1;
relerr_invfeatmra_1 = relative_error(x_true, x_est_1); % 2-norm, up to integer shifts
x_est_1 = align_to_reference(real(x_est_1), x_true);

%% real
tic_bispectrum_inversion_2= tic();
[Est_phase_2, problem] = phases_from_bispectrum_real(meanB,1); 
x_est_2 = get_recon(Emean,Y_power,Est_phase_2);
time_bispectrum_inversion_2=toc(tic_bispectrum_inversion_2);
time_invfeatmra_2=time_invariants+time_bispectrum_inversion_2;
x_est_2=real(x_est_2);
relerr_invfeatmra_2 = relative_error(x_true, x_est_2); % 2-norm, up to integer shifts
x_est_2 = align_to_reference(real(x_est_2), x_true);

%% FM_real
y = fft(x_true);
tic_bispectrum_inversion_3= tic();
z_est_3 = phases_from_bispectrum_FM_real(meanB,sign(y(1)), sign(y(2))); 
x_est_3  = get_recon(Emean,Y_power,z_est_3);
time_bispectrum_inversion_3=toc(tic_bispectrum_inversion_3);
time_invfeatmra_3=time_invariants+time_bispectrum_inversion_3;
relerr_invfeatmra_3 = relative_error(x_true, x_est_3); % 2-norm, up to integer shifts
x_est_3 = align_to_reference(real(x_est_3), x_true);

%% It is informative to compare against invariant  APS_real
tic_bispectrum_inversion_4= tic();
z_est_4=phases_from_bispectrum_APS_real(meanB); 
x_est_4  = get_recon(Emean,Y_power,Est_phase_2);
time_bispectrum_inversion_4=toc(tic_bispectrum_inversion_4);
time_invfeatmra_4=time_invariants+time_bispectrum_inversion_4;
x_est_4=real(x_est_4);
relerr_invfeatmra_4 = relative_error(x_true, x_est_4);
x_est_4 = align_to_reference(x_est_4, x_true);

%% relative error and time
fprintf('Proposed model. RMSE: %.4g; time: %.4g [s].\n',  relerr_proposed,time_proposed );
fprintf('Expectation maximization. RMSE: %.4g; time: %.4g [s]\n', relerr_em, time_em);
fprintf('Spectral M. largest spectral gap. RMSE: %.4g; time: %.4g [s].\n', relerr_invfeatmra_1, time_invfeatmra_1);
fprintf('Invariant features, optim. phase manifold. RMSE: %.4g; time: %.4g [s]\n', relerr_invfeatmra_2, time_invfeatmra_2);
fprintf('Invariant features, FM. RMSE: %.4g; time: %.4g [s]\n', relerr_invfeatmra_3, time_invfeatmra_3);
fprintf('Invariant features,iter, phase synch. RMSE: %.4g; time: %.4g [s]\n', relerr_invfeatmra_4, time_invfeatmra_4);
 
%% plot
figure(1),
hold all;
t = 0:(N-1);
plot(t, x_true, 'k-');
plot(t,x_proposed ,'r.-');
plot(t, x_original_em, 'g.-');
plot(t, x_est_1, 'b.-','linewidth',0.1);
plot(t, x_est_2, 'c--o');
plot(t, x_est_3, 'm-','linewidth',0.4);
plot(t, x_est_4, 'r-*','linewidth',0.7);
title(sprintf('Multireference alignment example, M = %d, sigma = %g', M, sigma));
h=legend('True signal', ...
    sprintf('Proposed model (RMSE: %.4g; time: %.4g [s])', relerr_proposed, time_proposed),...
    sprintf('Expectation maximization (RMSE: %.4g; time: %.4g [s])', relerr_em, time_em),...
    sprintf('Spectral M. largest spectral gap (RMSE: %.4g; time: %.4g [s])', relerr_invfeatmra_1, time_invfeatmra_1), ...
    sprintf('Invariant features, optim. phase manifold (RMSE: %.4g; time: %.4g [s])', relerr_invfeatmra_2, time_invfeatmra_2), ...
    sprintf('Invariant features, FM (RMSE: %.4g; time: %.4g [s])', relerr_invfeatmra_3, time_invfeatmra_3), ...
    sprintf('Invariant features,iter, phase synch(RMSE: %.4g; time: %.4g [s])', relerr_invfeatmra_4, time_invfeatmra_4));
   set(h,'fontsize',7);
xlim([0, N-1]);
ylim([0, 1.5]);
