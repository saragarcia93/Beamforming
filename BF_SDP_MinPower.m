clc; clear all; close all;

%% CONFIGURATION
% System parameters
M = 2;                       % Number of transmitting antennas
N = 2;                       % Number of receiving antennas
sinrmin = [0.5; 0.7];        % minimum SNR for each RX (according to the bitrate required)
sigma2 = ones(N,1);          % Noise variance
iterations = 1;              % Monte Carlo iterations

% Output parameters
w_ex = zeros(M,iterations);
SINR_ex = zeros(M,iterations);
w_sdp = zeros(M,iterations);
SINR_sdp = zeros(M,iterations);
SINR_sum_ex = 0;
SINR_sum_sdp = 0;

% Begin Monte Carlo simulation
for it = 1:iterations
    % Generate inputs
    Hrici = abs(random('Rician',0,1,N,M)+1i*random('Rician',0,1,N,M));   % Channel matrix, Rician distrib. LoS
    % Call exhaustive solver
    [w_ex(:,it), SINR_ex(:,it)] = exhaustive_solver(Hrici,M);
    % Call SDP solver
    [w_sdp(:,it), SINR_sdp(:,it)] = sdp_solver(Hrici,M);
end

% Parse results
SINR_av_ex = sum(SINR_ex,2) / iterations;  % Average SINR Exhaustive
SINR_av_sdp = sum(SINR_sdp,2) / iterations;  % Average SINR Exhaustive
w_av_ex = sum(w_ex,2) / iterations;  % Average SINR Exhaustive
w_av_sdp = sum(w_sdp,2) / iterations;  % Average SINR Exhaustive
ave_tx_power_sdp = sol1AVE + sol2AVE;  % Average transmitted power 



%% Power allocation. EXHAUSTIVE
function [w, SINR] = exhaustive_solver(Hrici,M)
    % To-do
end

%% Power allocation. SEMI-DEFINITE PROGRAMMING
function [w, SINR] = sdp_solver(Hrici, M)
    % Configure parameters
    SINR = zeros(1,size(Hrici,2));
    w = zeros(1,size(Hrici,2));

    R21 = Hrici(1,:)'*Hrici(1,:);  % Autocorrelation matrices of user 1 channel link
    R22 = Hrici(2,:)'*Hrici(2,:);  % Autocorrelation matrices of user 2 channel link

    X21 = zeros(M,M);              % Inizialize optimization variable for problem (2) 
    X22 = zeros(M,M);              % Inizialize optimization variable for problem (2) 
    
    cvx_begin quiet

        variables X21(m,m) X22(m,m)

        minimize(trace(X21)+trace(X22)) % Maximize total SNR quad_form(X,I)

        subject to
        diag(X21)+diag(X22)<=ones(M,1); % Total transmitted power equals to 1 (Change 1 by Pt or battery level as wished).
        trace(R21*X21)-(sinrmin(1)*trace(R21*X22))>= sigma2(1)*sinrmin(1);  %#ok
        trace(R22*X22)-(sinrmin(2)*trace(R22*X21))>= sigma2(2)*sinrmin(2);  %#ok
        X21 == hermitian_semidefinite(M);  %#ok
        X22 == hermitian_semidefinite(M);  %#ok
        % 0 <= X21 % Needed?
        % 0 <= X22 % Needed?

    cvx_end
   
    % Checking whether the solution is rank 1
    eval1=diag(svd(X21));
    eval2=diag(svd(X22));
    if((sum((eval1(2:end)<10e-6))-sum((eval1(2:end)<10e-6)))==0) % Checking that rank(X21) && rank (X22) are both 1.
    w(1) = sqrt(diag(X21)/sqrt(eval1(1)));                    % Weighs for user 1, for all transmitting antennas
    w(2) = sqrt(diag(X22)/sqrt(eval2(1)));                    % Weighs for user 2, for all transmitting antennas
    else 
      disp('Randomization required')                           % For not robust formulation rank = 1 [1].
    end

    % Monte Carlo
    SINR(1) = trace(R21*X21)/(trace(R21*X22)+sigma2(1));
    SINR(2) = trace(R22*X22)/(trace(R22*X21)+sigma2(2)); 
end
  
%% REF [1] Mats Bengtsson, Björn Ottersten, Optimum and Suboptimum Transmit Beamforming