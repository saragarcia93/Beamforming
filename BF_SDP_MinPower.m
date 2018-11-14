%% POWER ALLOCATION FOR CENTRALIZED BEAMFORMING ACCORDING TO SINR CONSTRAINTS
 % Minimizes the sum of the power transmitted 

clc
clear all;
close all;

% System parameters
   count=0;
   m=3;                         % Number of transmitting antennas
   n=2;                         % Number of receiving antennas
   sinrmin = [0.5; 0.7];        % minimum SNR for each RX (according to the bitrate required)
   sigma2 = ones(n,1);          % Noise variance
   
% Monte Carlo variables
   iterations=1;
   sol1=zeros(m,1);   
   sol2=zeros(m,1); 
   sinr1 = 0;
   sinr2 = 0;

%% Power allocation. SEMI-DEFINITE PROGRAMMING 

for ii=1:iterations
   Hrici = abs(random('Rician',0,1,n,m)+1i*random('Rician',0,1,n,m));   % Channel matrix, Rician distrib. LoS
   R21 = Hrici(1,:)'*Hrici(1,:);  % Autocorrelation matrices of user 1 channel link
   R22 = Hrici(2,:)'*Hrici(2,:);  % Autocorrelation matrices of user 2 channel link
   
   X21 = zeros(m,m);              % Inizialize optimization variable for problem (2) 
   X22 = zeros(m,m);              % Inizialize optimization variable for problem (2) 
    
 cvx_begin quiet
 
   variables X21(m,m) X22(m,m)
   
   minimize(trace(X21)+trace(X22)) % Maximize total SNR quad_form(X,I)
   
   subject to
   diag(X21)+diag(X22)<=ones(m,1); % Total transmitted power equals to 1 (Change 1 by Pt or battery level as wished).
   trace(R21*X21)-(sinrmin(1)*trace(R21*X22))>= sigma2(1)*sinrmin(1)
   trace(R22*X22)-(sinrmin(2)*trace(R22*X21))>= sigma2(2)*sinrmin(2)
   X21 == hermitian_semidefinite(m)
   X22 == hermitian_semidefinite(m)
   % 0 <= X21 % Needed?
   % 0 <= X22 % Needed?
    
 cvx_end
   
  % Checking whether the solution is rank 1
  eval1=diag(svd(X21));
  eval2=diag(svd(X22));
  if((sum((eval1(2:end)<10e-6))-sum((eval1(2:end)<10e-6)))==0) % Checking that rank(X21) && rank (X22) are both 1.
  sol1=sol1+sqrt(diag(X21)/sqrt(eval1(1)));                    % Weighs for user 1, for all transmitting antennas
  sol2=sol2+sqrt(diag(X22)/sqrt(eval2(1)));                    % Weighs for user 2, for all transmitting antennas
  else 
      disp('Randomization required')                           % For not robust formulation rank = 1 [1].
  end
 
  % Monte Carlo 
  sinr1 = sinr1 + trace(R21*X21)/(trace(R21*X22)+sigma2(1));
  sinr2 = sinr2 + trace(R22*X22)/(trace(R22*X21)+sigma2(2));
  
end

%% Monte Carlo
  sinr1AVE=sinr1/iterations;
  sinr2AVE=sinr2/iterations;
  sol1AVE=sol1./iterations;
  sol2AVE=sol2./iterations;
  ave_tx_power_sdp=sol1AVE+sol2AVE;   % Average transmitted power 
  
%% REF [1] Mats Bengtsson, Björn Ottersten, Optimum and Suboptimum Transmit Beamforming