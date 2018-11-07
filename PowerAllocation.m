%% POWER ALLOCATION FOR CENTRALIZED BEAMFORMING ACCORDING TO SINR CONSTRAINTS
 % (1) Maximizes the total SINR at the receivers
 % (2) Minimizes the sum of the power transmitted 
 % (3) Minimizes the difference between the SINR achieved and the SINR targeted 
 
clc
clear all;
close all;

% System parameters
   count=0;
   m=3;                         % Number of transmitting antennas
   n=2;                         % Number of receiving antennas
   sinrmin = [0.5; 0.7];        % minimum SNR for each RX (according to the bitrate required)
   
% Monte Carlo variables
   iterations=5;
   txpowerEach1=zeros(m,1);     % For approach (1)
   rxSINR1 = zeros(n,1);
   txpowerEach2=zeros(m,1);     % For approach (2)
   rxSINR2 = zeros(n,1);
   txpowerEach3=zeros(m,1);     % For approach (3)
   rxSINR3 = zeros(n,1);
   
% Power allocation
   
for ii=1:iterations
   Hrici = abs(random('Rician',0,1,n,m)+1i*random('Rician',0,1,n,m));   % Channel matrix, Rician distrib. LoS
   
   X1 = zeros(m,n);          % Inizialize optimization variable for problem (1)
   X2 = zeros(m,n);          % Inizialize optimization variable for problem (2)
   X3 = zeros(m,n);          % Inizialize optimization variable for problem (3)
   
  % (1) Maximizes the total SINR at the receivers 
 cvx_begin quiet
    variable X1(m,n) 
    maximize( trace(Hrici*X1) )  % Maximize total SNR
    subject to
    0 <= X1 <= 1;            % Box relaxation
    sum (X1') <= ones(1,m);  % Each tx antenna for only one rx
    diag(Hrici*X1)-(sinrmin.*(sum(((Hrici*X1)-diag(diag(Hrici*X1)))'))')>= sinrmin ;
 cvx_end
 
  % Monte Carlo 
  txpowerEach1=txpowerEach1+diag(X1*X1');
  rxSINR1=rxSINR1+ diag(Hrici*X1)-(sinrmin.*(sum(((Hrici*X1)-diag(diag(Hrici*X1)))'))');
  
  
  % (2) Minimizes the sum of the power transmitted 
 cvx_begin quiet
   variable X2(m,n) 
   minimize(sum(sum(X2.^2)))  % Maximize total SNR quad_form(X,I)
   subject to
   0 <= X2 <= 1;            % Box relaxation
   sum (X2') <= ones(1,m);  % Each tx antenna for only one rx
   diag(Hrici*X2)-(sinrmin.*(sum(((Hrici*X2)-diag(diag(Hrici*X2)))'))')>= sinrmin ;
 cvx_end
  
  % Monte Carlo 
  txpowerEach2=txpowerEach2+diag(X2*X2');
  rxSINR2=rxSINR2+ diag(Hrici*X2)-(sinrmin.*(sum(((Hrici*X2)-diag(diag(Hrici*X2)))'))');
  
  
  % (3) Minimizes the difference between the SINR achieved and the SINR targeted 
 cvx_begin quiet    
    variable X3(m,n) 
    minimize(sum((diag(Hrici*X3))-(sinrmin.*(sum(((Hrici*X3)-diag(diag(Hrici*X3)))'))'))) 
    subject to
    0 <= X3 <= 1;            % (C1) Box relaxation (in this case allow to assign one tx antenna to more than one rx)
    sum (X3') <= ones(1,m);  % (C2) Normalizes transmitted power
   (diag(Hrici*X3))-(sinrmin.*(sum(((Hrici*X3)-diag(diag(Hrici*X3)))'))') >= sinrmin ; % (C3) SINR constraint for each user
 cvx_end

% Monte Carlo 
  txpowerEach3=txpowerEach3+diag(X3*X3');
  rxSINR3=rxSINR3+ diag(Hrici*X3)-(sinrmin.*(sum(((Hrici*X3)-diag(diag(Hrici*X3)))'))');
    
end

% Monte Carlo
AVtxpowerEach1=txpowerEach1./iterations;   % For approach (1)
AVErxSINR1=rxSINR1./iterations;
totalAVPower1=sum(AVtxpowerEach1);
totalAVSINR1=sum(AVErxSINR1);

AVtxpowerEach2=txpowerEach2./iterations;   % For approach (2)
AVErxSINR2=rxSINR2./iterations;
totalAVPower2=sum(AVtxpowerEach2);
totalAVSINR2=sum(AVErxSINR2);

AVtxpowerEach3=txpowerEach3./iterations;   % For approach (3)
AVErxSINR3=rxSINR3./iterations;
totalAVPower3=sum(AVtxpowerEach3);
totalAVSINR3=sum(AVErxSINR3);

%% TESTS, for each iteration
% (1) 
% diag(Hrici*X1)-(snrmin.*(sum(((Hrici*X1)-diag(diag(Hrici*X1)))'))')>= snrmin  % The solution given satisfies all SINR constraints
% X1                                    % Computed power allocation
% powerCVX1 = sum(diag(X1*X1'))          % Transmitted power, computed solution
% sinrmin                               % Minimim SINR
% diag(Hrici*X1)-(sinrmin.*(sum(((Hrici*X1)-diag(diag(Hrici*X1)))'))')   % SINR, computed solution

% (2) 
% diag(Hrici*X2)-(snrmin.*(sum(((Hrici*X2)-diag(diag(Hrici*X2)))'))')>= snrmin  % The solution given satisfies all SINR constraints
% X2                                    % Computed power allocation
% powerCVX2 = sum(diag(X2*X2'))          % Transmitted power, computed solution
% sinrmin                               % Minimim SINR
% diag(Hrici*X2)-(sinrmin.*(sum(((Hrici*X2)-diag(diag(Hrici*X2)))'))')   % SINR, computed solution

% (3) 
% diag(Hrici*X1)-(snrmin.*(sum(((Hrici*X1)-diag(diag(Hrici*X1)))'))')>= snrmin  % The solution given satisfies all SINR constraints
% X1                                    % Computed power allocation
% powerCVX3 = sum(diag(X3*X3'))          % Transmitted power, computed solution
% sinrmin                               % Minimim SINR
% diag(Hrici*X3)-(sinrmin.*(sum(((Hrici*X3)-diag(diag(Hrici*X3)))'))')   % SINR, computed solution
