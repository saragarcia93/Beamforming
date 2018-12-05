clc; clear all; close all;        

%% CONFIGURATION
% System parameters
antenna_allocation=[1,0;0,1;0,1]; % Antenna allocation, input parameter, dim(MxN)
M = size(antenna_allocation,1);                            % Number of transmitting antennas
N = size(antenna_allocation,2);                            % Number of users
sinrmin = [0.1; 0.3];             % Minimum SINR for each user according to the bitrate required
sigma2 = ones(N,1);               % Noise variance
Pt_max = 1;                       % Maximum transmitted power per each radius

% Output parameters
w_sdp = zeros(M,N);    % Weighs provided by the SDP solver
SINR_sdp = zeros(N);   % SINR provided by the SDP solver

% Generate inputs
%H = [channelEst1;channelEst2];
H = (random('Rician',0,1,N,M)+1i*random('Rician',0,1,N,M)); 
H_abs = abs(H);
Hrici_angle=angle(H);

% Call SDP solver
[w_sdp, SINR_sdp] = sdp_solver(H_abs,M,N,Pt_max,sinrmin,sigma2,antenna_allocation);
ang_w = -1.*Hrici_angle;
weighs = w_sdp.*cos(ang_w')+1i.* w_sdp.*sin(ang_w');  % Each column is the transmitted weighs for each user

%% Power allocation. SEMI-DEFINITE PROGRAMMING
function [w, SINR] = sdp_solver(H_abs, M,N, Pt_max, sinrmin, sigma2,antenna_allocation)
    % Configure parameters
    SINR = zeros(1,size(H_abs,1)); 
    w = zeros(size(H_abs,2),size(H_abs,1));

    R_user1 = H_abs(1,:)'*H_abs(1,:);  % Autocorrelation matrices of user 1 channel link
    R_user2 = H_abs(2,:)'*H_abs(2,:);  % Autocorrelation matrices of user 2 channel link

    X_user1 = zeros(M,M);              % Inizialize optimization variable for problem (2) 
    X_user2 = zeros(M,M);              % Inizialize optimization variable for problem (2) 

    % Creating association matrices
    association_A = zeros(M,M,N);
    association_B = ones(M,M,N);
    for i=1:M
        for j=1:N
            if(antenna_allocation(i,j)==1)
              association_A(i,i,j)=1;
            elseif (antenna_allocation(i,j)==0)
              association_B(:,i,j)=0;   
              association_B(i,:,j)=0; 
            end
        end
    end
    
    cvx_begin %quiet

        variables X_user1(M,M) X_user2(M,M)

        minimize(trace(association_A(:,:,1).*X_user1)+trace(association_A(:,:,2).*X_user2)) % Maximize total SINR 

        subject to
        diag(association_A(:,:,1).*X_user1)+diag(association_A(:,:,2).*X_user2)<=Pt_max*ones(M,1);                                      %#ok. Constraint on total transmitted power per antenna 
        trace(R_user1*(association_B(:,:,1).*X_user1))-(sinrmin(1)*trace(R_user1*(association_B(:,:,2).*X_user2)))>= sigma2(1)*sinrmin(1);  %#ok. Contraint on minimum SINR for user 1
        trace(R_user2*(association_B(:,:,2).*X_user2))-(sinrmin(2)*trace(R_user2*(association_B(:,:,1).*X_user1)))>= sigma2(2)*sinrmin(2);  %#ok. Contraint on minimum SINR for user 2
        X_user1 == hermitian_semidefinite(M) ;                                              %#ok. X_user1=X_user1*, X_user1>=0
        X_user2 == hermitian_semidefinite(M);                                               %#ok. X_user2=X_user2*, X_user2>=0
        
    cvx_end
    
   X_user1 = X_user1.*association_B(:,:,1);
   X_user2 = X_user2.*association_B(:,:,2);
   
    % Checking whether the solution is rank 1
    eval1=diag(svd(X_user1));
    eval2=diag(svd(X_user2));
    if((sum((eval1(2:end)<10e-6))-sum((eval2(2:end)<10e-6)))==0) % Checking that rank(X21) && rank (X22) are both 1.
      sol_1 = sqrt(diag(X_user1));                                   % Weighs for user 1, for all transmitting antennas
      sol_2 = sqrt(diag(X_user2));                                   % Weighs for user 2, for all transmitting antennas
      w =[sol_1,sol_2];
      sinr_1test = ((H_abs(1,:)*sol_1)^2)/(((H_abs(1,:)*sol_2)^2)+sigma2(1));           % SINR achieved for user 1
      sinr_2test = ((H_abs(2,:)*sol_2)^2)/(((H_abs(2,:)*sol_1)^2)+sigma2(2));           % SINR achieved for user 2
      SINR = [sinr_1test,sinr_2test];
     else 
       disp('Randomization required')                             % For not robust formulation rank = 1 [1].
     end 
end

%% REF [1] Mats Bengtsson, Björn Ottersten, Optimum and Suboptimum Transmit Beamforming
