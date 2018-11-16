clc; clear all; close all;        %#ok

%% CONFIGURATION
% System parameters
M = 2;                            % Number of transmitting antennas
N = 2;                            % Number of users
sinrmin = [0.5; 0.7];             % Minimum SINR for each user according to the bitrate required
sigma2 = ones(N,1);               % Noise variance
Pt_max = 1;                       % Maximum transmitted power per each radius
iterations = 1;                   % Number of Monte Carlo iterations

% Output parameters
w_ex = zeros(M,N,iterations);     % Weighs provided by the exhaustive solver
SINR_ex = zeros(N,iterations);    % SINR provided by the exhaustive solver
w_sdp = zeros(M,N,iterations);    % Weighs provided by the SDP solver
SINR_sdp = zeros(N,iterations);   % SINR provided by the SDP solver
SINR_sum_ex = 0;
SINR_sum_sdp = 0;

% Begin Monte Carlo simulation
for it = 1:iterations
    % Generate inputs
    Hrici = abs(random('Rician',0,1,N,M)+1i*random('Rician',0,1,N,M));   % Channel matrix, Rician distribution. LoS. 
    % Call exhaustive solver
    %[w_ex(:,:,it), SINR_ex(:,it)] = exhaustive_solver(Hrici,M,N,Pt_max,sinrmin,sigma2);
    % Call SDP solver
    [w_sdp(:,:,it), SINR_sdp(:,it)] = sdp_solver(Hrici,M,Pt_max,sinrmin,sigma2);
end

%% Parse results  CHECKED
 SINR_av_ex = sum(SINR_ex,2) / iterations;    % Average SINR Exhaustive 
 SINR_av_sdp = sum(SINR_sdp,2) / iterations;  % Average SINR Exhaustive 
 w_av_ex = sum(w_ex,3) / iterations;          % Average SINR Exhaustive 
 w_av_sdp = sum(w_sdp,3) / iterations;        % Average SINR Exhaustive 
 ave_tx_power_sdp = sum(w_av_sdp.^2,2);       % Average transmitted power per radius 



%% Power allocation. EXHAUSTIVE
function [w, SINR] = exhaustive_solver(Hrici,M,N,Pt_max,sinrmin,sigma2)
    % Creo una matriz de dimension 3x2 que contiene los pesos de cada antena
% para cada user. Cada fila pesos de un user
exh = zeros(M,N,10);                                         % Where each column corresponds to a different user weighs.
cont=1;                         
step=0.01;                                                      % Echaustive search step
for a=0:step:Pt_max
    for b=0:step:Pt_max
        for c=0:step:Pt_max
            for d=0:step:Pt_max
%                 for e=0:step:Pt_max
%                     for f=0:step:Pt_max
                    exh(:,:,cont)=[a,b;c,d];               % Creating all solution matrices combinations. 2 radius
%                   exh(:,:,cont)=[a,b;c,d;e,f];           % Creating all solution matrices combinations. 3 radius
                    cont=cont+1;
%                     end 
%                 end
            end
        end
    end
end

% Checking maximum power constraint % BIEN
size_exh=size(exh);
exh2 = zeros(M,N,10); 
count2=1;
for i=1:size_exh(3)
    if (sum((sum(((exh(:,:,i)).^2),2)<=Pt_max*ones(M,1)))==M)
        exh2(:,:,count2)=exh(:,:,i);
        count2=count2+1;
    end
end

% Checking SINR constraints
size_exh2=size(exh2);
exh3 = zeros(M,N,10^2); 
count3=1;
Ha = Hrici(1,:);
Hb = Hrici(2,:);
for i=1:size_exh2(3)
    temp=exh2(:,:,i);
    if ((((Ha*temp(:,1))^2) / ((((Ha*temp(:,2))^2)+sigma2(1)))>=sinrmin(1))&&...
        (((Hb*temp(:,2))^2) / ((((Hb*temp(:,1))^2)+sigma2(2)))>=sinrmin(2)))
        exh3(:,:,count3)=exh2(:,:,i);
        count3=count3+1;
    end
 end

% Generating the objective function value for the set of feasible solutions
size_exh3=size(exh3);
exh4 = zeros(2,1); 
for i=1:size_exh3(3)
     temp=exh3(:,:,i);
    exh4(i)= sum(sum(temp.^2));
end
[min_val, ind]= min(exh4);
sol_exh= exh3(:,:,ind);
min_pow= min_val;
sinr_1exh = ((Ha*sol_exh(:,1))^2)/(((Ha*sol_exh(:,2))^2)+1);           % SINR achieved for user 1
sinr_2exh = ((Hb*sol_exh(:,2))^2)/(((Hb*sol_exh(:,1))^2)+1);           % SINR achieved for user 2
SINR = [sinr_1exh,sinr_2exh];
w = sol_exh;
end

%% Power allocation. SEMI-DEFINITE PROGRAMMING
function [w, SINR] = sdp_solver(Hrici, M, Pt_max, sinrmin, sigma2)
    % Configure parameters
    SINR = zeros(1,size(Hrici,1)); 
    w = zeros(size(Hrici,2),size(Hrici,1));

    R_user1 = Hrici(1,:)'*Hrici(1,:);  % Autocorrelation matrices of user 1 channel link
    R_user2 = Hrici(2,:)'*Hrici(2,:);  % Autocorrelation matrices of user 2 channel link

    X_user1 = zeros(M,M);              % Inizialize optimization variable for problem (2) 
    X_user2 = zeros(M,M);              % Inizialize optimization variable for problem (2) 
    
    cvx_begin quiet

        variables X_user1(M,M) X_user2(M,M)

        minimize(trace(X_user1)+trace(X_user2)) % Maximize total SINR 

        subject to
        diag(X_user1)+diag(X_user2)<=Pt_max*ones(M,1);                                      %#ok. Constraint on total transmitted power per antenna 
        trace(R_user1*X_user1)-(sinrmin(1)*trace(R_user1*X_user2))>= sigma2(1)*sinrmin(1);  %#ok. Contraint on minimum SINR for user 1
        trace(R_user2*X_user2)-(sinrmin(2)*trace(R_user2*X_user1))>= sigma2(2)*sinrmin(2);  %#ok. Contraint on minimum SINR for user 2
        X_user1 == hermitian_semidefinite(M);                                               %#ok. X_user1=X_user1*, X21>=0
        X_user2 == hermitian_semidefinite(M);                                               %#ok. X_user2=X_user2*, X_user2>=0

    cvx_end
   
    % Checking whether the solution is rank 1
    eval1=diag(svd(X_user1));
    eval2=diag(svd(X_user2));
    if((sum((eval1(2:end)<10e-6))-sum((eval2(2:end)<10e-6)))==0) % Checking that rank(X21) && rank (X22) are both 1.
      sol_1 = sqrt(diag(X_user1));                                   % Weighs for user 1, for all transmitting antennas
      sol_2 = sqrt(diag(X_user2));                                   % Weighs for user 2, for all transmitting antennas
      w =[sol_1,sol_2];
      Ha = Hrici(1,:);
      Hb = Hrici(2,:);
      sinr_1test = ((Ha*sol_1)^2)/(((Ha*sol_2)^2)+1);           % SINR achieved for user 1
      sinr_2test = ((Hb*sol_2)^2)/(((Hb*sol_1)^2)+1);           % SINR achieved for user 2
      SINR = [sinr_1test,sinr_2test];
      
    else 
      disp('Randomization required')                             % For not robust formulation rank = 1 [1].
    end 
end
  
%% REF [1] Mats Bengtsson, Björn Ottersten, Optimum and Suboptimum Transmit Beamforming