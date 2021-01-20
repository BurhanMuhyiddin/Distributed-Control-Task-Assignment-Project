clc; clear all; close all;

N = 3; % Number of agents
M = N; % Dimension of each agent

%[c, ~, ~, ~, ~, LB, UB] = get_constraint_matrices(N, M);
s1 = Simulation(N, 0, 40);
[s1, c, LB, UB] = s1.get_constraint_matrices();
A_eq = zeros(2*N,N*N);
ind = 0;
st = size(A_eq, 2) / N;
for ii=1:N
    A_eq(ii,ii+ind:ii*st) = ones(1, N);
    A_eq(N+1:end, ii+ind:ii*st) = eye(N);
    ind = ind + N-1;
end
b_eq = ones(N+M, 1);
%%
% Get central solution in order to calculate cost error
options = optimoptions('linprog','Display','none');
[~, fopt, exit_flag] = linprog(c,[],[],A_eq,b_eq,LB,UB,options);
%[~, fopt1, exit_flag] = linprog(c,[],[],A_eq,b_eq,LB,UB,options);

if exit_flag ~= 1
  fprintf(2,'A problem occurred in the centralized solution\n');
  return;
end

fprintf('Centralized optimal cost is %.4g\n',fopt);
%%

p = 0.1;
[AA_NW, AA] = binomialGraph(1, N, 'doubly');
%%

MAXITERS = 10e4;

% The size of cell array is NxMAXITERS. N is the number of agents
% and each cell is Mx1 vector. M is dimension of each agent.
% The cells are initiated with zero value
XX = cellmat(N,MAXITERS,M,1,0);
XX_RA = cellmat(N,MAXITERS,M,1,0);
LM = cellmat(N,MAXITERS,size(A_eq,1),1,0);
VV = zeros(size(A_eq,1), N);

primal_cost = zeros(MAXITERS, 1);
dual_cost   = zeros(MAXITERS, 1);
primal_cost_RA = zeros(MAXITERS, 1);
dual_cost_RA   = zeros(MAXITERS, 1);

consensus_err = zeros(MAXITERS, 1);

for tt = 1:MAXITERS-1
  if mod(tt,100)==0
      fprintf('Iteration n. %d\n',tt);
  end
  
  gamma_t = 0.1*(1/tt)^0.6; % 1/tt^alpha with alpha in (0.5, 1]

  for ii=1:N
    N_ii = find(AA_NW(:,ii) == 1)';

    % v_i =  sum_{k \in N_i U { i } } w_ik mu_k^t
    VV(:,ii) = AA(ii,ii) * LM{ii,tt};
    
    for kk = N_ii
      VV(:,ii) = VV(:,ii) + AA(ii,kk) * LM{kk,tt};
    end
  end

  % Primal Update
  % ii+ind:ii+ind+N takes into account dimension of each agent.
  % For example in this case this will result at each iteration...
  % 1 to 4; 5 to 8; 9 to 12
  ind = 0;
  for ii=1:N
     XX{ii,tt} = linprog(c(1,ii+ind:ii*st)+VV(:,ii)'*A_eq(:,ii+ind:ii*st), ...
         [],[],[],[],LB(ii+ind:ii*st,1),UB(ii+ind:ii*st,1),options);
     ind = ind + N - 1;
  end
  
  % Running average
  for ii=1:N
      if tt==1
          XX_RA{ii,tt} = XX{ii,tt};
      else
          XX_RA{ii,tt} = (1/tt)*((tt-1)*XX_RA{ii,tt-1}+XX{ii,tt});
      end
  end
  
  % Dual Update
  % mu^{t+1} = mu^t + gamma^t* ( sum_i grad_q_i(mu^t) )
  ind=0;
  for ii=1:N
    grad_ii = A_eq(:,ii+ind:ii*st) * XX{ii,tt}-b_eq/N;
    LM{ii,tt+1} = VV(:,ii) + gamma_t*grad_ii;  
    ind = ind + N - 1;
  end
  
  % Performance check
  LM_avg = mean(vertcat(LM{:,tt}),1);

  ind = 0;
  for ii=1:N
    ff_ii = c(1,ii+ind:ii*st) * XX{ii,tt};
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = c(1,ii+ind:ii*st) * XX_RA{ii,tt};
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = c(1,ii+ind:ii*st) * XX{ii,tt} ...
        + LM{ii,tt}'*(A_eq(:,ii+ind:ii*st) * XX{ii,tt}-b_eq/N);
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = c(1,ii+ind:ii*st) * XX_RA{ii,tt} ...
        + LM{ii,tt}'*(A_eq(:,ii+ind:ii*st) * XX_RA{ii,tt}-b_eq/N);
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;

    consensus_err(tt) = consensus_err(tt) + norm(vertcat(LM{ii,tt}) - LM_avg);
    ind = ind + N - 1;
  end
end
%%
% Last value [for plot]
tt = MAXITERS;
fprintf('Iteration n. %d\n',tt);
LM_avg = mean(vertcat(LM{:,tt}),1);

for ii=1:N
    N_ii = find(AA_NW(:,ii) == 1)';
  
    % v_i =  sum_{j \in N_i U { i } } w_ij mu_j^t
    VV(:,ii) = AA(ii,ii) * LM{ii,tt};
    
    for jj = N_ii
      VV(:,ii) = VV(:,ii) + AA(ii,jj) * LM{jj,tt};
    end
end

ind = 0;
for ii=1:N
    XX{ii,tt} = linprog(c(1,ii+ind:ii*st)+VV(:,ii)'*A_eq(:,ii+ind:ii*st),...
         [],[],[],[],LB(ii+ind:ii*st,1),UB(ii+ind:ii*st,1),options);
    XX_RA{ii,tt} = (1/tt)*((tt-1)*XX_RA{ii,tt-1}+XX{ii,tt});

    ff_ii = c(1,ii+ind:ii*st) * XX{ii,tt};
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = c(1,ii+ind:ii*st) * XX_RA{ii,tt};
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = c(1,ii+ind:ii*st) * XX{ii,tt} ...
        + LM{ii,tt}'*(A_eq(:,ii+ind:ii*st) * XX{ii,tt}-b_eq/N);
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = c(1,ii+ind:ii*st) * XX_RA{ii,tt} ...
        + LM{ii,tt}'*(A_eq(:,ii+ind:ii*st) * XX_RA{ii,tt}-b_eq/N);
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;
    
    ind = ind + N - 1;
end

%%
figure
  semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(primal_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  semilogy(1:MAXITERS,abs(dual_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('cost error')
  legend('primal cost','primal cost with running avg','dual cost','dual cost with running avg')
  
%%
% figure
%   plot(1:MAXITERS,LM, 'LineWidth',2);
%   hold on, grid on, zoom on
%   xlabel('t')
%   ylabel('lm_i^t')
% 
% %%
% figure
%   plot(1:MAXITERS,sum(XX,1)-bb_centr, 'LineWidth',2);
%   hold on, grid on, zoom on
%   plot(1:MAXITERS,sum(XX_RA,1)-bb_centr, 'LineWidth',2);
%   xlabel('t')
%   ylabel('x_1^t + ... + x_N^t - b')
%   legend('x','x from running avg')
%   
%   
% %%
% figure
%   semilogy(1:MAXITERS,consensus_err(1:MAXITERS), 'LineWidth',2);
%   hold on, grid on, zoom on
%   xlabel('t')
%   ylabel('consensus error')

figure
    s1.create_environment();
    hold on
    s1.draw_lines([XX{1, MAXITERS}; XX{2, MAXITERS}; XX{3, MAXITERS}]');