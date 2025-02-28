close all 
clear all

% tn = sym('tn');
% tn_1 = sym('tn_1');
% n = sym("n");
% tau = sym("tau");
% t_prime = sym('t_prime');
% q = sym('q', [1 2]);
% 
% phase = [1 (tn-t_prime); 0 1];
% Q = [q(1) 0; 0 q(2)];
% 
% expr = phase*Q*phase.';
% 
% expr_int = int(expr, tn_1, tn);
% 
% expr_int_sub = subs(expr_int,[tn, tn_1],[n*tau, (n-1)*tau]);
% simplify(expr_int_sub)

rng('default')                                % set the state of randn
T=1;N=500; dt=T/N; t = [dt:dt:1];

dW = sqrt(dt)*randn(1,N);   % increments
W = cumsum(dW);             % cumulative sum
plot([0:dt:T],[0,W],'r-')   % plot W against t
xlabel('t','FontSize',13)
ylabel('W(t)','FontSize',13)

% dX_t = mu(t)*dt + V(t)dWt
X1 = bm([0.1;0.1],[0.1;0.8]) % (A=mu, sigma)
nPeriods = 500;      % # of simulated observations
dt       =   1;      % time increment = 1 day

rng('default')                % make output reproducible
[X, T] = X1.simulate(nPeriods, 'DeltaTime', dt);

plot(T, X(:,:))
xlabel('Time (seconds)'), ylabel('X(t)')
			
%legend({'Primary Path' 'Antithetic Path'}, ...
%			'Location', 'Best')



%expr_int_full = int(expr, n*tau, (n-1)*tau)
%simplify(expr_int_full)
