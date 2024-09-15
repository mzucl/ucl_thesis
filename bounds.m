% Some definitions
sigma = @(a) exp(a)./(exp(a)+1);
lse   = @(a)log(exp(a)+1);

% A test to make sure the two answers are the same
a = 0.9;
d = 1;
log(sigma(a))*d + log(1-sigma(a))*(1-d)
a*d - lse(a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration of Laplace approximation and Bohnoing and Jaakola bounds
% Mean and variance of w*z
mu_a = -1.5; % <w>*<z>
sig2_a = 0.7;

a = (-5:0.1:5)';

% Laplace approximation
c  = lse(mu_a);
gl = exp(mu_a - lse(mu_a));    % Gradient
hl = 1/(2*(cosh(mu_a) + 1));  % for Laplace approximation
la = c + (a-mu_a)*gl + (a-mu_a).^2*hl/2;

% Bohning bound
c  = lse(mu_a);
gb = exp(mu_a - lse(mu_a));
hb = 0.25;
bb = c + (a-mu_a)*gb + (a-mu_a).^2*hb/2;

% Jaakola bound
lam = (1/(1+exp(-mu_a)) - 0.5)/(2*mu_a);
c   = lse(mu_a);
gj  = 0.5 + 2*lam*mu_a;
hj  = 2*lam;
jb  = c + (a-mu_a)*gj + (a-mu_a).^2*hj/2;

plot(a,[lse(a) la bb jb])
legend({'LSE(a)','Laplace', 'Bohning','Jaakola'})
% Notice that the Laplace approximation can fall below the function
% so it isn't a bound.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expectations computed in different ways
% E[lse(a)] \simeq E[lse(mu_a) + (a-mu_a)*g + (a-mu_a).^2*h/2]
% = lse(mu_a) + sig2_a*h/2;

Nsamp = 10000;
a = mu_a + randn(Nsamp,1)*sqrt(sig2_a); % Random samples
[mean(lse(a)), lse(mu_a), ...
 lse(mu_a) + sig2_a*hl/2, ...
 lse(mu_a) + sig2_a*hj/2, ...
 lse(mu_a) + sig2_a*hb/2]