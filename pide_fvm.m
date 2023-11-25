%% settings

clear; close all;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

%% params

measure = 'no'; % 'yes', 'no'

if strcmp(measure, 'yes')
    pics = 'no'; % 'yes', 'no'
    for n = [800, 1200, 1600, 2000, 5000]
        runs = 100;
        Ts = zeros(1, runs);
        for i = 1:runs
            start = tic;
            PIDE(pics, n);
            Ts(i) = toc(start);
        end
        fprintf('n=%d took %.4f secs\n', n, mean(Ts))
    end
else
    pics = 'yes'; % 'yes', 'no'
    n = 1000;
    [~, b1, g] = PIDE(pics, n);
end

function [p_eq, b1, g] = PIDE(pics, n)
gx1 = @(x) 4e-4+0*x;
km1 = 3.2e-3;
b1 = 16;
H11 = -4;
K1 = 45;
eps1 = 0.15;

tau = 50;
T = tau/gx1(0);

dirac = 1e-2; % scale of Gaussian function

L1 = 350;

c1 = @(x) (K1^H11+eps1*x.^H11)./(K1^H11+x.^H11);
beta1 = @(x) exp(-x/b1)/b1-exp(-(x/dirac).^2)/(dirac*sqrt(pi)/2);

xi1p2 = linspace(0, L1, n+1);
xi = conv(xi1p2, [1 1]/2, 'valid');

%%

dx = diff(xi);
dx = dx(1);
xim = xi-dx/2;
Ntrapz = 10;

b1 = trapz(beta1(xi-dx+linspace(0,dx,Ntrapz)'))*dx/(Ntrapz-1);
c = trapz(c1(xi-dx/2+eps+linspace(0,dx,Ntrapz)'))/(Ntrapz-1);
g = trapz(gx1(xi-dx/2+eps+linspace(0,dx,Ntrapz)'))/(Ntrapz-1).*xim/dx^2;
c = c/dx;

b1(1) = integral(beta1, 0, dx/2);
B = toeplitz(b1, [b1(1) zeros(1, n-1)]);
B(n, :) = 0;
B(n, :) = -sum(B);

Gamma = diag(g(2:end),1)-diag(g)+km1*B*diag(c);

%% equilibrium

Gamma_eq = Gamma;
Gamma_eq(1, :) = 1;
b = zeros(n, 1);
b(1) = n/L1;
p_eq = linsolve(Gamma_eq, b);

%% plot

if strcmp(pics, 'yes')
    figure;
    plot(xi, p_eq, 'LineWidth', 2);
    xlabel('$x$', 'FontSize', 14);
    ylabel('$\overline p(x)$', 'FontSize', 14);
    title('Stationary probability distribution function', 'FontSize', 14);
    xlim([0 L1]);
    ylim([0 1e-2]);
end
end