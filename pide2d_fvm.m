clear; close all;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

%% settings

measure = 'no'; % 'no', 'yes'

if strcmp(measure, 'yes')
    pics = 'no'; % 'no', 'yes'
    runs = 100;
    Ts = zeros(1, runs);
    types = {'act', 'rep', 'mix', 'baci'};
    for ix = 1:length(types)
        type = types{ix};
        for n = [100, 200, 300, 400]
            for i = 1:runs
                start = tic;
                PIDE(type, pics, n, n);
                Ts(i) = toc(start);
            end
            fprintf('%s with %dx%d took %.4f secs\n', type, n, n, mean(Ts))
        end
    end
else
    type = 'baci'; % 'act', 'rep', 'mix', 'baci'
    pics = 'yes'; % 'no', 'yes'
    n = 100;
    tic
    P_eq = PIDE(type, pics, n, n);
    toc
end

%% params

function [P_eq, nz] = PIDE(type, pics, n1, n2)
if strcmp(type, 'act')
    L1 = 350;
    L2 = 350;
    km1 = 3.4e-3;
    km2 = 3.4e-3;
    gx1 = @(x,y) 4e-4+0*x+0*y;
    gx2 = @(x,y) 4e-4+0*x+0*y;
    b1 = 18;
    b2 = 18;
    H12 = -4;
    H21 = -4;
    K1 = 70;
    K2 = 70;
    eps1 = 0.2;
    eps2 = 0.2;

    c1 = @(x, y) (K1^H12+eps1*y.^H12)./(K1^H12+y.^H12);
    c2 = @(x, y) (K2^H21+eps2*x.^H21)./(K2^H21+x.^H21);

    pdflim = 1e-4;
    mpdflim1 = 0.01;
    mpdflim2 = 0.01;
    barlim = [0 8e-5];

    tau = 50;
    T = tau/gx1(0,0);
elseif strcmp(type, 'rep')
    L1 = 400;
    L2 = 400;
    km1 = 3.2e-3;
    km2 = 3.2e-3;
    gx1 = @(x,y) 4e-4+0*x+0*y;
    gx2 = @(x,y) 4e-4+0*x+0*y;
    b1 = 16;
    b2 = 16;
    H12 = 4;
    H21 = 4;
    K1 = 45;
    K2 = 45;
    eps1 = 0.15;
    eps2 = 0.15;

    c1 = @(x, y) (K1^H12+eps1*y.^H12)./(K1^H12+y.^H12);
    c2 = @(x, y) (K2^H21+eps2*x.^H21)./(K2^H21+x.^H21);

    pdflim = 2e-4;
    mpdflim1 = 0.02;
    mpdflim2 = 0.02;
    barlim = [0 13e-5];

    tau = 50;
    T = tau/gx1(0,0);
elseif strcmp(type, 'mix')
    L1 = 150;
    L2 = 200;
    km1 = 4e-3;
    km2 = 8e-3;
    gx1 = @(x,y) 4e-4+0*x+0*y;
    gx2 = @(x,y) 4e-4+0*x+0*y;
    b1 = 10;
    b2 = 20;
    H11 = -4;
    H12 = 2;
    H21 = -6;
    H22 = 2;
    K11 = 45;
    K12 = 45;
    K21 = 70;
    K22 = 70;
    eps11 = 0.002;
    eps21 = 0.002;
    eps12 = 0.02;
    eps22 = 0.1;
    eps13 = 0.2;
    eps23 = 0.2;

    c1 = @(x, y) (eps11*x.^H11.*y.^H12+eps12*K11^H11*y.^H12+eps13*x.^H11*K12^H12+K11^H11*K12^H12)...
        ./(x.^H11.*y.^H12+K11^H11*y.^H12+x.^H11*K12^H12+K11^H11*K12^H12);
    c2 = @(x, y) (eps21*x.^H21.*y.^H22+eps23*K21^H21*y.^H22+eps22*x.^H21*K22^H22+K21^H21*K22^H22)...
        ./(x.^H21.*y.^H22+K21^H21*y.^H22+x.^H21*K22^H22+K21^H21*K22^H22);

    pdflim = 2e-3;
    mpdflim1 = 0.1;
    mpdflim2 = 0.04;
    barlim = [0 10e-4];

    tau = 50;
    T = tau/gx1(0,0);
elseif strcmp(type, 'baci')
    L1 = 200;
    L2 = 600;
    alpha = 2.8e-3;
    beta_k = 4.9e-2;
    beta_s = 5.7e-2;
    delta_k = 1.4e-3;
    delta_s = 1.4e-3;
    Gamma_k = 500;
    Gamma_s = 50;
    b1 = 2;
    b2 = 5;
    km1 = (alpha+beta_k)/b1;
    km2 = beta_s/b2;
    gx1 = @(x,y) delta_k*Gamma_k*Gamma_s./(Gamma_k*Gamma_s+Gamma_s*x+Gamma_k*y);
    gx2 = @(x,y) delta_s*Gamma_k*Gamma_s./(Gamma_k*Gamma_s+Gamma_s*x+Gamma_k*y);
    H11 = -2;
    H21 = 5;
    K1 = 100;
    K2 = 110;
    eps1 = alpha/(alpha+beta_k);
    eps2 = 0;

    c1 = @(x, y) (K1^H11+eps1*x.^H11)./(K1^H11+x.^H11);
    c2 = @(x, y) (K2^H21+eps2*x.^H21)./(K2^H21+x.^H21);

    pdflim = 4e-4;
    mpdflim1 = 0.1;
    mpdflim2 = 0.01;
    barlim = [0 3e-4];

    tau = 200;
    T = tau/delta_k;
end

thresh = 1e-8; % decrease if not precise, increase if too slow

maxiter = 5000;

xi1p2 = linspace(0, L1, n1+1);
xi = conv(xi1p2, [1 1]/2, 'valid');

yi1p2 = linspace(0, L2, n2+1);
yi = conv(yi1p2, [1 1]/2, 'valid');

dx = diff(xi);
dx = dx(1);
dy = diff(yi);
dy = dy(1);
xim = xi-dx/2;
yim = yi-dy/2;
h = dx*dy;

dirac = 1e-2;
beta1 = @(x) exp(-x/b1)/b1-exp(-(x/dirac).^2)/(dirac*sqrt(pi)/2);
beta2 = @(x) exp(-x/b2)/b2-exp(-(x/dirac).^2)/(dirac*sqrt(pi)/2);

%% gamma, C

Ntrapz = 10;
[Xi, Yi] = meshgrid(xi-dx/2+eps,yi-dy/2+eps);
[X, Y] = meshgrid(linspace(0,dx,Ntrapz),linspace(0,dy,Ntrapz));
g1 = squeeze(trapz(trapz(gx1(reshape(Xi,[1 1 n1 n2])+X,...
    reshape(Yi,[1 1 n1 n2])+Y'))))*dx*dy/(Ntrapz-1)^2/h;
g2 = squeeze(trapz(trapz(gx2(reshape(Xi,[1 1 n1 n2])+X,...
    reshape(Yi,[1 1 n1 n2])+Y'))))*dx*dy/(Ntrapz-1)^2/h;
C1 = reshape((squeeze(trapz(trapz(c1(reshape(Xi,[1 1 n1 n2])+X,...
    reshape(Yi,[1 1 n1 n2])+Y'))))*dx*dy/(Ntrapz-1)^2/h)', [n1*n2 1]);
C2 = reshape((squeeze(trapz(trapz(c2(reshape(Xi,[1 1 n1 n2])+X,...
    reshape(Yi,[1 1 n1 n2])+Y'))))*dx*dy/(Ntrapz-1)^2/h)', [n1*n2 1]);

%% b

Ntrapz = 100;
b1 = trapz(beta1(xi-dx+linspace(0,dx,Ntrapz)'))*dx/(Ntrapz-1);
% b1(1) = integral(beta1, 0, dx/2);
b1(1) = -sum(b1(2:end));
b2 = trapz(beta2(yi-dy+linspace(0,dy,Ntrapz)'))*dy/(Ntrapz-1);
% b2(1) = integral(beta2, 0, dy/2);
b2(1) = -sum(b2(2:end));

%% Gamma

% x degradation
tmp = reshape(g1'.*xim'/dx, [n1*n2 1]);
degx = spdiags([-tmp tmp], [0 1], n1*n2, n1*n2);

% y degradation
tmp = reshape((yim'.*g2/dy)', [n1*n2 1]);
degy = spdiags([-tmp tmp], [0 n1], n1*n2, n1*n2);

clear tmp;

% x burst
Bx = flip(sparse(hankel(flip(b1))));
Bx(end,:) = 0;
Bx = C1'.*kron(eye(n2),Bx);
Bx = km1*(Bx - kron(speye(n2),kron(ones(1,n1),flip(speye(n1,1)))).*sum(Bx));

% y burst
By = flip(sparse(hankel(flip(b2))));
By(end,:) = 0;
By = C2'.*kron(By,eye(n1));
By = km2*(By - kron(ones(1,n2),kron(flip(speye(n2,1)),speye(n1))).*sum(By));

% assembly
Gamma = Bx+By+degx+degy;
% clear Bx By degx degy;

nz = nnz(Gamma);

%% equilibrium

tildeGamma = [ones(1, n1*n2)*h; Gamma(2:end,:)];
%     [L, U] = ilu(tildeGamma);
%     [p_eq,flag,relres,iter] = gmres(tildeGamma, eye(n1*n2,1), [], thresh, maxiter, L, U);
[P, R, C] = equilibrate(tildeGamma);
tildeGamma = R*P*tildeGamma*C;
[L, U] = ilu(tildeGamma);
[p_eq,flag,relres,iter] = gmres(tildeGamma, R*P*eye(n1*n2,1), [], thresh, maxiter, L, U);
if flag ~= 0
    flag
    relres
    iter
end
P_eq = reshape(C*p_eq, [n1 n2])';

%% plot

if strcmp(pics, 'yes')
    close all;

    [X, Y] = meshgrid(xi, yi);
    figure;
    surf(X, Y, P_eq, 'EdgeColor', 'none');
    c = colorbar;
    clim(barlim);
    c.FontSize = 12;
    xlabel('$x_1$', 'FontSize', 14);
    ylabel('$x_2$', 'FontSize', 14);
    zlabel('$p(x_1,x_2)$', 'FontSize', 14);
    xlim([0 L1]);
    ylim([0 L2]);
    zlim([0 pdflim]);
    title('Stationary probability density function', 'FontSize', 14);
    view([-10 55]);

    figure;
    plot(xi, sum(P_eq)*L2/n2, 'LineWidth', 2);
    xlabel('$x_1$', 'FontSize', 14);
    ylabel('$\overline p(x_1)$', 'FontSize', 14);
    xlim([0 L1]);
    ylim([0 mpdflim1]);
    title('Stationary marginal probability density function', 'FontSize', 14);

    figure;
    plot(yi, sum(P_eq, 2)*L1/n1, 'LineWidth', 2);
    xlabel('$x_2$', 'FontSize', 14);
    ylabel('$\overline p(x_2)$', 'FontSize', 14);
    xlim([0 L2]);
    ylim([0 mpdflim2]);
    title('Stationary marginal probability density function', 'FontSize', 14);

    if strcmp(type, 'baci')
        figure;
        contour(X, Y, P_eq, linspace(0.5e-4, 4e-4, 50));
        hold on;
        contour(X, Y, P_eq, 4*logspace(-6.5,-5,10));
        xlabel('$x_1$', 'FontSize', 14);
        ylabel('$x_2$', 'FontSize', 14);
        zlabel('$\overline p(x_1,x_2)$', 'FontSize', 14);
        xlim([0 150]);
        ylim([0 500]);
        title('Contour of stationary probability density function', 'FontSize', 14);
    end
end
end