%  MIT License
%
%  Copyright (c) 2024 Roumen Popov
%
%  Permission is hereby granted, free of charge, to any person obtaining a copy
%  of this software and associated documentation files (the "Software"), to deal
%  in the Software without restriction, including without limitation the rights
%  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%  copies of the Software, and to permit persons to whom the Software is
%  furnished to do so, subject to the following conditions:
%
%  The above copyright notice and this permission notice shall be included in all
%  copies or substantial portions of the Software.
%
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%  SOFTWARE.

%  The code is meant primarily for research purposes and has not been optimised
%  or tested for numerical stability. Use it with caution!

clear;

graphics_toolkit('qt'); % remove if not using Octave

du = 0.01; u = (-40:du:40);

TN = 30; w = 0.5; ug = 2; vg = 1; uc = 0; vc = 10; up = 0; vp = 10^2;

rspp = 1/sqrt(2*vp*pi)*exp(-1/2*(u - up).^2/vp);

%N = 20; z = rand(N,1) < (1 - w); x = z.*(ug + sqrt(vg)*randn(N,1)) + (1 - z).*(uc + sqrt(vc)*randn(N,1));

%x = [-0.7108, 1.7009, -0.9289, -2.6814, -3.5422, 1.7143, 1.1686, 1.0208, -3.9754, 1.4664, ...
%     -0.5582, 2.9642,  2.5201,  1.9800,  1.9652, 1.0547, 1.2376, 1.4283,  1.2855, 0.5809]';  % normal

x = [-1.3252, -2.3011, -5.0710, -1.4873,  0.4896,  1.1387,  2.9932,  2.7633,  1.4019,  2.6568, ...
      0.5034,  1.6761,  1.0904, -3.4887,  1.3403,  2.3313, -0.4575, -1.1529, -1.7524,  0.3278]'; % skewed

rspg = 1/sqrt(2*vg*pi)*exp(-1/2*(repmat(u,[length(x),1]) - repmat(x,[1,length(u)])).^2/vg);

rspc = 1/sqrt(2*vc*pi)*exp(-1/2*(x - uc).^2/vc);

rspu = prod((1 - w)*rspg + w*repmat(rspc,[1,length(u)]),1).*rspp;

logrspu = sum(log((1 - w)*rspg + w*repmat(rspc,[1,length(u)])),1) - 1/2*log(2*vp*pi) - 1/2*(u - up).^2/vp;

loge_true = log(sum(rspu)*du);

uq0 = mean(x); vq0 = mean((x - uq0).^2) + vg;

[uq_lbnu,vq_lbnu] = elbo_numerical_em(u,logrspu,uq0,vq0,TN); [elbo_lbnu] = elbo_numerical(u,logrspu,uq_lbnu,vq_lbnu);
[uq_lap,vq_lap] = laplace(x,vg,up,vp,w,rspc,TN); [elbo_lap] = elbo_numerical(u,logrspu,uq_lap,vq_lap);
[uq_lbmf,vq_lbmf] = elbo_mean_field(x,vg,up,vp,w,rspc,TN); [elbo_lbmf] = elbo_numerical(u,logrspu,uq_lbmf,vq_lbmf);
[uq_eps,vq_eps] = ep_serial(x,vg,up,vp,w,rspc,TN); [elbo_eps] = elbo_numerical(u,logrspu,uq_eps,vq_eps);
[uq_lbga,vq_lbga] = elbo_gaa_em(x,vg,up,vp,w,rspc,uq0,vq0,TN); [elbo_lbga] = elbo_numerical(u,logrspu,uq_lbga,vq_lbga);

u = (-1:du:3);

rspg = 1/sqrt(2*vg*pi)*exp(-1/2*(repmat(u,[length(x),1]) - repmat(x,[1,length(u)])).^2/vg);

logrspu = sum(log((1 - w)*rspg + w*repmat(rspc,[1,length(u)])),1) - 1/2*log(2*vp*pi) - 1/2*(u - up).^2/vp - loge_true;

logrspu_lap = -1/2*((u - uq_lap(end)).^2/vq_lap(end) + log(2*vq_lap(end)*pi));
logrspu_lbmf = -1/2*((u - uq_lbmf(end)).^2/vq_lbmf(end) + log(2*vq_lbmf(end)*pi));
logrspu_eps = -1/2*((u - uq_eps(end)).^2/vq_eps(end) + log(2*vq_eps(end)*pi));
logrspu_lbga = -1/2*((u - uq_lbga(end)).^2/vq_lbga(end) + log(2*vq_lbga(end)*pi));
logrspu_lbnu = -1/2*((u - uq_lbnu(end)).^2/vq_lbnu(end) + log(2*vq_lbnu(end)*pi));

hfig = figure;

subplot(1,3,1); hold on;
plot(u,exp(logrspu_lap),'LineWidth',2.0,'color','blue');
plot(u,exp(logrspu_lbmf),'LineWidth',2.0,'color',[0,0.9,0]);
plot(u,exp(logrspu_eps),'LineWidth',2.0,'color','red');
plot(u,exp(logrspu_lbga),'LineWidth',2.0,'color',[0,0.8,0.8]);
plot(u,exp(logrspu),'LineWidth',2.0,'color','magenta');
plot(u,exp(logrspu_lbnu),'LineWidth',2.0,'--','color','black');

%xlabel ('$\boldsymbol{\mu}$','interpreter','latex'); ylabel ('$\boldsymbol{p(\mu|X), q(\mu)}$','interpreter','latex');
xlabel ('\bf \mu'); ylabel ('\bf p(\mu|X), q(\mu)');

set(gca,'FontSize',12,'LabelFontSizeMultiplier',1.2,'LineWidth',1.0);

subplot(1,3,2); hold on;
semilogy(loge_true - elbo_lap,'LineWidth',2.0,'color','blue');
semilogy(loge_true - elbo_lbmf,'LineWidth',2.0,'color',[0,0.9,0]);
semilogy(loge_true - elbo_eps,'LineWidth',2.0,'color','red');
semilogy(loge_true - elbo_lbga,'LineWidth',2.0,'color',[0,0.8,0.8]);
semilogy(loge_true - elbo_lbnu(end),'LineWidth',2.0,'color','magenta');
semilogy((loge_true - elbo_lbnu(end))*ones(1,length(elbo_lbnu)),'--','LineWidth',2.0,'color','black');

%xlabel ('\textbf{iteration}','interpreter','latex'); ylabel ('$\boldsymbol{\text{KL}(q(\mu)||p(\mu|X))}$','interpreter','latex');
xlabel ('\bf iteration'); ylabel ('\bf KL(q(\mu)||p(\mu|X))');

set(gca,'FontSize',12,'LabelFontSizeMultiplier',1.2,'LineWidth',1.0);

hlegend = legend('Laplace','MF','EP','ELBO GAA','Exact posterior','Numerical ELBO');

%set(hlegend,'FontSize',16,'LineWidth',1.0,'location','southoutside','orientation','horizontal','interpreter','latex');
set(hlegend,'FontSize',16,'LineWidth',1.0,'location','southoutside','orientation','horizontal');

subplot(1,3,3); hold on;
semilogy(abs(uq_lap - uq_lbnu(end)),'LineWidth',2.0,'color','blue');
semilogy(abs(uq_lbmf - uq_lbnu(end)),'LineWidth',2.0,'color',[0,0.9,0]);
semilogy(abs(uq_eps - uq_lbnu(end)),'LineWidth',2.0,'color','red');
semilogy(abs(uq_lbga - uq_lbnu(end)),'LineWidth',2.0,'color',[0,0.8,0.8]);

%xlabel ('\textbf{iteration}','interpreter','latex'); ylabel ('$\boldsymbol{|\mu_q - \bar\mu_q|}$','interpreter','latex');
xlabel ('\bf iteration'); ylabel ('\bf |\mu_q - \_\mu_q|');

set(gca,'FontSize',12,'LabelFontSizeMultiplier',1.2,'LineWidth',1.0);

set(hfig,'papersize',[12,4.5],'paperposition',[-1.0,0.1,14.0,4.5]);

%print('fig3.pdf');







