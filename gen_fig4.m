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

TN = 20; w = 0.5; ug = 2; vg = 1; uc = 0; vc = 10; up = 0; vp = 10^2;

rspp = 1/sqrt(2*vp*pi)*exp(-1/2*(u - up).^2/vp);

%N = 5; z1 = rand(N,1) < (1 - w); x1 = z1.*(ug + sqrt(vg)*randn(N,1)) + (1 - z1).*(uc + sqrt(vc)*randn(N,1));
%N = 10; z2 = rand(N,1) < (1 - w); x2 = z2.*(ug + sqrt(vg)*randn(N,1)) + (1 - z2).*(uc + sqrt(vc)*randn(N,1));
%N = 100; z3 = rand(N,1) < (1 - w); x3 = z3.*(ug + sqrt(vg)*randn(N,1)) + (1 - z3).*(uc + sqrt(vc)*randn(N,1));

x1 = [-3.6048,  1.1717,  2.9157,  0.3442,  2.2250]';
x2 = [ 3.7161,  1.8442,  1.5167, -1.7605,  1.7936, -2.7840,  1.0120,  1.6812,  1.1503,  1.1121]';
x3 = [-3.246070,  0.619798,  2.323727, -7.951559,  1.667643,  1.204077,  1.486479,  1.644630,  0.884827,  2.629964,...
       1.006960, -0.847871,  2.753053,  2.303124,  1.078538, -0.042687,  4.872329,  1.038326,  2.789549,  0.819441,...
       2.070223,  3.711730, -3.105390, -2.092618,  2.254361,  2.973100,  1.633393, -2.993167,  2.757961,  0.934959,...
       2.319979,  0.708533,  0.927233,  1.102534,  1.019905,  1.140254,  2.037188,  1.046088,  1.871257,  4.918475,...
      -2.913199,  0.514955,  2.413798, -2.349509,  1.803889,  3.730963,  2.631266, -1.342075,  1.421190,  2.495447,...
       1.038637, -5.680024,  1.891150,  2.104155,  1.442317, -1.476125,  0.884525,  0.624891,  1.149958,  1.840244,...
       3.807929, -5.403366,  2.826735,  0.602245, -5.997698,  0.097400, -0.143678,  1.413421,  0.172258,  1.718791,...
      -2.013564,  5.619807,  0.398621, -1.518153,  2.886925,  1.375935,  2.404110,  2.290330,  2.737141, -0.650202,...
       2.933960, -1.592842,  2.446137,  3.241905,  7.784414,  1.473917,  1.759823,  2.825907,  1.208771, -2.653372,...
       6.571994,  1.594241,  1.220390, -1.030845,  2.682349, -0.390644, -2.506642,  3.381637, -1.336332,  2.767432]';

X{1} = x1; X{2} = x2; X{3} = x3;

hfig = figure;

for n = 1:3
  x = X{n};

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

  subplot(1,3,n); hold on;
  semilogy(loge_true - elbo_lap,'LineWidth',2.0,'color','blue');
  semilogy(loge_true - elbo_lbmf,'LineWidth',2.0,'color',[0,0.9,0]);
  semilogy(loge_true - elbo_eps,'LineWidth',2.0,'color','red');
  semilogy(loge_true - elbo_lbga,'LineWidth',2.0,'color',[0,0.8,0.8]);
  semilogy((loge_true - elbo_lbnu(end))*ones(1,length(elbo_lbnu)),'--','LineWidth',2.0,'color','black');

%  xlabel ('\textbf{iteration}','interpreter','latex'); ylabel ('$\boldsymbol{\text{KL}(q(\mu)||p(\mu|X))}$','interpreter','latex');
  xlabel ('\bf iteration'); ylabel ('\bf KL(q(\mu)||p(\mu|X))');

  set(gca,'FontSize',12,'LabelFontSizeMultiplier',1.2,'LineWidth',1.0);

  if n == 2
    hlegend = legend('Laplace','MF','EP','ELBO GAA','Numerical ELBO');
%    set(hlegend,'FontSize',16,'LineWidth',1.0,'location','southoutside','orientation','horizontal','interpreter','latex');
    set(hlegend,'FontSize',16,'LineWidth',1.0,'location','southoutside','orientation','horizontal');
  end;
end;

set(hfig,'papersize',[12,4.5],'paperposition',[-1.0,0.1,14.0,4.5]);

%print('fig4.pdf');







