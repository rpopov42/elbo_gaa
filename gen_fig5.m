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

TN = 300; TN_eps = 50; w = 0.5; ug = 2; vg = 1; uc = 0; vc = 10; up = 0; vp = 10^2;

rspp = 1/sqrt(2*vp*pi)*exp(-1/2*(u - up).^2/vp);

x = [-1.3252, -2.3011, -5.0710, -1.4873,  0.4896,  1.1387,  2.9932,  2.7633,  1.4019,  2.6568, ...
      0.5034,  1.6761,  1.0904, -3.4887,  1.3403,  2.3313, -0.4575, -1.1529, -1.7524,  0.3278]'; % skewed

rspg = 1/sqrt(2*vg*pi)*exp(-1/2*(repmat(u,[length(x),1]) - repmat(x,[1,length(u)])).^2/vg);

rspc = 1/sqrt(2*vc*pi)*exp(-1/2*(x - uc).^2/vc);

rspu = prod((1 - w)*rspg + w*repmat(rspc,[1,length(u)]),1).*rspp;

logrspu = sum(log((1 - w)*rspg + w*repmat(rspc,[1,length(u)])),1) - 1/2*log(2*vp*pi) - 1/2*(u - up).^2/vp;

loge_true = log(sum(rspu)*du);

uq0 = mean(x); vq0 = mean((x - uq0).^2) + vg;

[uq_lbnu,vq_lbnu] = elbo_numerical_em(u,logrspu,uq0,vq0,TN); [elbo_lbnu] = elbo_numerical(u,logrspu,uq_lbnu,vq_lbnu);

tic(); [uq_lap,vq_lap] = laplace_serial(x,vg,up,vp,w,rspc,TN); time_laps = toc(); [elbo_laps] = elbo_numerical(u,logrspu,uq_lap,vq_lap);
tic(); [uq_lbmf,vq_lbmf] = elbo_mean_field_serial(x,vg,up,vp,w,rspc,TN); time_lbmfs = toc(); [elbo_lbmfs] = elbo_numerical(u,logrspu,uq_lbmf,vq_lbmf);
tic(); [uq_eps,vq_eps] = ep_serial(x,vg,up,vp,w,rspc,TN); time_eps = toc(); [elbo_eps] = elbo_numerical(u,logrspu,uq_eps,vq_eps);
tic(); [uq_lbga,vq_lbga] = elbo_gaa_em_serial(x,vg,up,vp,w,rspc,uq0,vq0,TN); time_lbgas = toc(); [elbo_lbgas] = elbo_numerical(u,logrspu,uq_lbga,vq_lbga);

uq = uq0; vq = vq0;

step = 1.0; step_decay = 0.97; M = 3; bufuq = zeros(1,TN); bufvq = zeros(1,TN);

%u0 = randn(TN,M);

u0 = get_rand_u0; % TN = 300, M = 3, step = 1.0; step_decay = 0.97

p = zeros(length(x),1);

tic();

for tn = 1:TN
  dspu = 0; dspv = 0;

  for m = 1:M
    for n = 1:length(x)
      rsp = 1/sqrt(2*vg*pi)*exp(-1/2*(u0(tn,m) - (x(n) - uq)/sqrt(vq))^2*vq/vg);

      p(n) = (1 - w)*rsp/((1 - w)*rsp + w*rspc(n));
    end;

    dspu = dspu + sum(p.*(x - (sqrt(vq)*u0(tn,m) + uq)))/vg + (up - uq)/vp;

    dspv = dspv - 1/2*(sum(p.*(u0(tn,m) - (x - uq)/sqrt(vq))*u0(tn,m))/vg - 1/vq + 1/vp);
  end;

  uq = uq + step*dspu/M; vq = min(max(vq + step*dspv/M,vq/2),2*vq);

  step = step*step_decay;

  bufuq(tn) = uq; bufvq(tn) = vq;
end;

time_stochs = toc();

[elbo_stochs] = elbo_numerical(u,logrspu,bufuq,bufvq);

hfig = figure;

subplot(1,2,1); hold on;
loglog(time_laps/TN*(1:TN),loge_true - elbo_laps,'LineWidth',2.0,'color','blue');
loglog(time_lbmfs/TN*(1:TN),loge_true - elbo_lbmfs,'LineWidth',2.0,'color',[0,0.9,0]);
loglog(time_eps/TN*(1:TN),loge_true - elbo_eps,'LineWidth',2.0,'color','red');
loglog(time_lbgas/TN*(1:TN),loge_true - elbo_lbgas,'LineWidth',2.0,'color',[0,0.8,0.8]);
loglog(time_stochs/TN*(1:TN),loge_true - elbo_stochs,'LineWidth',2.0,'color','magenta');
loglog(time_stochs/TN*(1:TN),loge_true - elbo_lbnu(end)*ones(1,TN),'--','color','black');

%xlabel ('\textbf{time,s}','interpreter','latex'); ylabel ('$\boldsymbol{\text{KL}(q(\mu)||p(\mu|X))}$','interpreter','latex');
xlabel ('\bf time,s'); ylabel ('\bf KL(q(\mu)||p(\mu|X))');

set(gca,'FontSize',12,'LabelFontSizeMultiplier',1.2,'LineWidth',1.0);

hlegend = legend('Laplace','MF','EP','ELBO GAA','SGD','Numerical ELBO');
%set(hlegend,'box','off','FontSize',16,'LineWidth',1.0,'location','southoutside','orientation','horizontal','interpreter','latex');
set(hlegend,'box','off','FontSize',16,'LineWidth',1.0,'location','southoutside','orientation','horizontal');

set(hlegend,'position',get(hlegend,'position') + [0.2,-0.02,0,0]);

tic(); [uq_lap,vq_lap] = laplace(x,vg,up,vp,w,rspc,TN); time_lap = toc(); [elbo_lap] = elbo_numerical(u,logrspu,uq_lap,vq_lap);
tic(); [uq_lbmf,vq_lbmf] = elbo_mean_field(x,vg,up,vp,w,rspc,TN); time_lbmf = toc(); [elbo_lbmf] = elbo_numerical(u,logrspu,uq_lbmf,vq_lbmf);
tic(); [uq_eps,vq_eps] = ep_serial(x,vg,up,vp,w,rspc,TN_eps); time_eps = toc(); [elbo_eps] = elbo_numerical(u,logrspu,uq_eps,vq_eps);
tic(); [uq_lbga,vq_lbga] = elbo_gaa_em(x,vg,up,vp,w,rspc,uq0,vq0,TN); time_lbga = toc(); [elbo_lbga] = elbo_numerical(u,logrspu,uq_lbga,vq_lbga);

uq = uq0; vq = vq0;

step = 1.0; step_decay = 0.97; M = 3; bufuq = zeros(1,TN); bufvq = zeros(1,TN);

%u0 = randn(TN,M);

u0 = get_rand_u0; % TN = 300, M = 3, step = 1.0; step_decay = 0.97

tic();

for tn = 1:TN
  dspu = 0; dspv = 0;

  for m = 1:M
    rsp = 1/sqrt(2*vg*pi)*exp(-1/2*(u0(tn,m) - (x - uq)/sqrt(vq)).^2*vq/vg);

    p = (1 - w)*rsp./((1 - w)*rsp + w*rspc);

    dspu = dspu + sum(p.*(x - (sqrt(vq)*u0(tn,m) + uq)))/vg + (up - uq)/vp;

    dspv = dspv - 1/2*(sum(p.*(u0(tn,m) - (x - uq)/sqrt(vq))*u0(tn,m))/vg - 1/vq + 1/vp);
  end;

  uq = uq + step*dspu/M; vq = min(max(vq + step*dspv/M,vq/2),2*vq);

  step = step*step_decay;

  bufuq(tn) = uq; bufvq(tn) = vq;
end;

time_stoch = toc();

[elbo_stoch] = elbo_numerical(u,logrspu,bufuq,bufvq);

subplot(1,2,2); hold on;
loglog(time_lap/TN*(1:TN),loge_true - elbo_lap,'LineWidth',2.0,'color','blue');
loglog(time_lbmf/TN*(1:TN),loge_true - elbo_lbmf,'LineWidth',2.0,'color',[0,0.9,0]);
loglog(time_eps/TN_eps*(1:TN_eps),loge_true - elbo_eps,'LineWidth',2.0,'color','red');
loglog(time_lbga/TN*(1:TN),loge_true - elbo_lbga,'LineWidth',2.0,'color',[0,0.8,0.8]);
loglog(time_stoch/TN*(1:TN),loge_true - elbo_stoch,'LineWidth',2.0,'color','magenta');
loglog(time_stoch/TN*(1:TN),loge_true - elbo_lbnu(end)*ones(1,TN),'--','color','black');

%xlabel ('\textbf{time,s}','interpreter','latex'); ylabel ('$\boldsymbol{\text{KL}(q(\mu)||p(\mu|X))}$','interpreter','latex');
xlabel ('\bf time,s'); ylabel ('\bf KL(q(\mu)||p(\mu|X))');

set(gca,'FontSize',12,'LabelFontSizeMultiplier',1.2,'LineWidth',1.0);

set(hfig,'papersize',[12,5.0],'paperposition',[-1.0,0.1,14.0,5.0]);

%print('fig4.pdf');




























