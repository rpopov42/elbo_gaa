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

u = (-7:0.01:7); uq = (-2.0:0.01:3.2); x = [-1.0,0.5,0.9,1.0,1.15,1.4]'; vg = 0.5^2; w = 0.3; Pc = 1/(u(end) - u(1)); vqt = vg/5; up = 0.0; vp = 10.0^2;

g = 1/sqrt(2*vg*pi)*exp(-1/2*(repmat(u,[length(x),1]) - repmat(x,[1,length(u)])).^2/vg); rsp0 = (1 - w)*g + w*Pc; p0 = (1 - w)*g./rsp0;

Guq = zeros(1,length(uq)); hGuq = zeros(1,length(uq));

for n = 1:length(uq)
  q = 1/sqrt(2*vqt*pi)*exp(-1/2* (u - uq(n)).^2/vqt); dq_duq = q.*(u - uq(n))/vqt;

  Guq(n) = sum(dq_duq.*sum(log(rsp0)))*0.01 + (up - uq(n))/vp;

  rspuq = 1/sqrt(2*vg*pi)*exp(-1/2*vg*(x - uq(n)).^2/(vg + vqt)^2); p = (1 - w)*rspuq./((1 - w)*rspuq + w*Pc);

  v = vg./((1 - p).*(p.*vg.*(x - uq(n)).^2/(vg + vqt)^2 + 1)*vqt + vg);

  A = exp(-1/2*vqt*(1 - p.^2.*v).*(x - uq(n)).^2/(vg + vqt)^2);

  hGuq(n) = sum(p.*sqrt(v).*A.*(vg + p.*v*vqt)/(vg + vqt).*(x - uq(n))/vg) + (up - uq(n))/vp;
end;

uqt = -0.2;

rspt = 1/sqrt(2*vg*pi)*exp(-1/2*vg*(x - uqt).^2/(vg + vqt)^2); pt = (1 - w)*rspt./((1 - w)*rspt + w*Pc);

vt = vg./((1 - pt).*(pt.*vg.*(x - uqt).^2/(vg + vqt)^2 + 1)*vqt + vg);

At = exp(-1/2*vqt*(1 - pt.^2.*vt).*(x - uqt).^2/(vg + vqt)^2);

Bt = pt.*sqrt(vt).*At.*(vg + pt.*vt*vqt)/(vg + vqt);

hGuq_uqt = sum(Bt.*(x - uqt))/vg + (up - uqt)/vp;

bGuq = -sum(Bt.*(uq - x))/vg - (uq - up)/vp; uqt1 = (sum(Bt.*x)/vg + up/vp)/(sum(Bt)/vg + 1/vp);

hfig = figure;

subplot(1,2,1); hold on;
plot(uq(uq < 1.3),0*uq(uq < 1.3),'LineWidth',1.0,'--','color','black','HandleVisibility','off');
plot(uq(uq < 2 | Guq < -5.0),Guq(uq < 2 | Guq < -5.0),'LineWidth',2.0,'color','blue');
plot(uq(uq < 2 | hGuq < -5.0),hGuq(uq < 2 | hGuq < -5.0),'LineWidth',2.0,'color','red');
plot(uq,bGuq,'LineWidth',2.0,'color',[0,0.9,0]);
plot(uqt,hGuq_uqt,'+','MarkerSize',12,'LineWidth',2.0,'color','black');
plot(uqt1,0,'x','MarkerSize',10,'LineWidth',2.0,'color','black');

%xlim([uq(1),uq(end)]); xlabel ('$\boldsymbol{\mu_q}$','interpreter','latex');
xlim([uq(1),uq(end)]); xlabel ('\bf \mu_q');

%ylim([1.1*min(Guq),1.1*max(Guq)]); ylabel ('$\boldsymbol{\partial \mathcal{L}(q(\mu))/\partial\mu_q}$','interpreter','latex');
ylim([1.1*min(Guq),1.1*max(Guq)]); ylabel ('\bf \partial L(q(\mu))/\partial\mu_q');

set(gca,'FontSize',14,'LabelFontSizeMultiplier',1.3,'LineWidth',1.0);

%hlegend = legend('$G_{\mu_q}$','$\widetilde G_{\mu_q}$','$\widehat G_{\mu_q}$', ...
%                 '$\widetilde G_{\mu_q}(\mu_q^{(t)})$','$\mu_q^{(t+1)}$');
hlegend = legend('G_{\mu_q}','~G_{\mu_q}','\^G_{\mu_q}','~G_{\mu_q}(\mu_q^{(t)})','\mu_q^{(t+1)}');

%set(hlegend,'FontSize',18,'LineWidth',1.0,'interpreter','latex');
set(hlegend,'FontSize',18,'LineWidth',1.0);

u = (-10:0.01:10); vq = (0.0:0.01:1.2); x = [0.7,1.2]'; vg = 0.5^2; w = 0.3; Pc = 1/(u(end) - u(1)); uqt = 0.9; up = 0.0; vp = 10.0^2;

g = 1/sqrt(2*vg*pi)*exp(-1/2*(repmat(u,[length(x),1]) - repmat(x,[1,length(u)])).^2/vg); rsp0 = (1 - w)*g + w*Pc; p0 = (1 - w)*g./rsp0;

Gvq = zeros(1,length(vq)); hGvq = zeros(1,length(vq)); hGvqvq = zeros(1,length(vq));

for n = 1:length(vq)
  q = 1/sqrt(2*vq(n)*pi)*exp(-1/2*(u - uqt).^2/vq(n)); dq_dvq = 1/2*q.*((u - uqt).^2/vq(n)^2 - 1/vq(n));

  Gvq(n) = sum(dq_dvq.*sum(log(rsp0)))*0.01 + 1/2*1/vq(n) - 1/2*1/vp;

  rspvq = 1/sqrt(2*vg*pi)*exp(-1/2*vg*(x - uqt).^2/(vg + vq(n))^2); p = (1 - w)*rspvq./((1 - w)*rspvq + w*Pc);

  v = vg./((1 - p).*(p.*vg.*(x - uqt).^2/(vg + vq(n))^2 + 1)*vq(n) + vg);

  A = exp(-1/2*vq(n)*(1 - p.^2.*v).*(x - uqt).^2/(vg + vq(n))^2);
  B = p.*sqrt(v).*A.*(vg + p.*v*vq(n))/(vg + vq(n));
  C = p.*sqrt(v).*A.*v; D = (1 - p.*v).*B;

  hGvq(n) = -1/2*(sum(C)/vg - sum(D.*(x - uqt).^2)/vg/(vg + vq(n)) - 1/vq(n) + 1/vp);

  hGvqvq(n) = -1/2*((sum(C)/vg + 1/vp)*vq(n) - sum(D.*(x - uqt).^2)/vg*vq(n)/(vg + vq(n)) - 1);
end;

vqt = 0.85;

rspt = 1/sqrt(2*vg*pi)*exp(-1/2*vg*(x - uqt).^2/(vg + vqt)^2); pt = (1 - w)*rspt./((1 - w)*rspt + w*Pc);

vt = vg./((1 - pt).*(pt.*vg.*(x - uqt).^2/(vg + vqt)^2 + 1)*vqt + vg);

At = exp(-1/2*vqt*(1 - pt.^2.*vt).*(x - uqt).^2/(vg + vqt)^2);
Bt = pt.*sqrt(vt).*At.*(vg + pt.*vt*vqt)/(vg + vqt);
Ct = pt.*sqrt(vt).*At.*vt; Dt = (1 - pt.*vt).*Bt;

hGvqvq_vqt = -1/2*((sum(Ct)/vg + 1/vp)*vqt - (sum(Dt.*(x - uqt).^2)/vg*vqt/(vg + vqt) + 1));

bGvq = -1/2*((sum(Ct)/vg + 1/vp)*vq - (sum(Dt.*(x - uqt).^2)/vg*vqt/(vg + vqt) + 1));

vqt1 = (sum(Dt.*(x - uqt).^2)/vg*vqt/(vg + vqt) + 1)/(sum(Ct)/vg + 1/vp);

subplot(1,2,2); hold on;
plot(vq,0*vq,'LineWidth',1.0,'--','color','black','HandleVisibility','off');
plot(vq(vq > 0.05),Gvq(vq > 0.05),'LineWidth',2.0,'color','blue');
plot(vq(vq > 0.05),hGvq(vq > 0.05),'LineWidth',2.0,'color','red');
plot(vq,hGvqvq,'LineWidth',2.0,'color','magenta');
plot(vq,bGvq,'LineWidth',2.0,'color',[0,0.9,0.0]);
plot(vqt1,0,'x','MarkerSize',10,'LineWidth',2.0,'color','black');
plot(vqt,hGvqvq_vqt,'+','MarkerSize',12,'LineWidth',2.0,'color','black');

set(hfig,'papersize',[7,5.0],'paperposition',[-0.2,0.1,7.8,5.1]);

%xlim([0,vq(end)]); xlabel ('$\boldsymbol{v_q}$','interpreter','latex');
xlim([0,vq(end)]); xlabel ('\bf v_q');

%ylim([1.4*min(hGvq),max(hGvq(vq > 0.05))]); ylabel ('$\boldsymbol{\partial \mathcal{L}(q(\mu))/\partial v_q}$','interpreter','latex');
ylim([1.4*min(hGvq),max(hGvq(vq > 0.05))]); ylabel ('\bf \partial L(q(\mu))/\partial v_q');

set(gca,'FontSize',14,'LabelFontSizeMultiplier',1.3,'LineWidth',1.0);

%hlegend = legend('$G_{v_q}$','$\widetilde G_{v_q}$','$\widetilde G_{v_q}v_q$','$\widehat G_{v_q}$',...
%                 '$v_q^{(t+1)}$','$\widetilde G_{v_q}(v_q^{(t)})v_q^{(t)}$');
hlegend = legend('G_{v_q}','~G_{v_q}','~G_{v_q}v_q','\^G_{v_q}','v_q^{(t+1)}','~G_{v_q}(v_q^{(t)})v_q^{(t)}');

%set(hlegend,'FontSize',18,'LineWidth',1.0,'location','northeast','orientation','horizontal','numcolumns',2,'interpreter','latex');
set(hlegend,'FontSize',18,'LineWidth',1.0,'location','northeast','orientation','horizontal','numcolumns',2);

set(hfig,'papersize',[12,4.5],'paperposition',[-0.7,0.1,13.5,4.5]);

%print('fig2.pdf');















