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

eps = (-3:0.01:14.0); hvg = 5; w = 0.3; pci = 1/20;

hx1 = 2.0;

L1 = (1 - w)*1/sqrt(2*hvg*pi)*exp(-1/2*(eps - hx1).^2/hvg) + w*pci;

eps1 = 4.5; v1 = 5/6;

Lu1 = 1/sqrt(2*hvg*pi)*exp(-1/2*(eps1 - hx1).^2/hvg); p1 = (1 - w)*Lu1/((1 - w)*Lu1 + w*pci);

logLu1 = log((1 - w)*Lu1 + w*pci); dlogLu1 = -p1*(eps1 - hx1)/hvg; d2logLu1 = p1*((1 - p1)*(eps1 - hx1)^2/hvg - 1)/hvg;

tL1 = exp(1/2*d2logLu1*(eps - eps1).^2 + dlogLu1*(eps - eps1) + logLu1);

g1 = 0.3*1/sqrt(2*v1*pi)*exp(-1/2*(eps - eps1).^2/v1);

hx2 = 3.0;

L2 = (1 - w)*1/sqrt(2*hvg*pi)*exp(-1/2*(eps - hx2).^2/hvg) + w*pci;

eps2 = 8.5; v2 = 5/6;

Lu2 = 1/sqrt(2*hvg*pi)*exp(-1/2*(eps2 - hx2).^2/hvg); p2 = (1 - w)*Lu2/((1 - w)*Lu2 + w*pci);

logLu2 = log((1 - w)*Lu2 + w*pci); dlogLu2 = -p2*(eps2 - hx2)/hvg; d2logLu2 = p2*((1 - p2)*(eps2 - hx2)^2/hvg - 1)/hvg;

tL2 = exp(1/2*d2logLu2*(eps - eps2).^2 + dlogLu2*(eps - eps2) + logLu2);

g2 = 0.2*1/sqrt(2*v2*pi)*exp(-1/2*(eps - eps2).^2/v2);

hfig = figure;

subplot(1,2,1); hold on;
plot(eps,L1,'LineWidth',2.0,'color','blue');
plot(eps,tL1,'LineWidth',2.0,'color','red');
plot(eps,g1,'LineWidth',2.0,'color','green');
plot(eps1,exp(logLu1),'x','MarkerSize',10,'LineWidth',2.0,'color','black');

%xlim([eps(1),eps(end)]); xlabel ('$\boldsymbol{\epsilon}$','interpreter','latex');
xlim([eps(1),eps(end)]); xlabel ('\bf \epsilon');

%ylim([-0.0*max(L1),1.2*max(L1)]); ylabel ('\textbf{Likelihood}','interpreter','latex');
ylim([-0.0*max(L1),1.2*max(L1)]); ylabel ('\bf Likelihood');

set(gca,'FontSize',14,'LabelFontSizeMultiplier',1.5,'LineWidth',1.0);

%hlegend = legend('$\mathcal{L}(\epsilon|\hat x_1)$','$\mathcal{\widetilde L}_1(\epsilon)$',...
%                 '$\rho_1(\epsilon)$','$\mathcal{L}(\epsilon_1|\hat x_1)$');
hlegend = legend('L(\epsilon|\^x_1)','~L_1(\epsilon)','\rho_1(\epsilon)','L(\epsilon_1|\^x_1)');

%set(hlegend,'FontSize',22,'LineWidth',1.0,'interpreter','latex');
set(hlegend,'FontSize',22,'LineWidth',1.0);

subplot(1,2,2); hold on;
plot(eps,L2,'LineWidth',2.0,'color','blue');
plot(eps,tL2,'LineWidth',2.0,'color','red');
plot(eps,g2,'LineWidth',2.0,'color','green');
plot(eps2,exp(logLu2),'+','MarkerSize',12,'LineWidth',2.0,'color','black');

%xlim([eps(1),eps(end)]); xlabel ('$\boldsymbol{\epsilon}$','interpreter','latex');
xlim([eps(1),eps(end)]); xlabel ('\bf \epsilon');

%ylim([-0.0*max(L2),1.2*max(L2)]); ylabel ('\textbf{Likelihood}','interpreter','latex');
ylim([-0.0*max(L2),1.2*max(L2)]); ylabel ('\bf Likelihood');

set(gca,'FontSize',14,'LabelFontSizeMultiplier',1.5,'LineWidth',1.0);

%hlegend = legend('$\mathcal{L}(\epsilon|\hat x_2)$','$\mathcal{\widetilde L}_2(\epsilon)$',...
%                 '$\rho_2(\epsilon)$','$\mathcal{L}(\epsilon_2|\hat x_2)$');
hlegend = legend('L(\epsilon|\^x_2)','~L_2(\epsilon)','\rho_2(\epsilon)','L(\epsilon_2|\^x_2)');

%set(hlegend,'FontSize',22,'LineWidth',1.0,'interpreter','latex');
set(hlegend,'FontSize',22,'LineWidth',1.0);

set(hfig,'papersize',[12,4.5],'paperposition',[0.0,0.1,12.5,4.5]);

%print('fig1.pdf');




















