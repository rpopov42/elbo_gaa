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

function [uq,vq] = elbo_gaa_em_serial(x,vg,up,vp,w,Pc,uq,vq,TN)

hvg = max(2*vq,vg); bufuq = zeros(1,TN); bufvq = zeros(1,TN);

B = zeros(length(x),1); C = zeros(length(x),1); D = zeros(length(x),1);

for t = 1:TN
  for n = 1:length(x)
    rsp = 1/sqrt(2*hvg*pi)*exp(-1/2*hvg*(x(n) - uq)^2/(hvg + vq)^2);

    p = (1 - w)*rsp/((1 - w)*rsp + w*Pc(n));

    v = hvg/((1 - p)*(p*hvg*(x(n) - uq)^2/(hvg + vq)^2 + 1)*vq + hvg);

    A = exp(-1/2*vq*(1 - p^2*v)*(x(n) - uq)^2/(hvg + vq)^2);

    B(n) = p*sqrt(v)*A*(hvg + p*v*vq)/(hvg + vq);

    C(n) = p*sqrt(v)*A*v; D(n) = (1 - p*v)*B(n);
  end;

  uq = (sum(B.*x)/hvg + up/vp)/(sum(B)/hvg + 1/vp);

  vq = (sum(D.*(x - uq).^2)/hvg*vq/(hvg + vq) + 1)/(sum(C)/hvg + 1/vp);

  hvg = max(min(2*vq,hvg/2),vg); vq = min(vq,max(vg,hvg/2));

  bufuq(t) = uq; bufvq(t) = vq;
end;

uq = bufuq; vq = bufvq;


