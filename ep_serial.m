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

function [uq,vq] = ep_serial(x,vg,up,vp,w,Pc,TN)

uq = up; vq = vp; uc = zeros(length(x),1); auc = zeros(length(x),1);

bufuq = up*ones(1,TN); bufvq = vp*ones(1,TN);

for t = 1:TN
    for n = 1:length(x)
        bvuc = 1/(1/vq - auc(n)); buc = uq + (uq - uc(n))*bvuc*auc(n);

        bu = (buc*vg + x(n)*bvuc)/(vg + bvuc); bv = 1/(1/vg + 1/bvuc); bmag = 1/sqrt(2*(vg + bvuc)*pi)*exp(-1/2*(x(n) - buc)^2/(vg + bvuc));

        p1 = (1 - w)*bmag/((1 - w)*bmag + w*Pc(n)); p2 = w*Pc(n)/((1 - w)*bmag + w*Pc(n));

        uq = p1*bu + p2*buc; vq = p1*(bv + (bu - uq)^2) + p2*(bvuc + (buc - uq)^2);

        auc(n) = ((1/vq >= 1/bvuc) - (1/vq < 1/bvuc))*max(abs(1/vq - 1/bvuc),0.1^10); uc(n) = uq + (uq - buc)/auc(n)/bvuc;
    end;

    if iscomplex(uq) | isnan(uq) | iscomplex(vq) | isnan(vq), t = max(t-1,1); break; end;

    bufuq(t) = uq; bufvq(t) = vq;
end;

uq = bufuq(1:t); vq = bufvq(1:t);

