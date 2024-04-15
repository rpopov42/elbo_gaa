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

function [uq,vq] = laplace(x,vg,up,vp,w,Pc,TN)

p = 0.5*ones(length(x),1); bufuq = zeros(1,TN); bufvq = zeros(1,TN);

for t = 1:TN
    uq = (sum(p.*x)/vg + up/vp)/(sum(p)/vg + 1/vp);

    rsp = 1/sqrt(2*vg*pi)*exp(-1/2*(x - uq).^2/vg);

    p = (1 - w)*rsp./((1 - w)*rsp + w*Pc);

    vq = 1/(max(sum(p.*(1 - (1 - p).*(x - uq).^2/vg)),0)/vg + 1/vp);

    bufuq(t) = uq; bufvq(t) = vq;
end;

uq = bufuq; vq = bufvq;

