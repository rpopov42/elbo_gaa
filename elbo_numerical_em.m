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

function [uq,vq] = elbo_numerical_em(u,logrspu,uq,vq,TN)

du = u(2) - u(1); bufuq = uq*ones(1,TN); bufvq = vq*ones(1,TN);

for tn = 2:TN
    dequ = sum(exp(-1/2*(u - uq).^2/vq).*(u - uq).*logrspu);

    d2equ = sum(exp(-1/2*(u - uq).^2/vq).*((u - uq).^2/vq - 1).*logrspu);

    uq = uq + sign(dequ)*abs(dequ)/max(abs(d2equ),1/sqrt(vq)*abs(dequ));

    deqv = sum(exp(-1/2*(u - uq).^2/vq).*((u - uq).^2/vq - 1).*logrspu)*du + sqrt(2*vq*pi);

    d2eqv = 1/2*sum(exp(-1/2*(u - uq).^2/vq).*(((u - uq).^2/vq).^2 - 6*(u - uq).^2/vq + 3).*logrspu)*du - sqrt(2*vq*pi);

    vq = vq*(1 + sign(deqv)*abs(deqv)/max(abs(d2eqv),2*abs(deqv)));

    bufuq(tn) = uq; bufvq(tn) = vq;
end;

uq = bufuq; vq = bufvq;

