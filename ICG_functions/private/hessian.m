function H = hessian(VI_dW_dp, N_p, w)
% HESSIAN - Compute Hessian
%   H = HESSIAN(VI_DW_DP, N_P, W)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: hessian.m,v 1.2 2007-05-03 12:16:51 ruether Exp $

if nargin<3 error('Not enough input arguments'); end

H = zeros(N_p, N_p);
for i=1:N_p
	h1 = VI_dW_dp(:,((i-1)*w)+1:((i-1)*w)+w);
	for j=1:N_p
		h2 = VI_dW_dp(:,((j-1)*w)+1:((j-1)*w)+w);
        product = h1 .* h2;
        sum1 = sum(product);
		H(j, i) = sum(sum1);
	end
end
