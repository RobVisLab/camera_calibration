function sd_delta_p = sd_update(VI_dW_dp, error_img, N_p, w)
% SD_UPDATE - Compute steepest descent parameter updates
%   SD_DELTA_P = SD_UPDATE(VI_DW_DP, ERROR_IMG, N_P, W)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: sd_update.m,v 1.2 2007-05-03 12:22:37 ruether Exp $

if nargin<4 error('Not enough input arguments'); end

sd_delta_p = zeros(N_p, 1);
for p=1:N_p
	h1 = VI_dW_dp(:,((p-1)*w)+1:((p-1)*w)+w);
    product = h1 .* error_img;
    sum1 = sum(product);
	sd_delta_p(p) = sum(sum1);
end
