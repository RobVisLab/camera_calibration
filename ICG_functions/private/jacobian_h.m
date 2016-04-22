function dW_dp = jacobian_h(nx, ny, warp_p);
% JACOBIAN_H - Compute Jacobian for projective warp
%   DW_DP = JACOBIAN_H(WIDTH, HEIGHT, P)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: jacobian_h.m,v 1.2 2007-06-08 14:07:18 ruether Exp $

% Easy bits
jac_x = kron([0:nx - 1],ones(ny, 1));
jac_y = kron([0:ny - 1]',ones(1, nx));
jac_zero = zeros(ny, nx);
jac_one = ones(ny, nx);

% Complicated bits are just homography of all image coordinates
xy = [repmat([1:ny]',nx,1) kron([1:nx]', ones(ny,1))];
xy = [xy ones(length(xy),1)]';
M = warp_p;
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;
uv = M * xy;
% uvc = uv ./ repmat(uv(3,:),3,1);
uvc = uv;
uvc(1,:) = uvc(1,:) ./ uv(3,:);
uvc(2,:) = uvc(2,:) ./ uv(3,:);
uvc(3,:) = uvc(3,:) ./ uv(3,:);

u_x = reshape(uvc(1,:),nx,ny)';
u_y = reshape(uvc(2,:),nx,ny)';
v = reshape(uv(3,:),nx,ny)';

% Divide each jacobian image by v
iv = 1 ./ v;
jac_x = iv .* jac_x;
jac_y = iv .* jac_y;
jac_one = iv .* jac_one;

[nr, nc] = size(jac_x);
dW_dp = zeros(2*nr, 8*nc);
dW_dp(1:nr, 1:nc) = jac_x;
dW_dp(1:nr, nc+1:2*nc) = jac_zero;
dW_dp(1:nr, 2*nc+1:3*nc) = -jac_x .* u_x;
dW_dp(1:nr, 3*nc+1:4*nc) = jac_y;
dW_dp(1:nr, 4*nc+1:5*nc) = jac_zero;
dW_dp(1:nr, 5*nc+1:6*nc) = -jac_y .* u_x;
dW_dp(1:nr, 6*nc+1:7*nc) = jac_one;
dW_dp(1:nr, 7*nc+1:8*nc) = jac_zero;

dW_dp(nr+1:2*nr, 1:nc) = jac_zero;
dW_dp(nr+1:2*nr, nc+1:2*nc) = jac_x;
dW_dp(nr+1:2*nr, 2*nc+1:3*nc) = -jac_x .* u_y;
dW_dp(nr+1:2*nr, 3*nc+1:4*nc) = jac_zero;
dW_dp(nr+1:2*nr, 4*nc+1:5*nc) = jac_y;
dW_dp(nr+1:2*nr, 5*nc+1:6*nc) = -jac_y .* u_y;
dW_dp(nr+1:2*nr, 6*nc+1:7*nc) = jac_zero;
dW_dp(nr+1:2*nr, 7*nc+1:8*nc) = jac_one;




% dW_dp = [jac_x, jac_zero, -jac_x .* u_x, jac_y, jac_zero, -jac_y .* u_x, jac_one, jac_zero;
%         jac_zero, jac_x, -jac_x .* u_y, jac_zero, jac_y, -jac_y .* u_y, jac_zero, jac_one];