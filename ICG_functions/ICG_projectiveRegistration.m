function fit = ICG_projectiveRegistration (img, tmplt, p_init, n_iters, verbose);
% HOMO_IC - Homography image alignment using inverse-compositional algorithm
%   FIT = HOMO_IC(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using a
%   projective warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p4, p7     paper_equivalent = [p1, p3, p5
%             p2, p5, p8                         p2, p4, p6
%             p3, p6, 1];                        p7, p8, 1];
%
%   This assumes greyscale images and rectangular templates.
%
%   c.f. Baker-Matthews
%
% ATTENTION: this is NOT a classical homography!!!!!!
%            If H(:) = [h1, h2, ..., h8, 1], is a homography, then 
%               p_init(i) = h(i) for i elem {2, 3, 4, 6, 7, 8} and h(i)-1
%               for i elem {1, 5}.

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: ICG_projectiveRegistration.m,v 1.2 2006-12-18 13:48:54 ruether Exp $

fit = homo_ic (img, tmplt, p_init, n_iters, verbose);