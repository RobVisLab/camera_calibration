function intervals_ndimsx2 = ICG_boundingBox (points);
% intervals_ndimsx2 = ICG_boundingBox (points);

ndims = size(points, 1);
intervals_ndimsx2 = zeros(ndims, 2);
for i=1:ndims,
    intervals_ndimsx2(i,1) = min(points(i,:));
    intervals_ndimsx2(i,2) = max(points(i,:));
end;