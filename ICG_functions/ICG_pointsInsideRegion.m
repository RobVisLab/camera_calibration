function mask = ICG_pointsInsideRegion (points, intervals_ndimsx2);
% mask = ICG_pointsInsideRegion (points, intervals_ndimsx2);

mask = logical(ones(1, size(points, 2)));
ndims = size(points, 1);
for i=1:ndims,
    mask = mask & ...
        points(i,:) >= intervals_ndimsx2(i,1) & ...
        points(i,:) <= intervals_ndimsx2(i,2);    
end;

