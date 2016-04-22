function [labeled_image,area] = ConnectedComponents(thresholded_image,min_region_area)
% [labeled_image,area] = ConnectedComponents(thresholded_image,min_region_area)
[num_rows num_cols] = size(thresholded_image);

    
% Add frame around image
thresholded_image(1,:)=1;
thresholded_image(num_rows,:)=1;
thresholded_image(:,1)=1;
thresholded_image(:,num_cols)=1;
inverted_image = ~thresholded_image;
% if(min(num_rows,num_cols)>200)
    inverted_image = imopen(inverted_image, strel('square', 3));
    inverted_image = imclose(inverted_image, strel('square', 3));
% else
%     inverted_image = imopen(inverted_image, strel('disk', 1));
%     inverted_image = imclose(inverted_image, strel('disk', 1));
% end
[temp_labeled_image num_regions] = bwlabel(inverted_image,8);
region_props = regionprops(temp_labeled_image, 'BoundingBox');
bb = cat(1, region_props(:).BoundingBox);
convex_areas = bb(:, 3) .* bb(:, 4);
good_label_mask = convex_areas >= min_region_area;
bad_label_mask = ~good_label_mask;
convex_areas(bad_label_mask) = [];
area = convex_areas;
labeled_image = zeros(size(temp_labeled_image));
good_labels = find(good_label_mask);
for i=1:length(good_labels),
    current_label = good_labels(i);
    mask = temp_labeled_image==current_label;
    labeled_image(mask) = current_label;
end;

% region_histogram = hist(reshape(thresholded_image,num_rows*num_cols,1),-0.5:1:num_regions+0.5);
% region_histogram = region_histogram(2:length(region_histogram));
% It's a tiny bit faster, when very large areas are dropped, too, but
% it's not as robust that way, especially when an unrealisticly small
% min_area is used.
% regions_to_keep = find(region_histogram>min_region_area); % & region_histogram<20*min_region_area);
% 
% labeled_image = zeros(num_rows,num_cols);
% final_region_index = 1;
% for current_region=regions_to_keep
%     [row col] = find(thresholded_image==current_region);
%     area(final_region_index) = length(row);
%     labeled_image((col-1)*num_rows+row) = final_region_index;
%     final_region_index = final_region_index+1;
% end
