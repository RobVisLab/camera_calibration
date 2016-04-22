function [corner_points_col, corner_points_row]= calcIntersection(B,corners)

skip_value = 5;

corner_points_col = [];
corner_points_row = [];

num_regions = length(B);
for region_index = 1:num_regions
    num_chain_points = size(B{region_index},1);

    line_start = [];
    line_direction = [];
    for corner_index = 1:4        %find vector to and vector along lines        
        first_point_on_line = corners(region_index,corner_index)+skip_value;
        last_point_on_line = corners(region_index,rem(corner_index,4)+1)-skip_value;

        if (first_point_on_line < 1)
            first_point_on_line = first_point_on_line+num_chain_points;
        end;
        if (first_point_on_line > num_chain_points)
            first_point_on_line = first_point_on_line-num_chain_points;
        end;
        if (last_point_on_line < 1)
            last_point_on_line = last_point_on_line+num_chain_points;
        end;
        if (first_point_on_line > last_point_on_line)
            temp = first_point_on_line;
            first_point_on_line=last_point_on_line;
            last_point_on_line = temp;
        end;

        start_point = [B{region_index}(first_point_on_line,2); B{region_index}(first_point_on_line,1)];
        line_start= [line_start start_point];
     
        points_to_fit_line = [B{region_index}(first_point_on_line:last_point_on_line,2)'; ...
                              B{region_index}(first_point_on_line:last_point_on_line,1)'];
        points_to_fit_line(1,:) = points_to_fit_line(1,:)-mean(points_to_fit_line(1,:),2);
        points_to_fit_line(2,:) = points_to_fit_line(2,:)-mean(points_to_fit_line(2,:),2);
        points_to_fit_line = points_to_fit_line*points_to_fit_line';
        [U L V] = svd(points_to_fit_line);
        line_direction = [line_direction U(:,1)];                            
    end;

    intersection_col = [];
    intersection_row = [];

    for corner_index = 1:4  %find co-ordinates of where lines meet by solve equation Ml = b
        line1 = corner_index;
        line2 = rem(corner_index,4)+1;
        b = line_start(:,line2)-line_start(:,line1);
        M = [line_direction(:,line1) -1*line_direction(:,line2)];
        l = M\b;
        intersection_col =[intersection_col line_start(1,line1)+l(1,1)*line_direction(1,line1)]; 
        intersection_row =[intersection_row line_start(2,line1)+l(1,1)*line_direction(2,line1)];           
    end;
    % correct point sequence
    intersection_points = [intersection_col; intersection_row];
    centroid = mean(intersection_points, 2);
    dirs = intersection_points - repmat(centroid, [1, 4]);
    for i=1:size(dirs, 2),
        dirs(:,i) = dirs(:,i) / norm(dirs(:,i));
    end;
    start_index = max([find(dirs(1,:)<0 & dirs(2,:)>0), 1]);
    intersection_points = circshift (intersection_points, [0 1-start_index]);
    
    
    corner_points_col = [corner_points_col; intersection_points(1,:)];
    corner_points_row = [corner_points_row; intersection_points(2,:)];
end;
