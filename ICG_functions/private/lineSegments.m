function [B,corner_list_3] = lineSegments(B,area,corner_thresh)
global corner_thresh_2;
corner_thresh_2 = corner_thresh;

corner_list_3 = [];

region_index = 1;
while region_index <= length(B)
    chain = B{region_index}';
    num_chain_points = size(chain,2);

    %% find a corner in the chain to start with
    start_point = chain(:,1);
    end_points  = chain(:,1:end-1);
    distances = sqrt(sum((end_points - repmat(start_point, [1, size(end_points,2)])).^2));
    %% the point with max distance to any points on the chain is always a
    %% corner
    [max_distance, max_distance_index] = max(distances);

    % rotate the chain around to have the corner on pos 1, end pos is per
    % convention = start pos
    chain_1_mm1 = chain(:,1:max_distance_index-1);
    chain_m_em1 = chain(:,max_distance_index:end-1);
    chain = [chain_m_em1 chain_1_mm1 chain(:,max_distance_index)];

    %% find the opposite corner
    start_point = chain(:,1);
    end_points  = chain(:,1:end-1);
    distances = sqrt(sum((end_points - repmat(start_point, [1, num_chain_points-1])).^2, 1));
    [max_distance, max_distance_index] = max(distances);

    [r corners] = get_corners(chain(1,:), chain(2,:),max_distance_index,area(region_index),num_chain_points);
    if (length(corners)~=4)
%         disp (sprintf('Object %d - not a square marker = found %d corners',region_index,length(corners)));
        B(region_index) = [];
        % use the same region_index again, as a different chain is at
        % that index now :-)
    else
%         disp (sprintf('Object %d - Possibly a marker - successfully found 4 corners',region_index));

        B{region_index} = chain';
        corner_list_3 = [corner_list_3; corners'];
        region_index = region_index + 1;
    end;
end;
% if (isempty(corner_list_3))
%     disp('#######################');
%     disp('MATLAB ARToolkit Error!');
%     disp('#######################');
%     disp('');
%     disp('Could not find an object with 4 corners');
%     disp('Reasons may include:');
%     disp('  1.    No marker is present');
%     disp('  2.    Luminance threshold set inappropriately');
%     disp('  3.    Corner Detection Threshold set inappropriately');
%     disp(' ');
% end;


function [r,corners] = get_corners(chain_col, chain_row, max_distance_index, area, chain_length)
global corner_thresh_2;

thresh = corner_thresh_2;
%thresh = (area/0.75)*0.01*0.1;
corners = zeros(2,1);
corners(1) = 1;

corners1 = zeros(1,1);
num_corners1 = 0;
%disp('First half search');
[r corners1 num_corners1] = get_corner(chain_col,chain_row,1,max_distance_index,thresh,corners1,num_corners1);
if ( r<0)
 %   disp('Failed to find corner');
end;

corners2 = zeros(1,1);
num_corners2 = 0;
%disp('Second half search');
[r corners2 num_corners2] = get_corner(chain_col,chain_row,max_distance_index,chain_length,thresh,corners2,num_corners2);
if (r <0)
    %disp('Failed to find corner');
end;

if (( num_corners1==1) && (num_corners2==1))
    %disp('Case 1:  Found one corner in both halves');
    corners(2) = corners1(1);
    corners(3) = max_distance_index;
    corners(4) = corners2(1);
elseif((num_corners1>1) && (num_corners2==0))
    %disp('Case 1:  Found more than one corner in first half, but none in second');
    half_index = round(max_distance_index/2);
    num_corners1 = 0;
    [r corners1 num_corners1] = get_corner(chain_col,chain_row,1,half_index,thresh,corners1,num_corners1);
    if ( r<0)
     %   disp('Failed to find corner');
    end;
    [r corners2 num_corners2] = get_corner(chain_col,chain_row,half_index,max_distance_index,thresh,corners1,num_corners1);
    if ( r<0)
     %   disp('Failed to find corner');
    end;
    if ((num_corners1==1)&&(num_corners2==1))
        corners(2)=  corners1(1);
        corners(3) = corners2(1);
        corners(4) = max_distance_index;
    else
     %   disp('Too many corners');
    end;
elseif((num_corners1==0)&&(num_corners2>1))
    %disp('Case 1:  Found more than one corner in second half, but none in first');
    half_index = round((max_distance_index+chain_length)/2);
    num_corners1 = 0;

    [r corners1 num_corners1] = get_corner(chain_col,chain_row,max_distance_index,half_index,thresh,corners1,num_corners1);
    if ( r<0)
      %  disp('Failed to find corner');
    end;

    [r corners2 num_corners2] = get_corner(chain_col,chain_row,half_index,chain_length,thresh,corners1,num_corners1);
    if ( r<0)
      %  disp('Failed to find corner');
    end;
    if ((num_corners1==1)&&(num_corners2==1))
        corners(2)=  max_distance_index;
        corners(3) = corners1(1);
        corners(4) = corners2(1);
    else
     %   disp('Too many corners');
    end;
end;
if (r>=0)
 %   corners
else
    %disp('Something wrong with corner finding...');
end;


function [r, corners, num_corners] = get_corner(chain_col, chain_row, start_index, end_index, thresh, corners, num_corners)
%disp(sprintf('Entering Routine: Searching Between %d and %d',start_index,end_index));
%corners
%num_corners

% Setup line between start- and endpoint
a = chain_row(end_index) - chain_row(start_index);
b = chain_col(start_index) - chain_col(end_index);
c = chain_col(end_index)*chain_row(start_index) -chain_row(end_index)*chain_col(start_index);

distances = a*chain_col(start_index+1 : end_index-1) + b*chain_row(start_index+1 : end_index-1) + c;
[max_distance, max_distance_index] = max(distances.^2);
max_distance_index = max_distance_index + start_index - 1;
% for (i = start_index+1:end_index-1)
%     dist = a*chain_col(i)+b*chain_row(i)+c;
%     if ((dist*dist)>max_distance)
%         max_distance=dist*dist;
%         max_distance_index = i;
%     end;
% end;


if ((max_distance/(a*a+b*b))>thresh)
  %  disp(sprintf ('Distance, %f greater than thresh, %f',max_distance/(a*a+b*b),thresh));
    [r corners num_corners] = get_corner(chain_col,chain_row,start_index,max_distance_index,thresh,corners,num_corners);
    if(r<0)
   %     disp(sprintf('Did not find 2nd Exiting'));
        r=-1;
        return;
    end;
    if (num_corners>5)
        r=-1;
        return;
    end;
    num_corners = num_corners+1;
    corners(num_corners) = max_distance_index;
   % disp(sprintf('Adding Vertex num_corners %d,  position %d',max_distance_index,num_corners));

    [r corners num_corners] = get_corner(chain_col,chain_row,max_distance_index,end_index,thresh,corners,num_corners);
    if (r<0)
        r=-1;return;
    end;
else
   % disp(sprintf ('Distance, %f not greater than thresh, %f',max_distance/(a*a+b*b),thresh));
end
%disp('Exiting from routine');
%corners
%num_corners
r=0;
return;
