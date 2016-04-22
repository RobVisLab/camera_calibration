function marker = ICG_createCalibrationMarker(black_width_pix, black_width_mm, type, verbose);
% marker = ICG_createCalibrationMarker(black_width_pix, black_width_mm, type, verbose);

white_border_width_pix = black_width_pix / 7 -1;
black_border_width_pix = black_width_pix / 7;
pattern_width_pix = black_width_pix * 5/7;

%dot
dot_radius_pix = black_width_pix/14;
dot = fspecial('disk', dot_radius_pix); % always odd number of rows/cols
dot = dot/max(dot(:));
dot = abs(dot-1);
dot_pattern_width = size(dot, 2);
dot_center_pix = [dot_pattern_width/2; dot_pattern_width/2]+0.5;

padsize_diameter =dot_pattern_width-1;
padsize_radius = padsize_diameter/2;
dot = padarray(dot, [padsize_radius padsize_radius], 1, 'both');
dot_center_pix = dot_center_pix+padsize_radius;

dot = dot(1:end-1, 1:end-1);
pattern = repmat(dot, [2 2]);
%     whitedot=ones(size(dot));
%     pattern = repmat([dot whitedot], [2 1]);
pattern = padarray(pattern, [1 1], 1, 'post');
dot_center_pix = [dot_center_pix dot_center_pix];
dot_center_pix(1,2) = dot_center_pix(1,1)+size(dot,2);

pattern = padarray(pattern, [padsize_radius padsize_radius]-1, 1, 'both');
dot_center_pix = dot_center_pix+padsize_radius-1;

pattern(floor(end/2):end, :) = 0;
pattern = padarray(pattern, [1 1], 0.5, 'both');
dot_center_pix = dot_center_pix+1;
pattern = padarray(pattern, [black_border_width_pix-1 black_border_width_pix-1], 0, 'both');
dot_center_pix = dot_center_pix+black_border_width_pix-1;

pattern = padarray(pattern, [1 1], 0.5, 'both');
dot_center_pix = dot_center_pix+1;
marker.pattern = padarray(pattern, [white_border_width_pix white_border_width_pix], 1, 'both');
dot_center_pix = dot_center_pix+white_border_width_pix;

ul = [white_border_width_pix white_border_width_pix]' + 1;
marker.corner_points = [ul(1) ul(1)+black_width_pix ul(1)+black_width_pix ul(1);...
                        ul(2) ul(2)                   ul(2)+black_width_pix ul(2)+black_width_pix];
marker.grid_points = [0 black_width_mm black_width_mm 0;...
                      0 0              black_width_mm black_width_mm];
marker.type = type;
marker.circle_center_corners = dot_center_pix;
marker.circle_center_gridpoints(1,:) = marker.circle_center_corners(1,:) - marker.corner_points(1,1);
marker.circle_center_gridpoints(2,:) = marker.circle_center_corners(2,:) - marker.corner_points(2,1);
marker.circle_center_gridpoints = marker.circle_center_gridpoints / black_width_pix * black_width_mm;

if exist('verbose', 'var') && verbose >= 2,
    figure, imshow(marker.pattern);
    hold on, plot(marker.corner_points(1,:), marker.corner_points(2,:), 'b.');
    hold on, plot(marker.circle_center_corners(1,:), marker.circle_center_corners(2,:), 'r.');
    if verbose >= 3,
        for i=1:size(marker.circle_center_corners, 2),
            hold on, text(marker.circle_center_corners(1,i)+2, marker.circle_center_corners(2,i)+2, ...
                sprintf('(%.2f, %.2f)', marker.circle_center_gridpoints(1,i), marker.circle_center_gridpoints(2,i)), ...
                'FontSize', black_width_mm, 'Color', 'r');
        end;
        for i=1:size(marker.corner_points, 2),
            hold on, text(marker.corner_points(1,i)+2, marker.corner_points(2,i)+2, ...
                sprintf('(%.2f, %.2f)', marker.grid_points(1,i), marker.grid_points(2,i)), ...
                'FontSize', black_width_mm, 'Color', 'b');
        end;
    end;
end;

