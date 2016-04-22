function template = make_target (grid_width_pixels, grid_width_mm, grid_coordinates_h, grid_coordinates_v, verbose)
% template = make_target (grid_width_pixels, grid_width_mm, grid_coordinates_h, grid_coordinates_v, verbose)
addpath ICG_functions/
if mod(grid_width_pixels, 4) ~= 0, % mod(grid_width_pixels, 16)
    warning('make_target: grid_width_pixels should be divisible by 4.');
end;

% grid_width_mm = 5;
% grid_coordinates_h = -20:20;
% grid_coordinates_v = -10:8;
% % % grid_coordinates_h = -5:5;
% % grid_coordinates_v = -18:16;
% % grid_coordinates_h = [-5:5];
% % grid_coordinates_v = [-7:7];
grid_coordinates_h_mm = grid_coordinates_h*grid_width_mm;
grid_coordinates_v_mm = grid_coordinates_v*grid_width_mm;

% elementary_template_size_mm = 2*grid_width_mm;
elementary_template_size_pix = grid_width_pixels*2;

registration_template = ICG_createRegistrationTemplate('circles', elementary_template_size_pix, elementary_template_size_pix, verbose);

%template.pattern = repmat(registration_template.pattern, [length(grid_coordinates_v)+1, length(grid_coordinates_h)+1]/2);
template.pattern = repmat(registration_template.pattern(1:end/2,1:end/2), [length(grid_coordinates_v)+1, length(grid_coordinates_h)+1]);

grid_coordinates_h_pix = (grid_coordinates_h-grid_coordinates_h(1))    * elementary_template_size_pix/2 + registration_template.center_xy(1);
grid_coordinates_v_pix = (-grid_coordinates_v+grid_coordinates_v(end)) * elementary_template_size_pix/2 + registration_template.center_xy(2);

[cx, cy] = meshgrid(grid_coordinates_h_pix, grid_coordinates_v_pix);
[gx, gy] = meshgrid(grid_coordinates_h_mm, grid_coordinates_v_mm);
template.corner_points = [cx(:) cy(:)]';
template.grid_points   = [gx(:) gy(:)]';

translation_grid_coord_circle_center_ul = registration_template.feature_points(:,1) - registration_template.center_xy;
circle_centers_h_pix = [grid_coordinates_h_pix grid_coordinates_h_pix(end)+diff(grid_coordinates_h_pix(1:2))]+translation_grid_coord_circle_center_ul(1);
circle_centers_v_pix = [grid_coordinates_v_pix(1)-diff(grid_coordinates_v_pix(1:2)) grid_coordinates_v_pix ]+translation_grid_coord_circle_center_ul(2);
%circle_centers_h_mm = [grid_coordinates_h_mm grid_coordinates_h_mm(end)+2]-1;
%circle_centers_v_mm = [grid_coordinates_v_mm(1)-2 grid_coordinates_v_mm]+1;
circle_centers_h_mm = [grid_coordinates_h_mm grid_coordinates_h_mm(end)+grid_width_mm]-grid_width_mm/2;
circle_centers_v_mm = [grid_coordinates_v_mm(1)-grid_width_mm grid_coordinates_v_mm]+grid_width_mm/2;
[cx, cy] = meshgrid(circle_centers_h_pix, circle_centers_v_pix);
[gx, gy] = meshgrid(circle_centers_h_mm, circle_centers_v_mm);

template.circle_centers_pix = [cx(:) cy(:)]';
template.circle_centers_mm  = [gx(:) gy(:)]';

template.marker = ICG_createCalibrationMarker(grid_width_pixels*3.5, 7, 'cirlces', verbose);
%left_marker_circle_index = find(template.circle_centers_mm(1,:)==-1 & template.circle_centers_mm(2,:)==-grid_width_mm+1); % template.circle_centers_mm(2,:)==-4
left_marker_circle_index = find(template.circle_centers_mm(1,:)<0 & template.circle_centers_mm(1,:)>=-grid_width_mm & template.circle_centers_mm(2,:)<-grid_width_mm & template.circle_centers_mm(2,:)>=-grid_width_mm*2); % template.circle_centers_mm(2,:)==-4
left_marker_circle_pix = template.circle_centers_pix(:, left_marker_circle_index);
left_marker_circle_mm = template.circle_centers_mm(:, left_marker_circle_index);
template.marker_offset_pix = left_marker_circle_pix - template.marker.circle_center_corners(:,1);
%template.marker_offset_mm = left_marker_circle_mm - template.marker.circle_center_gridpoints(:,1); % DONE: fix position!
template.marker_offset_mm = left_marker_circle_mm - ([1;-1] .* template.marker.circle_center_corners(:,1)*grid_width_mm/grid_width_pixels);
template.pattern(template.marker_offset_pix(2)+1:template.marker_offset_pix(2)+size(template.marker.pattern, 1),...
                 template.marker_offset_pix(1)+1:template.marker_offset_pix(1)+size(template.marker.pattern, 2)) = template.marker.pattern;

% Padding
template.pattern = padarray(template.pattern, [elementary_template_size_pix elementary_template_size_pix], 1, 'both');
template.corner_points = template.corner_points + elementary_template_size_pix;
template.circle_centers_pix = template.circle_centers_pix + elementary_template_size_pix;
template.marker_offset_pix = template.marker_offset_pix + elementary_template_size_pix;

% Visualization
if exist('verbose', 'var') && verbose >= 2,
    figure, imshow (template.pattern);
    
    corner_points = template.corner_points;
    bad = template.grid_points(1,:) >= -2 * grid_width_mm & template.grid_points(1,:) <= 2 * grid_width_mm & ...
          template.grid_points(2,:) >= -4 * grid_width_mm & template.grid_points(2,:) <= 0 * grid_width_mm;
    corner_points(:, bad) = [];
    circle_centers_pix = template.circle_centers_pix;
    bad = template.circle_centers_mm(1,:) >= -2 * grid_width_mm & template.circle_centers_mm(1,:) <= 2 * grid_width_mm & ...
          template.circle_centers_mm(2,:) >= -4 * grid_width_mm & template.circle_centers_mm(2,:) <= 0 * grid_width_mm;
    circle_centers_pix(:, bad) = [];
    
    hold on, plot(corner_points(1,:), corner_points(2,:), 'b.');
    hold on, plot(circle_centers_pix(1,:), circle_centers_pix(2,:), 'r.');
    
    mp = template.marker.corner_points;
    mp(1,:) = mp(1,:) + template.marker_offset_pix(1);
    mp(2,:) = mp(2,:) + template.marker_offset_pix(2);
    dmp = diff(mp, 1, 2);
    hold on, quiver(mp(1,1:end-1), mp(2,1:end-1), dmp(1,:), dmp(2,:), 0, 'Color', 'y');
    
    cc = template.marker.circle_center_corners;
    cc(1,:) = cc(1,:) + template.marker_offset_pix(1);
    cc(2,:) = cc(2,:) + template.marker_offset_pix(2);
    hold on, plot(cc(1,:), cc(2,:), 'g.');
    
    coord_center = template.corner_points(:, template.grid_points(1,:) == 0 & template.grid_points(2,:) == 0);
    coord_up = template.corner_points(:, template.grid_points(1,:) == 0 & template.grid_points(2,:) == 2*grid_width_mm); % template.grid_points(2,:) == 10
    coord_right = template.corner_points(:, template.grid_points(1,:) == 2*grid_width_mm & template.grid_points(2,:) == 0); % template.grid_points(1,:) == 10
    hold on, quiver(coord_center(1), coord_center(2), coord_up(1)-coord_center(1), coord_up(2)-coord_center(2), 0, 'LineWidth', 2, 'Color', 'm');
    hold on, quiver(coord_center(1), coord_center(2), coord_right(1)-coord_center(1), coord_right(2)-coord_center(2), 0, 'LineWidth', 2, 'Color', 'g');
    
    if verbose >= 3,
        for i=1:size(template.grid_points , 2),
            hold on, text(template.corner_points(1,i)+2, template.corner_points(2,i)+2, ...
                sprintf('(%.2f,%.2f)', template.grid_points(1,i), template.grid_points(2,i)), 'FontSize', 6);
        end;
    end;
end;
