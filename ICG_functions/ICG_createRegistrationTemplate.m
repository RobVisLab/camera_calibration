function template = ICG_createRegistrationTemplate(template_type, num_template_rows, num_template_cols, verbose);
% template = ICG_createRegistrationTemplate(template_type, num_template_rows, num_template_cols, verbose);

if mod(num_template_rows,8)~=0 || mod(num_template_cols,8)~=0,
    error('ICG_createRegistrationTemplate: num_template_rows and num_template_cols must be divisible by 8.');
end;

num_template_rows_eighth = num_template_rows/8;
num_template_cols_eighth = num_template_cols/8;
center_xy = [];
template = [];
if strcmp(template_type, 'checkerboard'),
    template.pattern = zeros(num_template_rows,num_template_cols);

    template.pattern(num_template_rows_eighth+1:4*num_template_rows_eighth, :) = xor(template.pattern(num_template_rows_eighth+1:4*num_template_rows_eighth, :),1);
    template.pattern(7*num_template_rows_eighth+1:8*num_template_rows_eighth, :) = xor(template.pattern(7*num_template_rows_eighth+1:8*num_template_rows_eighth, :),1);

    template.pattern(:, num_template_cols_eighth+1:4*num_template_cols_eighth) = xor(template.pattern(:, num_template_cols_eighth+1:4*num_template_cols_eighth),1);
    template.pattern(:, 7*num_template_cols_eighth+1:8*num_template_cols_eighth) = xor(template.pattern(:, 7*num_template_cols_eighth+1:8*num_template_cols_eighth),1);
    template.center_xy = flipdim(size(template.pattern)/2+0.5, 2)';
    grid_distance_pix = 4*num_template_rows_eighth;
    template.neighbors = [0     0    0 -0.75 0.75;...
                          0 -0.75 0.75     0    0];
    template.feature_points = [-0.75/2  0.75/2  0.75/2 -0.75/2;... %Inside checker to check color
                               -0.75/2 -0.75/2  0.75/2  0.75/2];
elseif strcmp(template_type, 'circles'),
    template.pattern = ones(num_template_rows,num_template_cols);
    disk = fspecial('disk', num_template_rows_eighth);
    disk = disk/max(disk(:));
    disk = abs(disk-1);
    template.pattern(num_template_rows_eighth:num_template_rows_eighth+size(disk,1)-1, ...
                     num_template_cols_eighth:num_template_cols_eighth+size(disk,2)-1) = disk;
    template.pattern(end-num_template_rows_eighth+1-size(disk,1):end-num_template_rows_eighth, ...
                     num_template_cols_eighth:num_template_cols_eighth+size(disk,2)-1) = disk;
    template.pattern(num_template_rows_eighth:num_template_rows_eighth+size(disk,1)-1, ...
                     end-num_template_cols_eighth+1-size(disk,1):end-num_template_cols_eighth) = disk;
    template.pattern(end-num_template_rows_eighth+1-size(disk,1):end-num_template_rows_eighth, ...
                     end-num_template_cols_eighth+1-size(disk,1):end-num_template_cols_eighth) = disk;
    template.center_xy = flipdim(size(template.pattern)/2, 2)';
    grid_distance_pix = 4*num_template_rows_eighth;
    template.neighbors = [0  0 0 -1 1;...
                          0 -1 1  0 0];
    template.feature_points = [-0.5  0.5  0.5 -0.5;... %Inside circles
                               -0.5 -0.5  0.5  0.5];
elseif strcmp(template_type, 'circles_gray'),
    template.pattern = ones(num_template_rows,num_template_cols);
    disk = fspecial('disk', num_template_rows_eighth);
    disk = disk/max(disk(:));
    disk = abs(disk-1);
    dist_disk = bwdist(logical(disk));
    dist_disk = double(dist_disk/max(dist_disk(:)));
    disk = disk+dist_disk;
    template.pattern(num_template_rows_eighth:num_template_rows_eighth+size(disk,1)-1, ...
                     num_template_cols_eighth:num_template_cols_eighth+size(disk,2)-1) = disk;
    template.pattern(end-num_template_rows_eighth+1-size(disk,1):end-num_template_rows_eighth, ...
                     num_template_cols_eighth:num_template_cols_eighth+size(disk,2)-1) = disk;
    template.pattern(num_template_rows_eighth:num_template_rows_eighth+size(disk,1)-1, ...
                     end-num_template_cols_eighth+1-size(disk,1):end-num_template_cols_eighth) = disk;
    template.pattern(end-num_template_rows_eighth+1-size(disk,1):end-num_template_rows_eighth, ...
                     end-num_template_cols_eighth+1-size(disk,1):end-num_template_cols_eighth) = disk;
    template.center_xy = flipdim(size(template.pattern)/2, 2)';
    grid_distance_pix = 4*num_template_rows_eighth;
    template.neighbors = [0  0 0 -1 1;...
                          0 -1 1  0 0];
    template.feature_points = [-0.5  0.5  0.5 -0.5;... %Inside circles
                               -0.5 -0.5  0.5  0.5];
end;

% if strcmp(template_type, 'circles'),
%
%     corners_pix = [4*num_template_cols_eighth+.5  4*num_template_cols_eighth+.5 4*num_template_cols_eighth+.5 0*num_template_cols_eighth+.5 8*num_template_cols_eighth+.5 ;...
%                    4*num_template_rows_eighth+.5  0*num_template_rows_eighth+.5 8*num_template_rows_eighth+.5 4*num_template_rows_eighth+.5 4*num_template_rows_eighth+.5];
%     template = template(1:end-1, 1:end-1);
% elseif strcmp(template_type, 'checkerboard'),
%     corners_pix = [4*num_template_cols_eighth+.5  4*num_template_cols_eighth+.5 4*num_template_cols_eighth+.5 1*num_template_cols_eighth+.5 7*num_template_cols_eighth+.5 ;...
%                    4*num_template_rows_eighth+.5  1*num_template_rows_eighth+.5 7*num_template_rows_eighth+.5 4*num_template_rows_eighth+.5 4*num_template_rows_eighth+.5];
% end;

template.feature_points = template.feature_points * grid_distance_pix;
template.feature_points(1,:) = template.feature_points(1,:) + template.center_xy(1);
template.feature_points(2,:) = template.feature_points(2,:) + template.center_xy(2);
template.feature_points(3,:) = 1;

template.neighbors = template.neighbors * grid_distance_pix;
template.neighbors(1,:) = template.neighbors(1,:) + template.center_xy(1);
template.neighbors(2,:) = template.neighbors(2,:) + template.center_xy(2);
template.neighbors(3,:) = 1;

template.center_xy(3)=1;

if exist('verbose', 'var') && verbose >= 2,
    figure, imshow(template.pattern);
    hold on, plot(template.neighbors(1,:), template.neighbors(2,:), 'b.');
    hold on, plot(template.center_xy(1), template.center_xy(2), 'r.');
    hold on, plot(template.feature_points(1,:), template.feature_points(2,:), 'y.');
end;