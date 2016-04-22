function matches = GetOrientationIdentity (thresh_image, template, corners_col, corners_row)
%cross correlates with stored version

[num_squares dummy]=size(corners_col);
matches = cell(1,num_squares);

for square_index = 1:num_squares
    square_corners_col = corners_col(square_index,:);
    square_corners_row = corners_row(square_index,:);
    source_points = [square_corners_col', square_corners_row'];
    target_points = [1 32; 1 1; 32 1; 32 32];
    try %%% the following might fail in case the source points are degenerate!
        projective_transformation = maketform('projective', source_points, target_points);
        rectified_pattern = imtransform(thresh_image, projective_transformation, 'bilinear', 'XData', [1, 32], 'YData', [1, 32]);

        corr = zeros(1,4);
        for temp_orientation = 1:4,

            %%% calculating the whole corralation is a waste! only one value is
            %%% needed
            %temp_correlation = normxcorr2 (template, rectified_pattern);              
            %corr(temp_orientation) = temp_correlation(32,32);
        
            corr(temp_orientation) = getSimilarity(template,rectified_pattern);
            template = rot90(template);
        end;
        [max_corr, max_ind] = max(corr);
        corner_sequence = circshift(1:4, [0 max_ind-1]);
        matches{square_index} = struct('similarity', max_corr, 'pattern_corners',...
                                       [corners_col(square_index,corner_sequence); corners_row(square_index,corner_sequence)]');
    catch
        matches{square_index} = [];
    end
end
