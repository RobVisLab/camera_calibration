function bw_img = adaptivethreshold(img, min_area)

    img = im2double(img); % just in case
    img = mat2gray(img);
    scale_factor = 0.3;
    sub_img = imresize(img, scale_factor);

    % Estimate the necessary window size (should be 4 checkers wide).
    % The min_area is the area of a square target which is 3 checkers wide.
    window_size = round(sqrt(min_area) * 4/3 * scale_factor);

    m_sub_img = imfilter(sub_img, fspecial('average', window_size), 'replicate');
    m_img = imresize(m_sub_img, size(img));
    bw_img = (img - m_img)>0;
