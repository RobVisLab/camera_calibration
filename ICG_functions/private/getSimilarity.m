function similarity =  getSimilarity(template,rectified_pattern)
    template = double(template(:));
    template = template-mean(template);
    rectified_pattern = double(rectified_pattern(:));
    rectified_pattern = rectified_pattern-mean(rectified_pattern);
        
    similarity = sum(template.*rectified_pattern)/(norm(template)*norm(rectified_pattern));
