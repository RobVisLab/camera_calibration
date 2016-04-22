function [marker_corners, marker_similarity, maxIdx, maxTemplateIdx] = selectBestMatch(matches, isTemplateArray)
    
    if nargin < 2
        isTemplateArray = 0;
    end
    
    if isTemplateArray ~= 0
        marker_similarity = -Inf;
        marker_corners = [];
        maxIdx = 0;
        for tIdx = 1:length(matches)
            [mc, msim, midx] = select(matches{tIdx});
            if msim > marker_similarity
                marker_similarity = msim;
                marker_corners = mc;
                maxIdx = midx;
                maxTemplateIdx = tIdx;
            end
        end
    else
        maxTemplateIdx = 1;
        [marker_corners, marker_similarity, maxIdx] = select(matches);
    end


function [marker_corners, marker_similarity, maxIdx] = select(matches)

    marker_similarity = -Inf;
    marker_corners = [];
    maxIdx = 0;
    for i = 1:length(matches),
        if ~isempty(matches{i}) && matches{i}.similarity > marker_similarity,
            marker_similarity = matches{i}.similarity;
            marker_corners = matches{i}.pattern_corners;
            maxIdx = i;
        end
    end
