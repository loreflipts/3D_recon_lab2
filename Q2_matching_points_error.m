%% without rectification
error_distance = [];

for k = 1:7
    name = "branch"+k+"";
    
    for j = 1:length(projection_LAT.(name))
        error_distance = [error_distance; ...
            ((projection_LAT.(name)(1,j)-matches_in_view2_using_view1.(name)(1,j))^2 + ...
            (projection_LAT.(name)(2,j)-matches_in_view2_using_view1.(name)(2,j))^2)];
    end
        
end

RMSE_points = sum(error_distance/mean(error_distance.^2))

%% with rectification
error_distance = [];

for k = 1:7
    name = "branch"+k+"";
    
    for j = 1:length(projection_LAT.(name))
        error_distance = [error_distance; ...
            ((projection_LAT_rect.(name)(1,j)-matches_in_view2_using_view1_rect.(name)(1,j))^2 + ...
            (projection_LAT_rect.(name)(2,j)-matches_in_view2_using_view1_rect.(name)(2,j))^2)];
    end
        
end

RMSE_points = sum(error_distance/mean(error_distance.^2))