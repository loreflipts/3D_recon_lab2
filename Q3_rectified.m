%% compute disparity

for k = 1:7
    name = "branch"+k+"";

    points_view1 = projection_AP_rect.(name);
    points_view2 = projection_LAT_rect.(name);

    matches_view1 = view1_rect.(name)(1:2, :);
    matches_view2 = matches_in_view2_using_view1_rect.(name)(1:2, :);

    diff_points_view1 = points_view1 - matches_view2;
    diff_points_view2 = points_view2 - matches_view1;

    disparity_points_view1_rect.(name) = vecnorm(diff_points_view1, 2, 1); % pixels
    disparity_points_view2_rect.(name) = vecnorm(diff_points_view2, 2, 1); 


    % compute depth
    focal_lenght_view1 = source_AP.f; % in pixels
    Z_points_view1_rect.(name) = focal_lenght_view1 * B ./ disparity_points_view1_rect.(name);

end



%% plot disparity and depth map

figure
hold on
for k = 1:7
    name = "branch"+k+"";
    plot_branch = projection_LAT.(name);
    scatter(plot_branch(1,:), plot_branch(2,:), [], disparity_points_view1_rect.(name))
    colorbar
end
title('Rectified image disparity map')

figure
hold on
for k = 1:7
    name = "branch"+k+"";
    plot_branch = projection_LAT.(name);
    scatter(plot_branch(1,:), plot_branch(2,:), [], Z_points_view1_rect.(name))
    colorbar
end
title('Rectified image depth map')
