%% Compare depth maps with and without rectification to groud truth
figure
tiledlayout(1,2)

% error without rectification
h(1) = nexttile;
hold on
for k = 1:7
    name = "branch"+k+"";
    plot_branch = projection_LAT.(name);
    % calculate error and normalize
    depth_error_view1.(name) = abs(Z_points_view1.(name) - abs(ground_truth.(name)(3,:)));
%     depth_error_view1.(name) = depth_error_view1.(name)/norm(depth_error_view1.(name));
    scatter(plot_branch(1,:), plot_branch(2,:), [], depth_error_view1.(name))
    colormap('jet')
    colorbar
end
title('a) Depth error without rectification')

% error with rectification
h(2) = nexttile;
hold on
for k = 1:7
    name = "branch"+k+"";
    plot_branch = projection_LAT.(name);
    % calculate error and normalize
    depth_error_view1_rect.(name) = abs(Z_points_view1_rect.(name) - abs(ground_truth.(name)(3,:)));
%     depth_error_view1_rect.(name) = depth_error_view1_rect.(name)/norm(depth_error_view1_rect.(name));
    scatter(plot_branch(1,:), plot_branch(2,:), [], depth_error_view1_rect.(name))
    colormap('jet')
    colorbar
end
title('b) Depth error with rectification')

% define colorbar
% set(h, 'Colormap', jet, 'CLim', [0 1200000])
% cbh = colorbar(h(end)); 
% cbh.Layout.Tile = 'east'; 