%% Ground truth calculation

[source_LAT, ~] = BuildViewGeom(1050, 750, 0.1953, -30, -25, [1024 1024]);
K_LAT = source_LAT.K;
T_LAT = source_LAT.T;
R_LAT = source_LAT.R;
A_LAT = [R_LAT T_LAT; 0 0 0 1];

% [source_AP, ~] = BuildViewGeom(1050, 750, 0.1953, 30, 25, [1024 1024]);
% K_AP = source_AP.K;
% T_AP = source_AP.T;
% R_AP = source_AP.R;
% A_AP = [R_AP T_AP; 0 0 0 1];


scan_3D = Coronary.skeleton;
branch1_scan = scan_3D(1);
branch1_scan = branch1_scan{1,1};

figure
hold on
for k = 1:7
    name = "branch"+k+"";
    branch_scan = scan_3D(k);
    branch_scan = branch_scan{1,1}';
    branch_scan = [branch_scan; ones(1,size(branch_scan,2))];
    ground_truth.(name) = branch_scan - [source_LAT.worldPos; 1]; %inv(A_AP) * branch_scan;          
    plot_branch = projection_LAT.(name);
    scatter(plot_branch(1,:), plot_branch(2,:), [], ground_truth.(name)(3,:))
    colormap('summer')
    colorbar
end
title('Ground truth depth map')

%% Calculate variance for view 1
temp = [];
for k = 1:7
    name = "branch"+k+"";
    temp = [temp, disparity_points_view1_rect.(name)];
end
disparity_variance_view1_rect = var(temp);
clear temp