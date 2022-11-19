%% Ground truth calculation

[source_LAT, ~] = BuildViewGeom(1050, 750, 0.1953, -30, -25, [1024 1024]);

K_LAT = source_LAT.K;
T_LAT = source_LAT.T;
R_LAT = source_LAT.R;

A_LAT = [R_LAT T_LAT; 0 0 0 1];

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
    ground_truth.(name) = inv(A_LAT) * branch_scan;
    plot_branch = projection_LAT.(name);
    scatter(plot_branch(1,:), plot_branch(2,:), [], ground_truth.(name)(3,:))
    colorbar
end
title('Ground truth depth map')
