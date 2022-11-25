load("Coronary.mat")
close all 

%% create projections

[source_AP, ~] = BuildViewGeom(1050, 750, 0.1953, 30, 25, [1024 1024]);
projection_AP = project_2D(source_AP,  Coronary, 1);

[source_LAT, ~] = BuildViewGeom(1050, 750, 0.1953, -30, -25, [1024 1024]);
projection_LAT = project_2D(source_LAT, Coronary, 1);


%% method 1

K_AP = source_AP.K;
T_AP = source_AP.T;
R_AP = source_AP.R;

K_LAT = source_LAT.K;
T_LAT = source_LAT.T;
R_LAT = source_LAT.R;

A_AP = [R_AP T_AP; 0 0 0 1];
A_LAT = [R_LAT T_LAT; 0 0 0 1];

As = A_LAT * inv(A_AP);

Ts = As(:,end);
B = norm(Ts,2); % distance between camera's, needed for Q3
Rs = As(1:3,1:3);

Ts_x = [0 -Ts(3) Ts(2); Ts(3) 0 -Ts(1); -Ts(2) Ts(1) 0];

E = Ts_x * Rs;

F_method1 = inv(K_LAT)' * E * inv(K_AP);
F_method1 = F_method1 / F_method1(3,3);


%% method 2

[points_AP_bif, points_LAT_bif] = choose_points_2D(projection_AP, projection_LAT);
F_method2 = FMatNorm8(points_AP_bif, points_LAT_bif);
F_method2 = F_method2 / F_method2(3,3);


%% method 3

[points_AP_all, points_LAT_all] = all_points_2D(projection_AP, projection_LAT);
F_method3 = FMatNorm8(points_AP_all, points_LAT_all);
F_method3 = F_method3 / F_method3(3,3);


%% Save 3D in branches

scan_3D = Coronary.skeleton;

for k = 1:7
    name = "branch"+k+"";
    branch3D = scan_3D(k);
    branch3D = branch3D{1,1}';
    %branch3D = [branch3D; ones(1,size(branch3D,2))];
    branches3D.(name) = branch3D;
end

%% functions
function [points_bifurcation1, points_bifurcation2] = choose_points_2D(image1, image2)

points_bifurcation1 = zeros(3,14);
points_bifurcation2 = zeros(3,14);

for k = 1:7
    name = "branch"+k+"";
    
    branch1 = image1.(name);
    branch2 = image2.(name);

    points_bifurcation1(:,2*k-1) = [branch1(:, 1); 1];
    points_bifurcation1(:,2*k) = [branch1(:, end); 1];
    points_bifurcation2(:,2*k-1) = [branch2(:, 1); 1];
    points_bifurcation2(:,2*k) = [branch2(:, end); 1];

end
points_bifurcation1 = unique(points_bifurcation1','rows','stable')';
points_bifurcation2 = unique(points_bifurcation2','rows','stable')';

% pick random points instead of bifurcations

% points1 = zeros(3,8);
% points2 = zeros(3,8);

% for k = 1:8
%     
%     branch_index = randi(7);
%     name = "branch"+branch_index+"";
%     
%     branch1 = image1.(name);
%     branch2 = image2.(name);
% 
%     point_index = randi(size(branch1,1));
%     points1(:,k) = [branch1(:, point_index)' 1];
%     points2(:,k) = [branch2(:, point_index)' 1];
% 
% end
end


function [points1, points2] = all_points_2D(image1, image2)
points1 = [];
points2 = [];
    for k = 1:7
        name = "branch"+k+"";
    
        branch1 = image1.(name);
        branch2 = image2.(name);

        points1 = [points1 [branch1; ones(1,size(branch1,2))]];
        points2 = [points2 [branch2; ones(1,size(branch2,2))]];

    end
end
