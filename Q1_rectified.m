load("Coronary.mat")
%close all 

%% create projections

[source_AP_rect, ~] = BuildViewGeom(1050, 750, 0.1953, 30, 25, [1024 1024]);
[source_LAT_rect, ~] = BuildViewGeom(1050, 750, 0.1953, -30, -25, [1024 1024]);

%% rectification 

[T1, T2, Pn1, Pn2] = Rectify(source_AP_rect.P, source_LAT_rect.P);

source_AP_rect.P = Pn1;
source_LAT_rect.P = Pn2;

projection_AP_rect = project_2D(source_AP_rect,  Coronary, 1);
projection_LAT_rect = project_2D(source_LAT_rect, Coronary, 1);

%% F calculation
[points_AP_rect_all, points_LAT_rect_all] = all_points_2D(projection_AP_rect, projection_LAT_rect);
F_method3_rect = FMatNorm8(points_AP_rect_all, points_LAT_rect_all);
F_method3_rect = F_method3_rect / F_method3_rect(3,3);

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
