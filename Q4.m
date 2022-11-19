%% Reconstruction
P1 = source_AP.P;
P1 = P1/P1(3,4);
M1 = reshape(P1',[12,1]);
M1 = M1(1:11);

P2 = source_LAT.P;
P2 = P2/P2(3,4);
M2 = reshape(P2',[12,1]);
M2 = M2(1:11);

figure
hold on
for k = 1:7
    name = "branch"+k+"";
    recon3D_original.(name) = reconstruct_in3D(projection_AP.(name)(1,:), projection_AP.(name)(2,:), ...
        projection_LAT.(name)(1,:), projection_LAT.(name)(2,:), M1, M2);

    recon3D_matchedpoints.(name) = reconstruct_in3D(matches_in_view1_using_view2.(name)(1,:), matches_in_view1_using_view2.(name)(2,:), ...
        matches_in_view2_using_view1.(name)(1,:), matches_in_view2_using_view1.(name)(2,:), M1, M2);

    scatter3(recon3D_original.(name)(1,:),recon3D_original.(name)(2,:), recon3D_original.(name)(3,:), MarkerEdgeColor='r',LineWidth=1,Marker='*')
    hold on
    scatter3(recon3D_matchedpoints.(name)(1,:),recon3D_matchedpoints.(name)(2,:), recon3D_matchedpoints.(name)(3,:),MarkerEdgeColor='g',LineWidth=1,Marker='*')
    hold on
    scatter3(branches3D.(name)(1,:),branches3D.(name)(2,:), branches3D.(name)(3,:),MarkerEdgeColor='b',LineWidth=1,Marker='.')
    view(3)

end



%% Backprojection in 2D

for k = 1:7
    name = "branch"+k+"";
    recon3D_matchedpoints_ones.(name) = [recon3D_matchedpoints.(name); ones(1,size(recon3D_matchedpoints.(name),2))];
end

projection_matched_view1 = backproject_2D(source_AP, recon3D_matchedpoints_ones, 0);
projection_matched_view2 = backproject_2D(source_LAT, recon3D_matchedpoints_ones, 0);

[matched_view1_all, matched_view2_all] = all_points_2D(projection_matched_view1, projection_matched_view2);
RMSE_2D_view1 = calc_RMSE(points_AP_all, matched_view1_all);
RMSE_2D_view2 = calc_RMSE(points_LAT_all, matched_view2_all);

%% Functions

function [p3D] = reconstruct_in3D(u1,v1,u2,v2,M1,M2)
% gives back reconstructed 3D points the 3D Points for ONE vertebra based on
% calibration matrices for both views (M_LAT & M_PA0) and the 
% 2D coordinates for the vertebra in view LAT and view PA0

N = size(u1,2);
p3D = zeros(3,N);

for k=1:N
    A = [M1(1) - M1(9) * u1(k)   M1(2) - M1(10) * u1(k)     M1(3) - M1(11) * u1(k);
        M1(5) - M1(9) * v1(k)    M1(6) - M1(10) * v1(k)     M1(7) - M1(11) * v1(k);
        M2(1) - M2(9) * u2(k)    M2(2) - M2(10) * u2(k)     M2(3) - M2(11) * u2(k);
        M2(5) - M2(9) * v2(k)    M2(6) - M2(10) * v2(k)     M2(7) - M2(11) * v2(k)];
    
    b = [-M1(4)+u1(k); -M1(8)+v1(k);
        -M2(4)+u2(k); -M2(8)+v2(k)];
    
    % reconstruct point k out of N in 3D for the vertebrae
    p3D(:,k) = pinv(A)*b;
end

end



function [error] = calc_RMSE(ideal, simple)
sum = 0;

for k = 1:length(ideal)
    sum = sum + (norm(ideal(:,k)-simple(:,k)))^2;
end
RMSE3 = sqrt(sum/length(ideal));

error = [RMSE3];

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