
%% Add noise to (SID, SOD, DP, alpha, beta) 
SID = 1050; SOD = 750; DP = 0.1953; alpha1 = 30; alpha2 = -30; beta1 = 25; beta2 = -25;

num_trials = 5;
num_levels = 50;
RMSE_3D_corr = [RMSE_3D(:,1) zeros(4,num_levels)];
RMSE_2D_corr_view1 = [RMSE_2D_view1(:,1) zeros(4,num_levels)];
RMSE_2D_corr_view2 = [RMSE_2D_view2(:,1) zeros(4,num_levels)];

for level = 1:num_levels
    level
    for trial = 1:num_trials
    trial_name = "trial"+trial+"";
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        SID1_noisy = SID*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        SOD1_noisy = SOD*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        DP1_noisy = DP*(1 + rnum*level/100);

        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        SID2_noisy = SID*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        SOD2_noisy = SOD*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        DP2_noisy = DP*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        alpha1_noisy = alpha1*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        alpha2_noisy = alpha2*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        beta1_noisy = beta1*(1 + rnum*level/100);
        
        r = rand;
        rnum = r + (sign(r-0.5)-1)/2;
        beta2_noisy = beta2*(1 + rnum*level/100);
        
        [source_view1_noisy, ~] = BuildViewGeom(SID1_noisy, SOD1_noisy, DP1_noisy, alpha1_noisy, beta1_noisy, [1024 1024]);
        [source_view2_noisy, ~] = BuildViewGeom(SID2_noisy, SOD2_noisy, DP2_noisy, alpha2_noisy, beta2_noisy, [1024 1024]);

        % define inputs for refineCamParam
        K_noisy_view1 = source_view1_noisy.K;
        R_noisy_view1 = source_view1_noisy.R;
        T_noisy_view1 = source_view1_noisy.T;

        K_noisy_view2 = source_view2_noisy.K;
        R_noisy_view2 = source_view2_noisy.R;
        T_noisy_view2 = source_view2_noisy.T;


        %% Reconstruction in 3D using noisy dicom parameters
        
        P1_noisy = source_view1_noisy.P;
        P1_noisy = P1_noisy/P1_noisy(3,4);
        M1_noisy = reshape(P1_noisy',[12,1]);
        M1_noisy = M1_noisy(1:11);
        
        P2_noisy = source_view2_noisy.P;
        P2_noisy = P2_noisy/P2_noisy(3,4);
        M2_noisy = reshape(P2_noisy',[12,1]);
        M2_noisy = M2_noisy(1:11);
        
        for k = 1:7
            name = "branch"+k+"";
            recon3D_newnoisy.(name) = reconstruct_in3D(projection_view1.(name)(1,:), projection_view1.(name)(2,:), ...
                matches_in_view2_using_view1.(name)(1,:), matches_in_view2_using_view1.(name)(2,:), M1_noisy, M2_noisy);
            recon3D_newnoisy_ones.(name) = [recon3D_newnoisy.(name); ones(1,size(recon3D_newnoisy.(name),2))];
        end
        
        %% Correction of DICOM parameters

        [projection_view1_all, projection_matched_view2_all] = all_points_2D(projection_view1, matches_in_view2_using_view1);
        recon3D_noisy_all = all_points_3D(recon3D_newnoisy);

        [K_corr_view1, R_corr_view1, T_corr_view1, rperr_view1] = ...
            RefineCamParam(projection_view1_all(1:2,:), recon3D_noisy_all, K_noisy_view1, R_noisy_view1, T_noisy_view1, 2^(-52), 100);
        
        [K_corr_view2, R_corr_view2, T_corr_view2, rperr_view2] = ...
            RefineCamParam(projection_matched_view2_all(1:2,:), recon3D_noisy_all, K_noisy_view2, R_noisy_view2, T_noisy_view2, 2^(-52), 100);

        [SID_corr_view1, SOD_corr_view1, alpha_corr_view1, beta_corr_view1] = ...
            GetDicomFromCalib(K_corr_view1, R_corr_view1, T_corr_view1, DP);

        [SID_corr_view2, SOD_corr_view2, alpha_corr_view2, beta_corr_view2] = ...
            GetDicomFromCalib(K_corr_view2, R_corr_view2, T_corr_view2, DP);

        %[source_view1_corr, ~] = BuildViewGeom(SID_corr_view1, SOD_corr_view1, DP, alpha_corr_view1, beta_corr_view1, [1024 1024]);
        %[source_view2_corr, ~] = BuildViewGeom(SID_corr_view2, SOD_corr_view2, DP, alpha_corr_view2, beta_corr_view2, [1024 1024]);

        %% Reconstruction in 3D using corrected dicom parameters
        
        %P1_corr = source_view1_corr.P;
        P1_corr = K_corr_view1 * [R_corr_view1 T_corr_view1];
        P1_corr = P1_corr/P1_corr(3,4);
        M1_corr = reshape(P1_corr',[12,1]);
        M1_corr = M1_corr(1:11);
        
        %P2_corr = source_view2_corr.P;
        P2_corr = K_corr_view2 * [R_corr_view2 T_corr_view2];
        P2_corr = P2_corr/P2_corr(3,4);
        M2_corr = reshape(P2_corr',[12,1]);
        M2_corr = M2_corr(1:11);
        
        for k = 1:7
            name = "branch"+k+"";
            recon3D_corr.(name) = reconstruct_in3D(projection_view1.(name)(1,:), projection_view1.(name)(2,:), ...
                matches_in_view2_using_view1.(name)(1,:), matches_in_view2_using_view1.(name)(2,:), M1_corr, M2_corr);
            recon3D_corr_ones.(name) = [recon3D_corr.(name); ones(1,size(recon3D_corr.(name),2))];
        end

        projection_corr_view1 = newbackproject_2D(P1_corr, recon3D_corr_ones, 0);
        projection_corr_view2 = newbackproject_2D(P2_corr, recon3D_corr_ones, 0);
    
        [corr_view1_all, corr_view2_all] = all_points_2D(projection_corr_view1, projection_corr_view2);
        RMSE_2D_corr_view1(:,level+1) = RMSE_2D_corr_view1(:,level+1) + calc_RMSE_total(points_AP_all, corr_view1_all);
        RMSE_2D_corr_view2(:,level+1) = RMSE_2D_corr_view2(:,level+1) + calc_RMSE_total(points_LAT_all, corr_view2_all);
    
        recon3D_corr_all = all_points_3D(recon3D_corr);
        RMSE_3D_corr(:,level+1) = RMSE_3D_corr(:,level+1) + calc_RMSE_total(branches3D_all, recon3D_corr_all);
    end
    RMSE_2D_corr_view1(:,level+1) = RMSE_2D_corr_view1(:,level+1) / num_trials;
    RMSE_2D_corr_view2(:,level+1) = RMSE_2D_corr_view2(:,level+1) / num_trials;
    RMSE_3D_corr(:,level+1) = RMSE_3D_corr(:,level+1) / num_trials;
end

figure
plot((0:num_levels), RMSE_3D_corr(1,:), Color='r')
hold on
plot((0:num_levels), RMSE_3D_corr(2,:), Color='b')
hold on
plot((0:num_levels), RMSE_3D_corr(3,:), Color='g')
hold on
plot((0:num_levels), RMSE_3D_corr(4,:), Color=[0.9294    0.6941    0.1255])
legend('Total','X', 'Y', 'Z')
xlabel('Percentage of added noise to DICOM parameters [%]')
ylabel('RMSE')
title('Trial-averaged refined 3D RMS errors with respect to reference 3D model')

figure
plot((0:num_levels), RMSE_2D_corr_view1(1,:), Color='r')
hold on
plot((0:num_levels), RMSE_2D_corr_view1(2,:), Color='b')
hold on
plot((0:num_levels), RMSE_2D_corr_view1(3,:), Color='g')
legend('Total','U', 'V')
xlabel('Percentage of added noise to DICOM parameters [%]')
ylabel('RMSE')
title('Trial-averaged refined 2D RMS errors with respect to initial image of view 1')

figure
plot((0:num_levels), RMSE_2D_corr_view2(1,:), Color='r')
hold on
plot((0:num_levels), RMSE_2D_corr_view2(2,:), Color='b')
hold on
plot((0:num_levels), RMSE_2D_corr_view2(3,:), Color='g')
legend('Total','U', 'V')
xlabel('Percentage of added noise to DICOM parameters [%]')
ylabel('RMSE')
title('Trial-averaged refined 2D RMS errors with respect to initial image of view 2')

figure
plot((0:num_levels), RMSE_3D_corr(1,:), Color='r')
hold on
plot((0:num_levels), RMSE_2D_corr_view1(1,:), Color='b')
hold on
plot((0:num_levels), RMSE_2D_corr_view2(1,:), Color='g')
legend('Total 3D','Total 2D view 1', 'Total 2D view 2')
xlabel('Percentage of added noise to DICOM parameters [%]')
ylabel('RMSE')
title('Trial-averaged Refined Total RMS error')


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



function [error] = calc_RMSE_total(ideal, simple)
sum = 0;
sum_dim = zeros(size(ideal,1),1);

for k = 1:length(ideal)
    sum = sum + (norm(ideal(:,k)-simple(:,k)))^2;
    for dim = 1:size(ideal,1)
        sum_dim(dim) = sum_dim(dim) + (ideal(dim,k)-simple(dim,k))^2;

    end
end
error_3D = sqrt(sum/length(ideal));
error_dim = sqrt(sum_dim/length(ideal));
error = [error_3D; error_dim];
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

function [points] = all_points_3D(scan3D)
points = [];
    for k = 1:7
        name = "branch"+k+"";
    
        branch = scan3D.(name);

        points = [points [branch]];

    end
end