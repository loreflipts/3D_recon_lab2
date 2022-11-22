%% resampling and calculate epipolar line and matches in view 2

view1 = points_AP_rect_all';
view2 = Interpole_Discretise(points_LAT_rect_all', 2 * length(points_LAT_rect_all));

for k = 1:7
    name = "branch"+k+"";

    % add ones
    view1 = add_ones_2D(projection_AP_rect.(name));
    view2 = add_ones_2D(projection_LAT_rect.(name));

    view1_rect.(name) = projection_AP_rect.(name);

    % resample
    %view1_resampled = Interpole_Discretise(view1', 2 * length(view1))';
    view2_resampled = Interpole_Discretise(view2', 2 * length(view2))';

    % epipolar lines
    %epi_view1_rect.(name) = F_method3_rect * view2;
    epi_view2_rect.(name) = F_method3_rect * view1;

    % matching points for each epipolar line
    [matches_in_view2_using_view1_rect.(name), distances_view2_rect] = matching(epi_view2_rect.(name), view2_resampled);
    %[matches_in_view1_using_view2_rect.(name), distances_view1_rect] = matching(epi_view1_rect.(name), view1_resampled);
end

%% plotting of epipolar lines
x_range = [-10000:10000];
name_branch = "branch1";
plot_epi_view = epi_view2_rect;
plot_branch = projection_LAT_rect.(name_branch);

figure
for k = 1:5:length(plot_epi_view.(name_branch))
    y_range = (-plot_epi_view.(name_branch)(1,k) / plot_epi_view.(name_branch)(2,k)) .* x_range ...
        - (plot_epi_view.(name_branch)(3,k) / plot_epi_view.(name_branch)(2,k));
    hold on
    plot(x_range,y_range, Color='r')
    hold on
    plot(plot_branch(1,:), plot_branch(2,:), LineWidth=2, Color='b')
    axis equal
end

%% functions

function [branch_ones] = add_ones_2D(branch)
% 2D coordinate in a column, 2 x n matrix as input

branch_ones = [branch; ones(1,size(branch,2))];
end


function [matches, distances] = matching(epi_line, points_2D)
    matches = zeros(length(epi_line),3);
    distances = zeros(length(epi_line),1);
    for i = 1:length(epi_line)
        D_temp = zeros(length(points_2D),1);
        for k = 1:length(points_2D)
            D_temp(k) = abs(epi_line(:,i)' * points_2D(:,k));

        end
        [dist, idx] = min(D_temp);
        distances(i) = dist;
        matches(i,:) = points_2D(:,idx);
    end
    matches = matches';
end
