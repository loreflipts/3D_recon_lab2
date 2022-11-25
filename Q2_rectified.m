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
    view1_resampled = Interpole_Discretise(view1', 2 * length(view1))'; % ---
    view2_resampled = Interpole_Discretise(view2', 2 * length(view2))';

    % epipolar lines
    epi_view1_rect.(name) = F_method3_rect * view2; % ---
    epi_view2_rect.(name) = F_method3_rect * view1;

    % matching points for each epipolar line
    [matches_in_view2_using_view1_rect.(name), distances_view2_rect] = matching(epi_view2_rect.(name), view2_resampled);
    [matches_in_view1_using_view2_rect.(name), distances_view1_rect] = matching(epi_view1_rect.(name), view1_resampled); % ---
end

%% plotting of epipolar lines for view 1
x_range = -10000:10000;
name_branch = "branch1";
plot_epi_view = epi_view1_rect;
plot_branch = projection_AP_rect.(name_branch);

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
title('a) Epipolar lines for branch 1 in view 1')
% xlim([-6000 1000])
ylim([-3000 2000])

%% plotting of epipolar lines for view 2
x_range = -10000:10000;
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
title('Epipolar lines for branch 1 in view 2')
xlim([-6000 1000])
ylim([-3000 2000])

%% plotting projection with selected epipolar lines for view 2
% select points and store them as [X; Y; number in branch; branchnumber]
selected_points = zeros(4,1);
selected_points(:,1) = [projection_LAT_rect.branch1(:,100); 100; 1];
selected_points(:,2) = [projection_LAT_rect.branch5(:,50); 50; 5];

projection_LAT_rect = project_2D(source_LAT_rect, Coronary, 1);
hold on

for k = 1:size(selected_points,2)
    % print epipolar lane with length +- 200 of points x-position
    x_range = (selected_points(1,k)-200) : (selected_points(1,k)+200);
    % take corresponding epipolar line for selected point
    name_branch = "branch"+selected_points(4,k)+"";
    point_number = selected_points(3,k);
    y_range = (-plot_epi_view.(name_branch)(1,point_number) / plot_epi_view.(name_branch)(2,point_number)) .* x_range ...
        - (plot_epi_view.(name_branch)(3,point_number) / plot_epi_view.(name_branch)(2,point_number));
    % plot epipolar line
    plot(x_range,y_range,'Color',[0 (1/size(selected_points,2))*k 1])
    % plot corresponding point
    plot(selected_points(1,k), selected_points(2,k),'+','Color',[0 (1/size(selected_points,2))*k 1],MarkerSize=8)
    % plot matching point from other view
    plot(matches_in_view2_using_view1_rect.(name_branch)(1,point_number),...
        matches_in_view2_using_view1_rect.(name_branch)(2,point_number),'r+')
end
xlabel('u')
ylabel('v')
title('Projection of view 2 with selected epipolar lines')
ax = gca;
% zoomed section
axes('position',[.62 .175 .25 .3])
box on 
hold on
set(gca,'xtick',[],'ytick',[]);
x_range = (selected_points(1,k)-200) : (selected_points(1,k)+200);
y_range = (-plot_epi_view.branch1(1,100) / plot_epi_view.branch1(2,100)) .* x_range ...
        - (plot_epi_view.branch1(3,100) / plot_epi_view.branch1(2,100));
plot(x_range,y_range,'Color',[0 (1/size(selected_points,2))*1 1])
plot(projection_LAT_rect.branch1(1,:),projection_LAT_rect.branch1(2,:),'Color','#37a820')
plot(selected_points(1,1), selected_points(2,1),'+','Color',[0 (1/size(selected_points,2))*1 1],MarkerSize=8)
plot(matches_in_view2_using_view1_rect.branch1(1,100),matches_in_view2_using_view1_rect.branch1(2,100),'r+')
xlim([-3807.14, -3798.11])
ylim([190.985, 206.188])
yline(196.271,'b')

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
