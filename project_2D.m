function[projections] = project_2D(source, scan_3D, plotting)

RT = [source.R,source.T];
K = source.K;

branches = scan_3D.skeleton;

if plotting == 1
    figure
end

for k = 1:size(branches,1)
    branch = branches(k);
    branch = branch{1,1}';
    branch = [branch; ones(1,size(branch,2))];
    
    name = "branch"+k+"";
    %projection = K * RT * branch;
    projection = source.P * branch;
    
    for i = 1:size(projection,2)
        projection(:,i) = projection(:,i)./projection(3,i);
    end
    
    projections.(name) = [projection(1,:); projection(2,:)];

    if plotting == 1
        plot(projection(1,:),projection(2,:),'Color','#37a820')
        axis equal
        hold on
    end
end

end 