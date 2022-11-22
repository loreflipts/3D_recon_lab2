function[projections] = newbackproject_2D(P, recon_3D, plotting)


if plotting == 1
    figure
end

for k = 1:7
    name = "branch"+k+"";
    %projection = K * RT * branch;
    projection = P * recon_3D.(name);
    
    for i = 1:size(projection,2)
        projection(:,i) = projection(:,i)./projection(3,i);
    end
    
    projections.(name) = [projection(1,:); projection(2,:)];

    if plotting == 1
        plot(projection(1,:),projection(2,:))
        axis equal
        hold on
    end
end

end 