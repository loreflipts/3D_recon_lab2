%% Save 3D in branches

scan_3D = Coronary.skeleton;

for k = 1:7
    name = "branch"+k+"";
    branch3D = scan_3D(k);
    branch3D = branch3D{1,1}';
    %branch3D = [branch3D; ones(1,size(branch3D,2))];
    branches3D.(name) = branch3D;
end