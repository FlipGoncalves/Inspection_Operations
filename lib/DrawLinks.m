function h = DrawLinks(Org)
    % Draws graphic representing the robot's links given the coordinates of the
    % each link's origin
    % Org - matrix of link origins, with each column being a link
    
    h = plot3(Org(1,:), Org(2,:), Org(3,:), '-o', 'LineWidth',2);
end

