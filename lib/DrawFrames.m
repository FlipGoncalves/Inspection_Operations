function H = DrawFrames(AA,P,F)
    % Draw an object defined by points P and faces F on each link
    % AA - hipermatrix of link transformations
    % P - matrix of object points
    % F - matrix of object faces
    % H - matrix of graphic handles
    
    H = cell(1, size(AA,3));
    
    % object at origin
    patch('Vertices', P(1:3,:)', 'Faces', F, 'Facecolor', 'w');
    
    Ak = eye(4);
    for i = 1:size(AA,3)
        Ak = Ak*AA(:,:,i);
        obj = Ak*P;
        H{i} = patch('Vertices', obj(1:3,:)', 'Faces', F, 'Facecolor', 'g');
    end
end

