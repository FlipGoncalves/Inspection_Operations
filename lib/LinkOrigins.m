function Org = LinkOrigins(AA)
    % Returns matrix of coordinate systems' origins given the transformations
    % hipermatrix of links
    % AA - Transformations hipermatrix
    % Org - Matrix of origins, with each origin given as a column, starting
    % with [0 0 0]'
    
    Org = zeros(3, size(AA,3)+1);
    
    Ak = eye(4);
    
    for i = 1:size(AA,3)
        Ak = Ak*AA(:,:,i);
        Org(:,i+1) = Ak(1:3,4);
    end
end
