function AAA = CalculateRobotMotion(MDH)
    N = size(MDH, 1);  % number of links
    M = size(MDH, 3);  % number of positions
    
    AAA = zeros(4,4,N,M);
    
    for i = 1:M
        AAA(:,:,:,i) = Tlinks(MDH(:,:,i));
    end

end

