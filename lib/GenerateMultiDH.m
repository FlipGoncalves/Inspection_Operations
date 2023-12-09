function MDH = GenerateMultiDH(DH,MQ,t)
    % Returns resolved DH matrices for several robot links' positions
    % DH - Base Denavit-Hartenberg matrix for the robot's zero hardware position
    % MQ - LinspaceVect(LinkVectorsI, LinkVectorsF, N)
    % t - binary vector of size equal to number of rows in DH. 0 means the
    % corresponding line in DH is a rotational link, 1 means it's prismatic.
    % Defaults to 0s if not provided
    % MDH - Hipermatrix of DH matrices defined for each MQ column
    
    if nargin < 3
        t = zeros(1, size(DH, 1));
    end
    
    nLinks = size(DH, 1);
    
    MDH = zeros(nLinks, 4, size(MQ, 2));

    for i = 1:size(MQ,2)
        MDH(:,:,i) = DH;
    
        for k = 1:nLinks
            if t(k)
                MDH(k,3,i) = MDH(k,3,i) + MQ(k,i);
            else
                MDH(k,1,i) = MDH(k,1,i) + MQ(k,i);
            end
        end
    end
end

