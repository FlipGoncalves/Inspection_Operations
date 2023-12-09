function AA = Tlinks(DH)
    % Returns a tranformation for each link in DH, as a 3 dimention hipermatrix
    % DH - Denavit-Harternberg matrix, where each row is a link [teta, l, d, alpha]

    AA = zeros(4, 4, size(DH, 1));
    for i = 1:size(DH,1)
        AA(:,:,i) = Tlink(DH(i,1), DH(i,2), DH(i,3), DH(i,4));
    end
end
