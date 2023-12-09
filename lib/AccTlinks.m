function AAc = AccTlinks(DH)
    % DH - robot's Denavit-Hartenberg matrix
    % AAc - accumulated geometric transformations

    AA = Tlinks(DH);

    AAc = zeros(4, 4, size(DH, 1));
    AAc(:,:,1) = AA(:,:,1);
    for i=2:size(DH, 1)
        AAc(:,:,i) = AAc(:,:,i-1)*AA(:,:,i);
    end
end

