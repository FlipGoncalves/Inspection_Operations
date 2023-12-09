function JJ = jacobianGeom(DH, MQ, jTyp)
    % DH - Matriz base de DH do robˆo
    % MQ - Matriz de todos os vetores de juntas a percorrer
    % jTyp - Tipos de juntas (0-rotational, 1-linear)
    % JJ - Hipermatriz com todos os jacobianos por MQ (6xN)

    MDH = GenerateMultiDH(DH, MQ, jTyp);
    N = numel(jTyp);

    JV = zeros(3, N);
    JW = JV;
    JJ = zeros(6, N, size(MDH, 3));
    for k=1:size(MDH, 3)
        AAc = AccTlinks(MDH(:,:,k));
        Org = [[0;0;0] squeeze(AAc(1:3, 4, :))];
        Zis = [[0;0;1] squeeze(AAc(1:3, 3, :))];
        ON = Org(:, end);

        for i=1:N
            Oi1 = Org(:,i);
            zi1 = Zis(:,i);
            
            if jTyp(i) == 0
                Jvi = cross(zi1, ON - Oi1);
                Jwi = zi1;
            else
                Jvi = zi1;
                Jwi = 0;
            end
            JV(:,i) = Jvi;
            JW(:,i) = Jwi;
        end

        JJ(:,:,k) = [JV; JW];
    end
end

