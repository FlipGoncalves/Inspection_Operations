function [H, h, P, AAA] = InitRobot(QQ, NN, DH, jTypes, sScale)
    % QQ - matriz de colunas das posições de junta (pelo menos duas)
    % NN - número de pontos de cada segmento do movimento
    % DH - Matriz-base dos parâmetros cinemáticos
    % jTypes - vetor com o tipo de juntas (0=rot, 1=prism; opcional)
    % sScale - Fator de escala dos seixos3() (opcional: 1 por defeito)
    
    if nargin < 4
        jTypes = zeros(1, size(DH, 1));
    end
    if nargin < 5
        sScale = 1;
    end
    
    [P, F] = seixos3(sScale);
    
    % first AAA
    MQ = LinspaceVect(QQ(:,1), QQ(:,2), NN);
    MDH = GenerateMultiDH(DH,MQ,jTypes);
    AAA = CalculateRobotMotion(MDH);
    
    Org = LinkOrigins(AAA(:,:,:,1));
    h = DrawLinks(Org);
    H = DrawFrames(AAA(:,:,:,1), P, F);
    
    for i = 2:(size(QQ, 2)-1)
        MQ = LinspaceVect(QQ(:,i), QQ(:,i+1), NN);
        MDH = GenerateMultiDH(DH,MQ,jTypes);
        AAA = cat(4, AAA, CalculateRobotMotion(MDH));
    end
end

