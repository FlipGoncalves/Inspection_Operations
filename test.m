clear 
close

%%% Variaveis Defeito
LA = 940;
LB = 1850;
LC = 400;
LD = 1600;
LE = 300;
LF = 1800;
LG = 970;
LH = 1150;
LI = 1000;
TB = 3400;
TH = 4000;
PH = 1000;
PD = 400;
D = 500;

%%% Ficheiro de Variáveis
try
    tfid = fopen("tp2.txt");
    tdata = textscan(tfid,"%s=%s");
    fclose(tfid);

    if( numel(tdata{1}) ~= numel(tdata{2}))
        disp("Error reading file. Missing = !")
        clear tdata tfid
    else
        ndata={ tdata{1} repmat("=", size(tdata{1})) tdata{2}};
        sdata=strcat(ndata{1},ndata{2},ndata{3});
        for i=1:numel(sdata)
            try
                eval(sdata{i});
            catch
                sprintf("Bad format in line %d of data file!",i)
            end
        end
        clear i tfid ndata tdata sdata
    end

catch
    disp("Cannot open file.");
end

% Variaveis com valores mais próximos de 0
LA = LA/1000;
LB = LB/1000;
LC = LC/1000;
LD = LD/1000;
LE = LE/1000;
LF = LF/1000;
LG = LG/1000;
LH = LH/1000;
LI = LI/1000;
TB = TB/1000;
TH = TH/1000;
PH = PH/1000;
PD = PD/1000;
D = D/1000;

% Tamanho Maximo da Junta Prismatica
DMAX = LF*1.7;

%%% Árvore de Natal
TB = TB/2;
PD = PD/2;
TD = LF + (LB-LD) + LG + LH + 0.2;
tree_distance = TD + TB;

% Cone de Folhas
[xTree, yTree, zTree] = cylinder([0, TB]);
zTree = zTree * -TH + TH + PH;
xTree = xTree + tree_distance;

% Cilindro do Tronco
[xLog, yLog, zLog] = cylinder([PD, PD]);
zLog = zLog * PH;
xLog = xLog + tree_distance;

%%% DH
DH = [
    0  0  LA  pi/2
    0  LB  0  0

    pi/2  LC  0  0
    pi/2  LD  0  0

    -pi/2  LE  0  0
    0  0  0  pi/2
    0  0   LF  -pi/2
    -pi/2  LG 0  -pi/2
    0  0   -LI  0            % virtual
    0  LH  0   0             % virtual
    0  0   LI  0             % virtual
];


%%% Movimento

N = 300;

% Juntas
jtypes = [0 0    0 0    0 0 1 0    0 0 0];

treetop = [TD+TB 0 TH+PH];
Arvore = treetop(1:2);

% Pontos
Points = [
    % I     A                B                     
      0     treetop(1)-D        treetop(1)
      0     0                   D
      0     treetop(3)          treetop(3)-1  
];

% Cinematica Inversa - 1º Segmento
Q = zeros([9, size(Points, 2)]);
for i = 2:size(Points, 2)
    Qi = invkinPROJETO(Points(1,i), Points(2,i), Points(3,i), LA, LB, LC, LD, LE, LF, DMAX, LG, LH, Arvore);
    Q(:,i) = Qi(:,1);  % cotovelo em baixo
end

% normalizar com juntas virtuais
QQ = [Q; zeros(2, size(Points, 2))];
QQ = LinspaceVect(QQ(:,1), QQ(:,3), N);

% Robo e Arvore de Natal
axis equal;
axis([-5 10 -5 5 0 6]);
view(120,30)
hold on;                    
grid on;                                
xlabel('X');
ylabel('Y');
zlabel('Z');

% Arvore de Natal
surf(xTree, yTree, zTree, 'FaceAlpha', 1, 'EdgeColor', 'black', 'FaceColor', 'green');
a = surf(xLog, yLog, zLog, 'FaceAlpha', 1, 'EdgeColor', 'black');
a.FaceColor = "#D95319";

% Iniciar o Robo
[H,h,P,AAA] = InitRobot(QQ, 1, DH, jtypes, 0.2);

% Animar o Robo
pause();
AnimateRobot(H, AAA, P, h, 0.01, 1);
