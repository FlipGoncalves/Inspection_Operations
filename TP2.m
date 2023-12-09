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
      0     treetop(1)-D
      0     0                             
      0     treetop(3)             
];

% Cinematica Inversa - 1º Segmento
Q = zeros([9, size(Points, 2)]);
for i = 2:size(Points, 2)
    Qi = invkinPROJETO(Points(1,i), Points(2,i), Points(3,i), LA, LB, LC, LD, LE, LF, DMAX, LG, LH, Arvore);
    Q(:,i) = Qi(:,1);  % cotovelo em baixo
end

% normalizar com juntas virtuais
QQ = [Q; zeros(2, size(Points, 2))];
QQ = LinspaceVect(QQ(:,1), QQ(:,2), N);

% Cinemática Diferencial - Segmentos Lineares e Descida em Zig Zag
Points_Linear = [
    %     A                B                   C     
      treetop(1)-D      treetop(1)-TB-D   treetop(1)-D
      0                 0                 0     
      treetop(3)        PH                treetop(3)
];

% Descida em Zig-Zag
step = TH/5;
theta1 = linspace(pi/7+pi, pi/3.5+pi, 5)';
theta2 = linspace(-pi/7+pi, -pi/3.5+pi, 5)';
theta = [theta1 theta2];

Points_ZZ = zeros(3, 10);

for i=1:5
    if mod(i, 2) == 0
        theta(i, :) = flip(theta(i, :));
    end
    x = (i*TB/5+D)*cos(theta(i, :)) + tree_distance;
    y = (i*TB/5+D)*sin(theta(i, :));
    z = TH+PH-i*step;
    pts = [x(1) x(2)
           y(1) y(2)
           z    z   ];
    Points_ZZ(:,i*2-1:i*2) = pts;
end

% Todos os Pontos
dq = zeros(11, N*((size(Points_Linear, 2) + size(Points_ZZ, 2))-1));
init = [invkinPROJETO(Points_Linear(1,1), Points_Linear(2,1), Points_Linear(3,1), LA, LB, LC, LD, LE, LF, DMAX, LG, LH, Arvore); 0; 0];

% Movimento Linear Inicial
ind = 1;
for pt=2:size(Points_Linear, 2)
    dr = [(Points_Linear(:,pt) - Points_Linear(:,pt-1))/N; 0; 0; 0];
    
    for i=1:N
        JJ = jacobianGeom(DH, init, jtypes);
        q = pinv(JJ) * dr;
        init = init + q;
        dq(:,ind) = init;
        ind = ind + 1;
    end
end

Points_ZZ = [Points_Linear(:,end) Points_ZZ];

% Movimento Zig Zag
for pt=2:size(Points_ZZ, 2)
    if mod(pt-1, 2) == 0
        C = [TD+TB 0];
        index = round(pt/2)-1;
        stp = theta(index, 2)/N;
        if theta(index, 2) < theta(index, 1)
            stp = -stp;
        end
        teta = linspace(theta(index, 1), theta(index, 2), N);
        x = C(1) + (index*TB/5+D)*cos(teta);
        y = C(2) + (index*TB/5+D)*sin(teta);
        z = linspace(TH+PH-index*step, TH+PH-index*step, N);
        dr = diff([[x x(end)]; [y y(end)]; [z z(end)]; zeros(1,N+1); zeros(1,N+1); zeros(1,N+1)]');

        for i=1:N
            JJ = jacobianGeom(DH, init, jtypes);
            q = pinv(JJ) * dr(i,:)';
            init = init + q;
            dq(:,ind) = init;
            ind = ind + 1;
        end
    else
        dr = [(Points_ZZ(:,pt) - Points_ZZ(:,pt-1))/N; 0; 0; 0];

        for i=1:N
            JJ = jacobianGeom(DH, init, jtypes);
            q = pinv(JJ) * dr;
            init = init + q;
            dq(:,ind) = init;
            ind = ind + 1;
        end
    end
end

% Cinematica Inversa - Último Segmento
AA = Tlinks(DH);
T = eye(4);
for A=1:size(DH, 1)
    T = T * AA(:,:,A);
end

Last_Point = invkinPROJETO(T(1,4), T(2,4), T(3,4), LA, LB, LC, LD, LE, LF, DMAX, LG, LH, Arvore);

% normalizar com juntas virtuais
Last_Point = [Last_Point; 0; 0];
Last_Point = LinspaceVect(dq(:,end), Last_Point, N);

QQ = [QQ dq Last_Point];

%%% Gráficos
% Representar LF num gráfico
LF_graph = QQ(7,:) + LF;
subplot(1,2,1);
hold on;
grid on;
plot(LF_graph);
plot(linspace(LF, LF, size(LF_graph, 2)));
plot(linspace(DMAX, DMAX, size(LF_graph, 2)));
legend(["LF", "LF(min)", "LF(max)"]);

% Robo e Arvore de Natal
subplot(1,2,2);
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

% Desenhar o Path
M = size(AAA, 4);
for m = 1:M
    Org = LinkOrigins(AAA(:,:,:,m));

    Ak = eye(4);
    for n = 1:size(AAA, 3)
        Ak = Ak*AAA(:,:,n,m);
        Pn = Ak*P;
    end

    plot3(Ak(1,4), Ak(2,4), Ak(3,4), '.r');
end

% Animar o Robo
pause(1);
AnimateRobot(H, AAA, P, h, 0.01);
