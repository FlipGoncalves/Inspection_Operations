function Q = invkinPROJETO(x,y,z,LA,LB,LC,LD,LE,DMIN,DMAX,LG,LH,A)

    % Verificar se o ponto pertence ao Espaço de Trabalho
    xymax = LB+LD+DMAX+LG+LH;
    zmax = LA+LB+LC+LD+LE+DMAX;
    inside_max = (x^2 / xymax^2) + (y^2 / xymax^2) + (z^2 / zmax^2) <= 1;

    if inside_max == 0
        Q = Nan;
        return
    end
    
    % theta 1 -> RR a 3D
    theta1 = atan2(y, x);

    % V2 = [x y] - A;
    % V2 = V2/norm(V2);
    % PLH = [x y] + V2*LH;
    % x = PLH(1);
    % y = PLH(2);

    % theta 9
    % P2 = A - [0 0];
    % P1 = [x y] - [0 0];
    % theta9 = -acos(dot(P1,P2) / (norm(P1) * norm(P2)));
    theta9 = 0;

    % Verificar se é preciso usar as Juntas 2/3 e 4/5
    V = [x, y, z] - [(LB-LD)*cos(theta1) (LB-LD)*sin(theta1) LA+LC+LE];

    % theta 2
    theta2 = 0;

    % theta 4
    theta4 = 0;

    if norm(V)-LG-LH > DMAX
        theta2 = -pi/7.5;
        theta4 = -pi/2;
    end

    % theta 6
    z6 = z-LA-LC-LE-(LB*sin(theta2)-LD*sin(theta4));
    x6 = (x-LG-LH-(LB*cos(theta2)-LD*cos(theta4)))*cos(theta1);
    y6 = (y-LG-LH-(LB*sin(theta2)-LD*sin(theta4)))*sin(theta1);
    theta6 = atan2(z6, x6);

    % LF
    LF = sqrt(z6.^2 + x6.^2 + y6.^2);
    d = LF - DMIN;
    
    % Return
    Q = [
        theta1
        theta2
        -theta2
        theta4
        -theta4
        theta6
        d
        -theta6
        theta9
        ];

end

