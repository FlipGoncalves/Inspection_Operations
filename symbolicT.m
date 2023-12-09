syms t1 t2 t4 t6 d DMIN t8 t9 LA LB LC LD LE LG LH LI

DH = [
    t1  0  LA  pi/2
    t2  LB  0  0
    -t2+pi/2  LC  0  0
    -t4+pi/2  LD  0  0
    t4-pi/2  LE  0  0
    t6-pi/2  0  0  0
    0  DMIN   0  0
    t8  LG 0  0
    t9  LH  0  0
];

T = eye(4);

for i = 1:size(DH, 1)
    T = T * (rotz(DH(i, 1))*trans(DH(i, 2),0,DH(i, 3))*rotx(DH(i, 4)));
end

% DH = [
%     0  0  LA  pi/2
%     0  LB  0  0
% 
%     pi/2  LC  0  0
%     pi/2  LD  0  0
% 
%     -pi/2  LE  0  0
%     0  0  0  pi/2
%     0  0   DMIN  -pi/2
%     -pi/2  LG 0  -pi/2
%     0  LH  0  0
% ];

%% T1
clear
close

syms t1 t2 LA LB

DH = [
    t1  0  LA  pi/2
    t2  LB  0  0
];

T1 = eye(4);
for i = 1:size(DH, 1)
    T1 = T1 * (rotz(DH(i, 1))*trans(DH(i, 2),0,DH(i, 3))*rotx(DH(i, 4)));
end

T1

%% T2
clear
close

syms t2 t4 LC LD

DH = [
    pi/2-t2  LC  0  0
    pi/2+t4  LD  0  0
];

T2 = eye(4);
for i = 1:size(DH, 1)
    T2 = T2 * (rotz(DH(i, 1))*trans(DH(i, 2),0,DH(i, 3))*rotx(DH(i, 4)));
end

T2

%% T3
clear
close

syms t4 t6 t9 LE LF LG LH

DH = [
    -t4-pi/2  LE  0  0         % -4
    t6  0  0  pi/2             % 6
    0  0   LF  -pi/2           % LF
    -t6-pi/2  LG 0  -pi/2      % -6
    t9  LH  0  0               % 9
];

T3 = eye(4);
for i = 1:size(DH, 1)
    T3 = T3 * (rotz(DH(i, 1))*trans(DH(i, 2),0,DH(i, 3))*rotx(DH(i, 4)));
end

T3

%% Tgeral
clear
close

syms t1 t2 LA LF

DH = [
    t1  0  LA  pi/2
    t2+pi/2  0  0  pi/2
    0  0  LF 0
];

T4 = eye(4);
for i = 1:size(DH, 1)
    T4 = T4 * (rotz(DH(i, 1))*trans(DH(i, 2),0,DH(i, 3))*rotx(DH(i, 4)));
end

T4

%% T5
clear
close

syms a b LA LB LC

DH = [
    a  0  LA  pi/2
    b  LB  0  0
    pi/2-b LC  0  0
];

T5 = eye(4);
for i = 1:size(DH, 1)
    T5 = T5 * (rotz(DH(i, 1))*trans(DH(i, 2),0,DH(i, 3))*rotx(DH(i, 4)));
end

T5

%% T6
clear
close

syms a b LA LB LC LD LE

DH = [
    b  0  LA  pi/2
    0  LB  0  0

    pi/2  LC  0  0
    pi/2+a  LD  0  0

    -pi/2-a  LE  0  0
];

T6 = eye(4);
for i = 1:size(DH, 1)
    T6 = T6 * (rotz(DH(i, 1))*trans(DH(i, 2),0,DH(i, 3))*rotx(DH(i, 4)));
end

T6




