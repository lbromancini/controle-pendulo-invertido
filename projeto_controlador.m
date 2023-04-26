clc, clear

%% Linearização

% variáveis auxiliares
syms x1, syms x2, syms x3, syms x4

% variáveis do sistema

% pênduloS
syms I, syms m, syms l, syms bp, syms g

% carrinho
syms kt, syms vi, syms ra, syms r, syms bc, syms M

I = (m*l^2)/3; % Momento de inércia

x1o = x3;
x2o = x4;
x4o = (((kt*vi)/(ra*r))-(bc*x4))/(M+m);
x3o = ((-m*l*x4o*cos(x1))-(bp*x3)+(m*l*g*sin(x1)))/(I+m*l^2);

% elementos da matriz jacobiana

a11 = diff(x1o,x1);
a12 = diff(x1o,x2);
a13 = diff(x1o,x3);
a14 = diff(x1o,x4);

a21 = diff(x2o,x1);
a22 = diff(x2o,x2);
a23 = diff(x2o,x3);
a24 = diff(x2o,x4);

a31 = diff(x3o,x1);
a32 = diff(x3o,x2);
a33 = diff(x3o,x3);
a34 = diff(x3o,x4);

a41 = diff(x4o,x1);
a42 = diff(x4o,x2);
a43 = diff(x4o,x3);
a44 = diff(x4o,x4);

b11 = diff(x1o,vi);
b21 = diff(x2o,vi);
b31 = diff(x3o,vi);
b41 = diff(x4o,vi);

% matriz A linear
A_nl = [a11 a12 a13 a14 ; a21 a22 a23 a24 ; 
     a31 a32 a33 a34 ; a41 a42 a43 a44];
 
% matriz B linear
B_nl = [b11 b21 b31 b41]';

%% Variáveis numéricas

% Transforma a função de matemática simbólica em uma função normal      
A_in = inline(A_nl,'x1','x2','x3','x4','I','m','l','bp','g','kt','vi','ra','r','bc','M');
B_in = inline(B_nl,'x1','x2','x3','x4','I','m','l','bp','g','kt','vi','ra','r','bc','M');

% variáveis do pêndulo

m = 0.2065; % massa do pêndulo;
l = 0.295; % distância do eixo do pêndulo ao centro de massa
bp = 0.0016733; % coeficiente de atrito do pêndulo
g = 9.807; % aceleração da gravidade

% variáveis do carro

ra = 3.5; % resistência da armadura do motor
r = 0.0065; % raio da polia
M = 0.5705; % massa do carro
kt = 0.03865; % constante de torque do motor
bc = 1.6408; % coeficiente de atrito do carrinho

% condições iniciais

x1 = 0;
x2 = 0;
x3 = 0;
x4 = 0;
vi = 0; % tensão no motor

% matrizes com os valores numéricos
     
A = A_in(x1,x2,x3,x4,I,m,l,bp,g,kt,vi,ra,r,bc,M)
B = B_in(x1,x2,x3,x4,I,m,l,bp,g,kt,vi,ra,r,bc,M)
C = [1 0 0 0 ; 0 1 0 0]

%% estabilidade
eig(A)

%% teste de controlabilidade
Co = [B A*B (A^2)*B (A^3)*B]; 
determinante_Co = det(Co)

%% teste de observabilidade
Ob = [C; C*A; C*(A^2); C*(A^3)];
n = rank(Ob)

%% projeto de controlador 

%expansão das matrizes
Ca = [0 1 0 0];
A_exp = [ A zeros(4,1); -Ca 0 ];
B_exp = [B; 0];

%controlador por LQR

%Ajustes
x1_max = pi/18; %10; %Variação da Posição Angular em radianos
x2_max = 0.1;  %Variação da posição linear em mm
x3_max = 1; %Atualizar o Valor
x4_max = 1; %Atualizar o Valor
x5_max = 1; %Atualizar o Valor
u_max = 12;

q1 = 1/((x1_max)^2);
q2 = 1/((x2_max)^2);
q3 = 1/((x3_max)^2);
q4 = 1/((x4_max)^2);
q5 = 1/((x5_max)^2);
r1 = 1/((u_max)^2);

Q = [q1 0 0 0 0; 0 q2 0 0 0; 0 0 q3 0 0; 0 0 0 q4 0; 0 0 0 0 q5];
R = r1;

[K,P,E] = lqr(A_exp,B_exp,Q,R)

Kr = [K(1) K(2) K(3) K(4)]
Ki = -K(5)

%Condição Inicial
x0 = [0 0.1 0 0]';

%% observador

% Observador de ordem Plena
Ta = 0.015 %(1/17.1223)/5;%Período de amostragem, considerado 5x mais rápido que a       
               %dinâmica do observador
Ad = (eye(4)+Ta*A);
Bd = Ta*B;
 
%Verificar Observabilidade par [Ad,C]
Od = [C; C*Ad; C*(Ad^2); C*(Ad^3)];
n2 = rank(Od)

autov = -3.4;
m = 3;

%Polos Observador
so1 = 1.1*m*autov; %
so2 = 1.2*m*autov;
so3 = 1.3*m*autov;
so4 = 1.4*m*autov;

%Projeto do observador
Pod = [exp(so1*Ta) exp(so2*Ta) exp(so3*Ta) exp(so4*Ta)];
kod = place(Ad',C',Pod);
ld = kod'