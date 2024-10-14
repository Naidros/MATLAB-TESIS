clear all
clear classes
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
clc
close all

%% CARGO EL ARCHIVO
load('C:\Users\adria\Desktop\DOCUMENTOS TESIS\Experimentos2Fugas\Lado_Derecho\DatosDosFugasEx4_HOUTBAJANDO_RUIDO.mat')

%% AÑADO YALMIP Y MOSEK
addpath(genpath('C:\Users\adria\Desktop\DOCUMENTOS TESIS\solvers\YALMIP-master'));
addpath(genpath('C:\Program Files\Mosek'));


%% BASE DE DATOS

QinFiltrado = out.QinReal.signals.values((1000:30000));
QoutFiltrado = out.QoutReal.signals.values((1000:30000));

HinFiltrado = out.HinReal.signals.values((1000:30000));
HoutFiltrado = out.HoutReal.signals.values((1000:30000));


iteraciones = length(QinFiltrado);

Qrest=(QinFiltrado-QoutFiltrado);


dt = 0.01;

tic

%% CARGO LOS PARÁMETROS DEL SIMULADOR
D=0.1048;
Ar=3.1416*(D/2)^2;
L=225;
b=721.6531;
f1=0.0276;
f2=0.0276;
f3=0.0276;
miu1=f1/(2*D*Ar);
miu2=f2/(2*D*Ar);
miu3=f3/(2*D*Ar);
g=9.81;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DESARROLLO PARA LA PRIMERA FUGA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CONDICIONES INICIALES
Q0 = mean(QinFiltrado(200:300));
H0 = mean(HinFiltrado(200:300));
z0 = 0.3*L;
H20=H0-f1*z0*Q0^2/(2*D*Ar^2*g);
% H20 = 14;

x1 = Q0;
x2 = H20;
x3 = Q0;
x4 = z0;
x5 = 0;
tic
x1iF1 = x1;
x2iF1 = x2;
x3iF1 = x3;
x4iF1 = x4;
x5iF1 = x5;
tic
%% METO ECUACIONES PARA LMI
P       = sdpvar(5,5);
Y       = sdpvar(5,5,'symmetric');
W1    = sdpvar(2,5);
W2    = sdpvar(2,5);
W3    = sdpvar(2,5);
W4    = sdpvar(2,5);
W5    = sdpvar(2,5);
W6    = sdpvar(2,5);
W7    = sdpvar(2,5);
W8    = sdpvar(2,5);
W9    = sdpvar(2,5);
W10    = sdpvar(2,5);
W11    = sdpvar(2,5);
W12    = sdpvar(2,5);
W13    = sdpvar(2,5);
W14    = sdpvar(2,5);
W15    = sdpvar(2,5);
W16    = sdpvar(2,5);

x1max=max(QinFiltrado);
x1min=min(QinFiltrado);

x2max=17;
x2min=5.5;


x3max=max(QoutFiltrado);
x3min=min(QoutFiltrado);

x4max=190;
x4min=20;

tic

Cd=[1 0 0 0 0;0 0 1 0 0];
Dd=0;

A1 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min)
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];

B1=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];

SYS1=ss(A1,B1,Cd,Dd);
SYSd1=c2d(SYS1,dt);
[adis1,bdis1,cdis1,ddis1]=ssdata(SYSd1);

%%           
A2 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min)
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0]; 
B2=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS2=ss(A2,B2,Cd,Dd);
SYSd2=c2d(SYS2,dt);
[adis2,bdis2,cdis2,ddis2]=ssdata(SYSd2);

%%
A3 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B3=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS3=ss(A3,B3,Cd,Dd);
SYSd3=c2d(SYS3,dt);
[adis3,bdis3,cdis3,ddis3]=ssdata(SYSd3);

%%              
A4 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0]; 
              
B4=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS4=ss(A4,B4,Cd,Dd);
SYSd4=c2d(SYS4,dt);
[adis4,bdis4,cdis4,ddis4]=ssdata(SYSd4);
%%
A5 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B5=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS5=ss(A5,B5,Cd,Dd);
SYSd5=c2d(SYS5,dt);
[adis5,bdis5,cdis5,ddis5]=ssdata(SYSd5);
%%
A6 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B6=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS6=ss(A6,B6,Cd,Dd);
SYSd6=c2d(SYS6,dt);
[adis6,bdis6,cdis6,ddis6]=ssdata(SYSd6);
%%              
A7 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B7=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS7=ss(A7,B7,Cd,Dd);
SYSd7=c2d(SYS7,dt);
[adis7,bdis7,cdis7,ddis7]=ssdata(SYSd7);
%%              
A8 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B8=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS8=ss(A8,B8,Cd,Dd);
SYSd8=c2d(SYS8,dt);
[adis8,bdis8,cdis8,ddis8]=ssdata(SYSd8);
%%              
A9 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B9=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS9=ss(A9,B9,Cd,Dd);
SYSd9=c2d(SYS9,dt);
[adis9,bdis9,cdis9,ddis9]=ssdata(SYSd9);
%%
A10 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B10=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS10=ss(A10,B10,Cd,Dd);
SYSd10=c2d(SYS10,dt);
[adis10,bdis10,cdis10,ddis10]=ssdata(SYSd10);
%%
A11 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B11=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS11=ss(A11,B11,Cd,Dd);
SYSd11=c2d(SYS11,dt);
[adis11,bdis11,cdis11,ddis11]=ssdata(SYSd11);
%%
A12 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B12=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS12=ss(A12,B12,Cd,Dd);
SYSd12=c2d(SYS12,dt);
[adis12,bdis12,cdis12,ddis12]=ssdata(SYSd12);
%%
A13 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
B13=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS13=ss(A13,B13,Cd,Dd);
SYSd13=c2d(SYS13,dt);
[adis13,bdis13,cdis13,ddis13]=ssdata(SYSd13);
%%
A14 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B14=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS14=ss(A14,B14,Cd,Dd);
SYSd14=c2d(SYS14,dt);
[adis14,bdis14,cdis14,ddis14]=ssdata(SYSd14);
%%
A15 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B15=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS15=ss(A15,B15,Cd,Dd);
SYSd15=c2d(SYS15,dt);
[adis15,bdis15,cdis15,ddis15]=ssdata(SYSd15);
%%
A16 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B16=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS16=ss(A16,B16,Cd,Dd);
SYSd16=c2d(SYS16,dt);
[adis16,bdis16,cdis16,ddis16]=ssdata(SYSd16);
              
Ad1 = adis1;
Ad2=adis2;
Ad3=adis3;
Ad4=adis4;
Ad5=adis5;
Ad6=adis6;
Ad7=adis7;
Ad8=adis8;
Ad9=adis9;
Ad10=adis10;
Ad11=adis11;
Ad12=adis12;
Ad13=adis13;
Ad14=adis14;
Ad15=adis15;
Ad16=adis16;
C= [1 0 0 0 0; 0 0 1 0 0]

% Qk=[1*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-2 0 0;0 0 0 2*10^3 0;0 0 0 0 1*10^-7];

%Qk=[1*10^-5 0 0 0 0 0 0;0 1*10^-2 0 0 0 0 0;0 0 1*10^-5 0 0 0 0;0 0 0 1*10^-2 0 0 0;0 0 0 0 1*10^-5 0 0; 0 0 0 0 0 2*10^3 0; 0 0 0 0 0 0 1*10^-7];
Qk=[1*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-6 0 0;0 0 0 20000 0;0 0 0 0 1*10^-7];

Rk=[1*10^-5 0;0 1*10^-5]; %

R_obs=Rk;
H = Qk^0.5;

sdpvar gamma_LQR_obs;

Trace_cond  = [gamma_LQR_obs*eye(5) eye(5); eye(5) Y];
%Model = [Trace_cond >= t*eye(2)];

        %Pole_1 = [kron(L_reg,Y)+kron(M_reg,(Ad1'*Y-C'*W1)')+kron(M_reg',(Ad1'*Y-C'*W1))];
H1     = [      -Y       ,  Y*Ad1-W1'*C  ,   Y*H'     ,      W1';
                  Ad1'*Y-C'*W1   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W1        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                    %Condition set
                   
H2     = [      -Y       ,  Y*Ad2-W2'*C  ,   Y*H'     ,      W2';
                  Ad2'*Y-C'*W2   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W2        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H3     = [      -Y       ,  Y*Ad3-W3'*C  ,   Y*H'     ,      W3';
                  Ad3'*Y-C'*W3   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W3        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H4     = [      -Y       ,  Y*Ad4-W4'*C  ,   Y*H'     ,      W4';
                  Ad4'*Y-C'*W4   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W4        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H5     = [      -Y       ,  Y*Ad5-W5'*C  ,   Y*H'     ,      W5';
                  Ad5'*Y-C'*W5   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W5        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H6     = [      -Y       ,  Y*Ad6-W6'*C  ,   Y*H'     ,      W6';
                  Ad6'*Y-C'*W6   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W6        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H7     = [      -Y       ,  Y*Ad7-W7'*C  ,   Y*H'     ,      W7';
                  Ad7'*Y-C'*W7   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W7        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H8     = [      -Y       ,  Y*Ad8-W8'*C  ,   Y*H'     ,      W8';
                  Ad8'*Y-C'*W8   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W8        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H9     = [      -Y       ,  Y*Ad9-W9'*C  ,   Y*H'     ,      W9';
                  Ad9'*Y-C'*W9   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W9        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H10     = [      -Y       ,  Y*Ad10-W10'*C  ,   Y*H'     ,      W10';
                  Ad10'*Y-C'*W10   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W10        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H11     = [      -Y       ,  Y*Ad11-W11'*C  ,   Y*H'     ,      W11';
                  Ad11'*Y-C'*W11   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W11        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H12     = [      -Y       ,  Y*Ad12-W12'*C  ,   Y*H'     ,      W12';
                  Ad12'*Y-C'*W12   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W12        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H13     = [      -Y       ,  Y*Ad13-W13'*C  ,   Y*H'     ,      W13';
                  Ad13'*Y-C'*W13   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W13        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H14     = [      -Y       ,  Y*Ad14-W14'*C  ,   Y*H'     ,      W14';
                  Ad14'*Y-C'*W14   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W14        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H15     = [      -Y       ,  Y*Ad15-W15'*C  ,   Y*H'     ,      W15';
                  Ad15'*Y-C'*W15   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W15        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H16     = [      -Y       ,  Y*Ad16-W16'*C  ,   Y*H'     ,      W16';
                  Ad16'*Y-C'*W16   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W16        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
        F = [Y>=0];
        F = F + [Trace_cond <= 0];
        F = F + [H1<=0];


% Configurar opciones del solucionador
ops = sdpsettings('solver', 'mosek', 'verbose', 0);

% Optimizar
       % [opt_info]  = optimize([[H1<=0]+[H2<=0]+[H3<=0]+[H4<=0]+[H5<=0]+[H6<=0]+[H7<=0]+[H8<=0]+[H9<=0]+[H10<=0]+[H11<=0]+[H12<=0]+[H13<=0]+[H14<=0]+[H15<=0]+[H16<=0]+[Y>=0]+[Trace_cond>=0]],[gamma_LQR_obs],ops)
        [opt_info]  = optimize([[H1<=0]+[H2<=0]+[H3<=0]+[H4<=0]+[H5<=0]+[H6<=0]+[H7<=0]+[H8<=0]+[H9<=0]+[H10<=0]+[H11<=0]+[H12<=0]+[H13<=0]+[H14<=0]+[H15<=0]+[H16<=0]+[Y>=0]+[Trace_cond>=0],[gamma_LQR_obs>=0]],[],ops)

Y_sol  = value(Y);%  Y solution -> WARNING: LQR solves for -Y!!

        W1_sol = value(W1); % W1 solution
        W2_sol = value(W2); % W2 solution
        W3_sol = value(W3); % W3 solution
        W4_sol = value(W4); % W4 solution
        
        W5_sol = value(W5); % W5 solution
        W6_sol = value(W6); % W6 solution
        W7_sol = value(W7); % W7 solution
        W8_sol = value(W8); % W8 solution
        
        W9_sol = value(W9); % W9 solution
        W10_sol = value(W10); % W10 solution
        W11_sol = value(W11); % W11 solution
        W12_sol = value(W12); % W12 solution
        
        W13_sol = value(W13); % W13 solution
        W14_sol = value(W14); % W14 solution
        W15_sol = value(W15); % W15 solution
        W16_sol = value(W16); % W16 solution


    
        L1F1 = inv(Y_sol)*W1_sol'
        tic
        L2F1 = inv(Y_sol)*W2_sol'
        L3F1 = inv(Y_sol)*W3_sol'
        L4F1 = inv(Y_sol)*W4_sol'
        
        L5F1 = inv(Y_sol)*W5_sol'
        L6F1 = inv(Y_sol)*W6_sol'
        L7F1 = inv(Y_sol)*W7_sol'
        L8F1 = inv(Y_sol)*W8_sol'
        
        L9F1 = inv(Y_sol)*W9_sol'
        L10F1 = inv(Y_sol)*W10_sol'
        L11F1 = inv(Y_sol)*W11_sol'
        L12F1 = inv(Y_sol)*W12_sol'
        
        L13F1 = inv(Y_sol)*W13_sol'
        L14F1 = inv(Y_sol)*W14_sol'
        L15F1 = inv(Y_sol)*W15_sol'
        L16F1 = inv(Y_sol)*W16_sol'

n=1
enclave=0;
while (n <= iteraciones)
    % Filtrar entradas en cada iteración
    Qin = QinFiltrado(n, 1);
    Qout = QoutFiltrado(n, 1);
    Hin = HinFiltrado(n, 1);
    Hout = HoutFiltrado(n, 1);


    if (Qin - Qout) > 0.5e-03 && enclave ==0
        enclave = 1;
    end

    if enclave ==1
            Q1max = max(QinFiltrado);
            Q1min = min(QinFiltrado);
            H2max = 17;
            H2min = 5;
            Q2max = max(QoutFiltrado);
            Q2min = min(QoutFiltrado);
            z1max = 225;
            z1min = 15;

            % Normalización (escalado) de las variables
            Q1Sch = (Q1max - x1) / (Q1max - Q1min);
            H2Sch = (H2max - x2) / (H2max - H2min);
            Q2Sch = (Q2max - x3) / (Q2max - Q2min);
            z1Sch = (z1max - x4) / (z1max - z1min);
            tic
            % Cálculo de estados del sistema
            Q1 = x1;
            H2 = x2;
            Q2 = x3;
            z1 = x4;
            lam = x5;
            Hk = [1 0 0 0 0; 0 0 1 0 0];

            % Parámetros y cálculos para dinámica del sistema
            A = Ar;
            promediolongitud = L;
            Q11=(-g*A*(H2-Hin)/(z1)-f1*Q1^2/(2*D*A))*dt+Q1;
            H21=(-b^2*(Q2-Q1+lam*sqrt(H2))/(g*A*z1))*dt+H2; 
            Q21=(-g*A*(Hout-H2)/(promediolongitud-z1)-f2*Q2^2/(2*D*A))*dt+Q2;
            tic
            %% Matrices para la dinámica del sistema
            Aact1 = [-miu1 * Q11, (-g * A / z1) * 0.5, 0, (-g * A * H21 / z1^2) * 0.5, 0;
                     b^2 / (g * A * z1), 0, -b^2 / (g * A * z1), 0, (-b^2 * sqrt(abs(H21))) / (g * A * z1);
                     0, g * A / (promediolongitud - z1), -miu2 * Q21, 0, 0;
                     0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0];
            xact1 = [Q11, H21, Q21, z1, lam]' * (dt / 2);
            Bact1 = [(g * A) / z1, 0;
                     0, 0;
                     0, (-g * A) / (promediolongitud - z1);
                     0, 0;
                     0, 0];
            u = [Hin * (dt / 2); Hout * (dt / 2)];
            Ax1 = Aact1 * xact1;
            Bu1 = Bact1 * u;

            %% Segunda parte de la dinámica del sistema
            Aact2 = [-miu1 * Q1, (-g * A / z1) * 0.5, 0, (-g * A * H2 / z1^2) * 0.5, 0;
                     b^2 / (g * A * z1), 0, -b^2 / (g * A * z1), 0, (-b^2 * sqrt(abs(H2))) / (g * A * z1);
                     0, g * A / (promediolongitud - z1), -miu2 * Q2, 0, 0;
                     0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0];
            xact2 = [Q1, H2, Q2, z1, lam]' * (dt / 2);
            Bact2 = [(g * A) / z1, 0;
                     0, 0;
                     0, (-g * A) / (promediolongitud - z1);
                     0, 0;
                     0, 0];
            Ax2 = Aact2 * xact2;
            Bu2 = Bact2 * u;

            %% Actualización de estados
            xx = [Q1; H2; Q2; z1; lam];
            dXmat = (Ax1 + Ax2) + (Bu1 + Bu2) + xx;
            dQ1 = dXmat(1);
            dH2 = dXmat(2);
            dQ2 = dXmat(3);
            dXm = [dQ1; dH2; dQ2; z1; lam];



            %% Cálculo de los coeficientes de pertenencia `mu`
            Q1SchMax= (1 - Q1Sch);
            Q1SchMin= Q1Sch;
            H2SchMax= (1 - H2Sch);
            H2SchMin= H2Sch;
            Q2SchMax= (1 - Q2Sch);
            Q2SchMin= Q2Sch;
            z1SchMax= (1-z1Sch);
            z1SchMin = z1Sch;

            %% FUNCIÓN PARA OBTENER LA COMBINATORIA DE MU'S
mu_1(n)= Q1Sch * H2Sch * Q2Sch * z1Sch; %todos MIN
mu_2(n)= (1-Q1Sch) * H2Sch * Q2Sch * z1Sch;
mu_3(n)= Q1Sch * (1-H2Sch) * Q2Sch * z1Sch;
mu_4(n)= Q1Sch * H2Sch * (1-Q2Sch) * z1Sch;

mu_5(n)= Q1Sch * H2Sch * Q2Sch * (1-z1Sch);
mu_6(n)= (1-Q1Sch) * (1-H2Sch) * (1-Q2Sch) * (1-z1Sch); %todos MAX
mu_7(n)= Q1Sch * (1-H2Sch) * (1-Q2Sch) * (1-z1Sch);
mu_8(n)= (1-Q1Sch) * H2Sch * (1-Q2Sch) * (1-z1Sch);

mu_9(n)= (1-Q1Sch) * (1-H2Sch) * Q2Sch * (1-z1Sch);
mu_10(n)= (1-Q1Sch) * (1-H2Sch) * (1-Q2Sch) * z1Sch;
mu_11(n)= (1-Q1Sch) * (1-H2Sch) * Q2Sch * z1Sch;
mu_12(n)= (1-Q1Sch) * H2Sch * (1-Q2Sch) * z1Sch;
 
mu_13(n)= (1-Q1Sch) * H2Sch * Q2Sch * (1-z1Sch);
mu_14(n)= Q1Sch * H2Sch * (1-Q2Sch) * (1-z1Sch);
mu_15(n)= Q1Sch * (1-H2Sch) * Q2Sch * (1-z1Sch);
mu_16(n)= Q1Sch * (1-H2Sch) * (1-Q2Sch) * z1Sch;           

            %% Interpolación y cálculo de L_interp (corregido)
%            L_interp = 0;
musL1(n)=mu_1(n)+mu_2(n)+mu_3(n)+mu_4(n)+mu_5(n)+mu_6(n)+mu_7(n)+mu_8(n)+mu_9(n)+mu_10(n)+mu_11(n)+mu_12(n)+mu_13(n)+mu_14(n)+mu_15(n)+mu_16(n);

 L_interp            = mu_1(n)*L1F1 + mu_2(n)*L2F1 + mu_3(n)*L3F1 + mu_4(n)*L4F1...
                       + mu_5(n)*L5F1 + mu_6(n)*L6F1 + mu_7(n)*L7F1 + mu_8(n)*L8F1...
                       + mu_9(n)*L9F1 + mu_10(n)*L10F1 + mu_11(n)*L11F1 + mu_12(n)*L12F1...
                       + mu_13(n)*L13F1 + mu_14(n)*L14F1 + mu_15(n)*L15F1 + mu_16(n)*L16F1;
                   
          Kk=L_interp;
            dX = dXm + Kk * ([Qin; Qout] - [dQ1; dQ2]);

            %% Actualización de las variables dinámicas
            x1 = dX(1);
            x2 = dX(2);
            x3 = dX(3);
            x4 = dX(4);
            x5 = dX(5);

    else
        x1 = x1iF1;
        x2 = x2iF1;
        x3 = x3iF1;
        x4 = x4iF1;
        x5 = x5iF1;

        % Inicialización de vectores
        dX = [x1; x2; x3; x4; x5];
        mus = 0;
        mu_vals = zeros(16, 1);  % Inicializamos todos los `mu_i` a 0
        tic
    end

    Q1Fuga1Observador(n)=dX(1);
    H2Fuga1Observador(n)=dX(2);
    Q2Fuga1Observador(n)=dX(3);
    z1Fuga1Observador(n)=dX(4);
    lamda1Fuga1Observador(n)=dX(5);

    if (Qin-Qout>=7.8e-04 && n>=20000)
        n=50000;
    end
    n=n+1
end
tic
% Vector original
z1Fuga1Observador = z1Fuga1Observador(:)'; % Asegurarte de que el vector es fila (1x26501)
H2Fuga1Observador = H2Fuga1Observador(:)';
lamda1Fuga1Observador = lamda1Fuga1Observador(:)';
muL1 = musL1(:)'
% Obtener el último valor del vector original
ultimo_valor = z1Fuga1Observador(end);
ultimo_valorH2 = H2Fuga1Observador(end);
ultimo_valorLAM1 = lamda1Fuga1Observador(end);
ultimo_valorMUL1 = musL1(end);
% Crear un vector con los valores adicionales (del 26502 al 49000)
valores_adicionales = repmat(ultimo_valor, 1, 49002 - 29001);
valores_adicionalesH2 = repmat(ultimo_valorH2, 1, 49002 - 29001);
valores_adicionalesLAM1 = repmat(ultimo_valorLAM1, 1, 49002 - 29001);
valores_adicionalesMUSL1 = repmat(ultimo_valorMUL1, 1, 49002 - 29001);


% Concatenar el vector original con los valores adicionales
z1Fuga1Observador_extendido = [z1Fuga1Observador, valores_adicionales];
H2Fuga1Observador_extendido = [H2Fuga1Observador, valores_adicionalesH2];
LAM1Fuga1Observador_extendido = [lamda1Fuga1Observador, valores_adicionalesLAM1];
MUSL1Fuga1Observador_extendido = [muL1, valores_adicionalesMUSL1];


z1 = z1Fuga1Observador(1,end)
lam1 = lamda1Fuga1Observador(1,end)
tic


%% BASE DE DATOS
QinFiltrado = out.QinReal.signals.values((1000:end));
QoutFiltrado = out.QoutReal.signals.values((1000:end));

HinFiltrado = out.HinReal.signals.values((1000:end));
HoutFiltrado = out.HoutReal.signals.values((1000:end));

%%  Agregar ruido
t=(out.tout(1000:end));

tic

iteraciones = length(QinFiltrado);

Qrest=(QinFiltrado-QoutFiltrado);

dt = 0.01;
tic
%% CONDICIONES INICIALES

Q0 = mean(QinFiltrado(20000:25000),1);
H0 = mean(HinFiltrado(200:300),1);
Qo0 = mean(QoutFiltrado(2000:2500),1);
z0 = 108;
H20o=H0-f1*z0*Q0^2/(2*D*Ar^2*g);

x1 = Q0;
x2 = 15;
x3 = 0.01693;
x4 = 12;
x5 = 0.01693;
x6 = z0;
x7 = 0;
tic

x1i = x1;
x2i = x2;
x3i = x3;
x4i = x4;
x5i = x5;
x6i = x6;
x7i = x7;

%% METO ECUACIONES PARA LMI
P       = sdpvar(7,7);
Y       = sdpvar(7,7,'symmetric');
W1  = sdpvar(2,7);
W2  = sdpvar(2,7);
W3  = sdpvar(2,7);
W4  = sdpvar(2,7);
W5  = sdpvar(2,7);
W6  = sdpvar(2,7);
W7  = sdpvar(2,7);
W8  = sdpvar(2,7);
W9  = sdpvar(2,7);
W10 = sdpvar(2,7);
W11 = sdpvar(2,7);
W12 = sdpvar(2,7);
W13 = sdpvar(2,7);
W14 = sdpvar(2,7);
W15 = sdpvar(2,7);
W16 = sdpvar(2,7);
W17 = sdpvar(2,7);
W18 = sdpvar(2,7);
W19 = sdpvar(2,7);
W20 = sdpvar(2,7);
W21 = sdpvar(2,7);
W22 = sdpvar(2,7);
W23 = sdpvar(2,7);
W24 = sdpvar(2,7);
W25 = sdpvar(2,7);
W26 = sdpvar(2,7);
W27 = sdpvar(2,7);
W28 = sdpvar(2,7);
W29 = sdpvar(2,7);
W30 = sdpvar(2,7);
W31 = sdpvar(2,7);
W32 = sdpvar(2,7);
W33 = sdpvar(2,7);
W34 = sdpvar(2,7);
W35 = sdpvar(2,7);
W36 = sdpvar(2,7);
W37 = sdpvar(2,7);
W38 = sdpvar(2,7);
W39 = sdpvar(2,7);
W40 = sdpvar(2,7);
W41 = sdpvar(2,7);
W42 = sdpvar(2,7);
W43 = sdpvar(2,7);
W44 = sdpvar(2,7);
W45 = sdpvar(2,7);
W46 = sdpvar(2,7);
W47 = sdpvar(2,7);
W48 = sdpvar(2,7);
W49 = sdpvar(2,7);
W50 = sdpvar(2,7);
W51 = sdpvar(2,7);
W52 = sdpvar(2,7);
W53 = sdpvar(2,7);
W54 = sdpvar(2,7);
W55 = sdpvar(2,7);
W56 = sdpvar(2,7);
W57 = sdpvar(2,7);
W58 = sdpvar(2,7);
W59 = sdpvar(2,7);
W60 = sdpvar(2,7);
W61 = sdpvar(2,7);
W62 = sdpvar(2,7);
W63 = sdpvar(2,7);
W64 = sdpvar(2,7);


%% EN ESTA PRUEBA LO ESTOY CALANDO CON EL QINFILTRADO DE TODO EL EXPEIRMENTO
%% **** EXPERIMENTAR SOLO CON LA MUESTRA DE 20000:END

%% FUNCIONA PARA 2000 Y 20000
x1max=max(QinFiltrado(2000:end));
x1min=min(QinFiltrado(2000:end));

x2max=14;
x2min=10;

tic
%% ESTAS TAMBIÉN FUNCIONAN, PERO EN TEORÍA DESCONOZCO EL FLUJO EN Q2 (SACADO DE SIMULINK)
%%x3max = 0.020
%%x3min = 0.015

x3max=x1max-(lam1*sqrt(x2max));
x3min= x1min-(lam1*sqrt(x2min));

x4max=9.5;
x4min=2.5;

x5max = max(QoutFiltrado(2000:end));
x5min = min(QoutFiltrado(2000:end));

x6max = 160;
x6min = 140;

Cd=[1 0 0 0 0 0 0;0 0 0 0 1 0 0];
Dd=0;

Combinatoria = matrizVerdad(6)

[adis,bdis,cdis,ddis,matricesA,matricesB] = evaluarF2(z1,lam1,x1max,x1min,x2max,x2min,x3max,x3min,x4max,x4min,x5max,x5min,x6max,x6min,Combinatoria,6,Cd,Dd,dt)

Ad1 = adis{1,1}; 
Ad2 = adis{2,1}; 
Ad3 = adis{3,1}; 
Ad4 = adis{4,1}; 
Ad5 = adis{5,1}; 
Ad6 = adis{6,1}; 
Ad7 = adis{7,1}; 
Ad8 = adis{8,1}; 
Ad9 = adis{9,1}; 
Ad10 = adis{10,1}; 
Ad11 = adis{11,1}; 
Ad12 = adis{12,1}; 
Ad13 = adis{13,1}; 
Ad14 = adis{14,1}; 
Ad15 = adis{15,1}; 
Ad16 = adis{16,1}; 
Ad17 = adis{17,1}; 
Ad18 = adis{18,1}; 
Ad19 = adis{19,1}; 
Ad20 = adis{20,1}; 
Ad21 = adis{21,1}; 
Ad22 = adis{22,1}; 
Ad23 = adis{23,1}; 
Ad24 = adis{24,1}; 
Ad25 = adis{25,1}; 
Ad26 = adis{26,1}; 
Ad27 = adis{27,1}; 
Ad28 = adis{28,1}; 
Ad29 = adis{29,1}; 
Ad30 = adis{30,1}; 
Ad31 = adis{31,1}; 
Ad32 = adis{32,1}; 
Ad33 = adis{33,1}; 
Ad34 = adis{34,1}; 
Ad35 = adis{35,1}; 
Ad36 = adis{36,1}; 
Ad37 = adis{37,1}; 
Ad38 = adis{38,1}; 
Ad39 = adis{39,1}; 
Ad40 = adis{40,1}; 
Ad41 = adis{41,1}; 
Ad42 = adis{42,1}; 
Ad43 = adis{43,1}; 
Ad44 = adis{44,1}; 
Ad45 = adis{45,1}; 
Ad46 = adis{46,1}; 
Ad47 = adis{47,1}; 
Ad48 = adis{48,1}; 
Ad49 = adis{49,1}; 
Ad50 = adis{50,1}; 
Ad51 = adis{51,1}; 
Ad52 = adis{52,1}; 
Ad53 = adis{53,1}; 
Ad54 = adis{54,1}; 
Ad55 = adis{55,1}; 
Ad56 = adis{56,1}; 
Ad57 = adis{57,1}; 
Ad58 = adis{58,1}; 
Ad59 = adis{59,1}; 
Ad60 = adis{60,1}; 
Ad61 = adis{61,1}; 
Ad62 = adis{62,1}; 
Ad63 = adis{63,1}; 
Ad64 = adis{64,1};
C=Cd;

% Qk=[1*10^-6 0 0 0 0 0 0 ; 0 0.1 0 0 0 0 0 ; 0 0 1*10^-6 0 0 0 0 ; 0 0 0 1*10^-3 0 0 0; 0 0 0 0 1*10^-7 0 0; 0 0 0 0 0 20000 0; 0 0 0 0 0 0 1*10^-5];
% Rk=[1*10^-8 0;0 1*10^-8]; %

%Qk=[1*10^-5 0 0 0 0 0 0;0 0.01 0 0 0 0 0;0 0 1*10^-5 0 0 0 0;0 0 0 2*10^-2 0 0 0;0 0 0 0 1*10^-7 0 0; 0 0 0 0 0 1*10^-1 0; 0 0 0 0 0 0 5*10^-10];

%% ESTA Qk FUNCIONA BIEN PARA DOS FUGAS
%Qk=[1*10^-5 0 0 0 0 0 0;0 0.01 0 0 0 0 0;0 0 1*10^-5 0 0 0 0;0 0 0 2*10^-2 0 0 0;0 0 0 0 1*10^-7 0 0; 0 0 0 0 0 2*10^3 0; 0 0 0 0 0 0 5*10^-10];

%% QK EXPERIMENTAL
Qk=[1*10^-6 0 0 0 0 0 0;0 0.01 0 0 0 0 0;0 0 1*10^-6 0 0 0 0;0 0 0 2*10^-3 0 0 0;0 0 0 0 1*10^-7 0 0; 0 0 0 0 0 2*10^3 0; 0 0 0 0 0 0 5*10^-10];



Qk=[1*10^-8 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1*10^-8 0 0 0 0;0 0 0 2*10^-6 0 0 0;0 0 0 0 1*10^-10 0 0; 0 0 0 0 0 8*10^3 0; 0 0 0 0 0 0 5*10^-12];


Qk=[1*10^-5 0 0 0 0 0 0;0 1*10^-2 0 0 0 0 0;0 0 1*10^-5 0 0 0 0;0 0 0 1*10^-2 0 0 0;0 0 0 0 1*10^-5 0 0; 0 0 0 0 0 2*10^3 0; 0 0 0 0 0 0 1*10^-7];

Rk=[1*10^-5 0;0 1*10^-5]; %funciona bien para muestro de 100 Hz

R_obs=Rk;
H = Qk.^0.5;
tic
sdpvar gamma_LQR_obs;

Trace_cond  = [gamma_LQR_obs*eye(7) eye(7); eye(7) Y];

H1  = [ -Y , Y*Ad1-W1'*C , Y*H' , W1';
        Ad1'*Y-C'*W1 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W1 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H2  = [ -Y , Y*Ad2-W2'*C , Y*H' , W2';
        Ad2'*Y-C'*W2 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W2 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H3  = [ -Y , Y*Ad3-W3'*C , Y*H' , W3';
        Ad3'*Y-C'*W3 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W3 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H4  = [ -Y , Y*Ad4-W4'*C , Y*H' , W4';
        Ad4'*Y-C'*W4 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W4 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H5  = [ -Y , Y*Ad5-W5'*C , Y*H' , W5';
        Ad5'*Y-C'*W5 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W5 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H6  = [ -Y , Y*Ad6-W6'*C , Y*H' , W6';
        Ad6'*Y-C'*W6 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W6 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H7  = [ -Y , Y*Ad7-W7'*C , Y*H' , W7';
        Ad7'*Y-C'*W7 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W7 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H8  = [ -Y , Y*Ad8-W8'*C , Y*H' , W8';
        Ad8'*Y-C'*W8 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W8 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H9  = [ -Y , Y*Ad9-W9'*C , Y*H' , W9';
        Ad9'*Y-C'*W9 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W9 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H10 = [ -Y , Y*Ad10-W10'*C , Y*H' , W10';
        Ad10'*Y-C'*W10 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W10 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H11 = [ -Y , Y*Ad11-W11'*C , Y*H' , W11';
        Ad11'*Y-C'*W11 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W11 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H12 = [ -Y , Y*Ad12-W12'*C , Y*H' , W12';
        Ad12'*Y-C'*W12 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W12 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H13 = [ -Y , Y*Ad13-W13'*C , Y*H' , W13';
        Ad13'*Y-C'*W13 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W13 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H14 = [ -Y , Y*Ad14-W14'*C , Y*H' , W14';
        Ad14'*Y-C'*W14 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W14 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H15 = [ -Y , Y*Ad15-W15'*C , Y*H' , W15';
        Ad15'*Y-C'*W15 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W15 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H16 = [ -Y , Y*Ad16-W16'*C , Y*H' , W16';
        Ad16'*Y-C'*W16 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W16 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H17 = [ -Y , Y*Ad17-W17'*C , Y*H' , W17';
        Ad17'*Y-C'*W17 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W17 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H18 = [ -Y , Y*Ad18-W18'*C , Y*H' , W18';
        Ad18'*Y-C'*W18 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W18 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H19 = [ -Y , Y*Ad19-W19'*C , Y*H' , W19';
        Ad19'*Y-C'*W19 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W19 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H20 = [ -Y , Y*Ad20-W20'*C , Y*H' , W20';
        Ad20'*Y-C'*W20 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W20 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H21 = [ -Y , Y*Ad21-W21'*C , Y*H' , W21';
        Ad21'*Y-C'*W21 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W21 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H22 = [ -Y , Y*Ad22-W22'*C , Y*H' , W22';
        Ad22'*Y-C'*W22 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W22 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H23 = [ -Y , Y*Ad23-W23'*C , Y*H' , W23';
        Ad23'*Y-C'*W23 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W23 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H24 = [ -Y , Y*Ad24-W24'*C , Y*H' , W24';
        Ad24'*Y-C'*W24 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W24 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H25 = [ -Y , Y*Ad25-W25'*C , Y*H' , W25';
        Ad25'*Y-C'*W25 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W25 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H26 = [ -Y , Y*Ad26-W26'*C , Y*H' , W26';
        Ad26'*Y-C'*W26 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W26 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H27 = [ -Y , Y*Ad27-W27'*C , Y*H' , W27';
        Ad27'*Y-C'*W27 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W27 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H28 = [ -Y , Y*Ad28-W28'*C , Y*H' , W28';
        Ad28'*Y-C'*W28 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W28 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H29 = [ -Y , Y*Ad29-W29'*C , Y*H' , W29';
        Ad29'*Y-C'*W29 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W29 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H30 = [ -Y , Y*Ad30-W30'*C , Y*H' , W30';
        Ad30'*Y-C'*W30 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W30 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H31 = [ -Y , Y*Ad31-W31'*C , Y*H' , W31';
        Ad31'*Y-C'*W31 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W31 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H32 = [ -Y , Y*Ad32-W32'*C , Y*H' , W32';
        Ad32'*Y-C'*W32 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W32 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H33 = [ -Y , Y*Ad33-W33'*C , Y*H' , W33';
        Ad33'*Y-C'*W33 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W33 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H34 = [ -Y , Y*Ad34-W34'*C , Y*H' , W34';
        Ad34'*Y-C'*W34 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W34 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H35 = [ -Y , Y*Ad35-W35'*C , Y*H' , W35';
        Ad35'*Y-C'*W35 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W35 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H36 = [ -Y , Y*Ad36-W36'*C , Y*H' , W36';
        Ad36'*Y-C'*W36 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W36 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H37 = [ -Y , Y*Ad37-W37'*C , Y*H' , W37';
        Ad37'*Y-C'*W37 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W37 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H38 = [ -Y , Y*Ad38-W38'*C , Y*H' , W38';
        Ad38'*Y-C'*W38 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W38 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H39 = [ -Y , Y*Ad39-W39'*C , Y*H' , W39';
        Ad39'*Y-C'*W39 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W39 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H40 = [ -Y , Y*Ad40-W40'*C , Y*H' , W40';
        Ad40'*Y-C'*W40 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W40 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];
H41 = [ -Y , Y*Ad41-W41'*C , Y*H' , W41';
        Ad41'*Y-C'*W41 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W41 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H42 = [ -Y , Y*Ad42-W42'*C , Y*H' , W42';
        Ad42'*Y-C'*W42 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W42 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H43 = [ -Y , Y*Ad43-W43'*C , Y*H' , W43';
        Ad43'*Y-C'*W43 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W43 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H44 = [ -Y , Y*Ad44-W44'*C , Y*H' , W44';
        Ad44'*Y-C'*W44 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W44 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H45 = [ -Y , Y*Ad45-W45'*C , Y*H' , W45';
        Ad45'*Y-C'*W45 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W45 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H46 = [ -Y , Y*Ad46-W46'*C , Y*H' , W46';
        Ad46'*Y-C'*W46 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W46 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H47 = [ -Y , Y*Ad47-W47'*C , Y*H' , W47';
        Ad47'*Y-C'*W47 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W47 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H48 = [ -Y , Y*Ad48-W48'*C , Y*H' , W48';
        Ad48'*Y-C'*W48 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W48 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H49 = [ -Y , Y*Ad49-W49'*C , Y*H' , W49';
        Ad49'*Y-C'*W49 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W49 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H50 = [ -Y , Y*Ad50-W50'*C , Y*H' , W50';
        Ad50'*Y-C'*W50 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W50 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H51 = [ -Y , Y*Ad51-W51'*C , Y*H' , W51';
        Ad51'*Y-C'*W51 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W51 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H52 = [ -Y , Y*Ad52-W52'*C , Y*H' , W52';
        Ad52'*Y-C'*W52 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W52 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H53 = [ -Y , Y*Ad53-W53'*C , Y*H' , W53';
        Ad53'*Y-C'*W53 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W53 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H54 = [ -Y , Y*Ad54-W54'*C , Y*H' , W54';
        Ad54'*Y-C'*W54 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W54 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H55 = [ -Y , Y*Ad55-W55'*C , Y*H' , W55';
        Ad55'*Y-C'*W55 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W55 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H56 = [ -Y , Y*Ad56-W56'*C , Y*H' , W56';
        Ad56'*Y-C'*W56 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W56 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H57 = [ -Y , Y*Ad57-W57'*C , Y*H' , W57';
        Ad57'*Y-C'*W57 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W57 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H58 = [ -Y , Y*Ad58-W58'*C , Y*H' , W58';
        Ad58'*Y-C'*W58 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W58 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H59 = [ -Y , Y*Ad59-W59'*C , Y*H' , W59';
        Ad59'*Y-C'*W59 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W59 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H60 = [ -Y , Y*Ad60-W60'*C , Y*H' , W60';
        Ad60'*Y-C'*W60 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W60 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H61 = [ -Y , Y*Ad61-W61'*C , Y*H' , W61';
        Ad61'*Y-C'*W61 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W61 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H62 = [ -Y , Y*Ad62-W62'*C , Y*H' , W62';
        Ad62'*Y-C'*W62 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W62 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H63 = [ -Y , Y*Ad63-W63'*C , Y*H' , W63';
        Ad63'*Y-C'*W63 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W63 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H64 = [ -Y , Y*Ad64-W64'*C , Y*H' , W64';
        Ad64'*Y-C'*W64 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W64 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];



tic
ops= sdpsettings('solver','mosek','verbose',0);
%F=[Y>=0] + [Trace_cond>=0]+[H1<=0] + [H2<=0] + [H3<=0] + [H4<=0] + [H5<=0] + [H6<=0] + [H7<=0] + [H8<=0] + [H9<=0] + [H10<=0] + [H11<=0] + [H12<=0] + [H13<=0] + [H14<=0] + [H15<=0] + [H16<=0] + ...
                       % [H17<=0] + [H18<=0] + [H19<=0] + [H20<=0] + [H21<=0] + [H22<=0] + [H23<=0] + [H24<=0] + [H25<=0] + [H26<=0] + [H27<=0] + [H28<=0] + [H29<=0] + [H30<=0] + [H31<=0] + [H32<=0] + ...
                       % [H33<=0] + [H34<=0] + [H35<=0] + [H36<=0] + [H37<=0] + [H38<=0] + [H39<=0] + [H40<=0] + [H41<=0] + [H42<=0] + [H43<=0] + [H44<=0] + [H45<=0] + [H46<=0] + [H47<=0] + [H48<=0] + ...
                       % [H49<=0] + [H50<=0] + [H51<=0] + [H52<=0] + [H53<=0] + [H54<=0] + [H55<=0] + [H56<=0] + [H57<=0] + [H58<=0] + [H59<=0] + [H60<=0] + [H61<=0] + [H62<=0] + [H63<=0] + [H64<=0];


%[opt_info] = optimize( F,gamma_LQR_obs,ops);

% [opt_info]=optimize([Y>=0] + [Trace_cond>=0]+[H1<=0] + [H2<=0] + [H3<=0] + [H4<=0] + [H5<=0] + [H6<=0] + [H7<=0] + [H8<=0] + [H9<=0] + [H10<=0] + [H11<=0] + [H12<=0] + [H13<=0] + [H14<=0] + [H15<=0] + [H16<=0] + ...
%                        [H17<=0] + [H18<=0] + [H19<=0] + [H20<=0] + [H21<=0] + [H22<=0] + [H23<=0] + [H24<=0] + [H25<=0] + [H26<=0] + [H27<=0] + [H28<=0] + [H29<=0] + [H30<=0] + [H31<=0] + [H32<=0] + ...
%                        [H33<=0] + [H34<=0] + [H35<=0] + [H36<=0] + [H37<=0] + [H38<=0] + [H39<=0] + [H40<=0] + [H41<=0] + [H42<=0] + [H43<=0] + [H44<=0] + [H45<=0] + [H46<=0] + [H47<=0] + [H48<=0] + ...
%                        [H49<=0] + [H50<=0] + [H51<=0] + [H52<=0] + [H53<=0] + [H54<=0] + [H55<=0] + [H56<=0] + [H57<=0] + [H58<=0] + [H59<=0] + [H60<=0] + [H61<=0] + [H62<=0] + [H63<=0] + [H64<=0],[gamma_LQR_obs>=0])

[opt_info] = optimize([[H1<=0] + [H2<=0] + [H3<=0] + [H4<=0] + [H5<=0] + [H6<=0] + [H7<=0] + [H8<=0] + [H9<=0] + [H10<=0] + [H11<=0] + [H12<=0] + [H13<=0] + [H14<=0] + [H15<=0] + [H16<=0] + ...
                       [H17<=0] + [H18<=0] + [H19<=0] + [H20<=0] + [H21<=0] + [H22<=0] + [H23<=0] + [H24<=0] + [H25<=0] + [H26<=0] + [H27<=0] + [H28<=0] + [H29<=0] + [H30<=0] + [H31<=0] + [H32<=0] + ...
                       [H33<=0] + [H34<=0] + [H35<=0] + [H36<=0] + [H37<=0] + [H38<=0] + [H39<=0] + [H40<=0] + [H41<=0] + [H42<=0] + [H43<=0] + [H44<=0] + [H45<=0] + [H46<=0] + [H47<=0] + [H48<=0] + ...
                       [H49<=0] + [H50<=0] + [H51<=0] + [H52<=0] + [H53<=0] + [H54<=0] + [H55<=0] + [H56<=0] + [H57<=0] + [H58<=0] + [H59<=0] + [H60<=0] + [H61<=0] + [H62<=0] + [H63<=0] + [H64<=0] + ...
                       [Y>=0] + [Trace_cond>=0],gamma_LQR_obs>=0], [],ops);

tic
tic
Y_sol  = value(Y);%  Y solution -> WARNING: LQR solves for -Y!!
W1_sol  = value(W1); % W1 solution
W2_sol  = value(W2); % W2 solution
W3_sol  = value(W3); % W3 solution
W4_sol  = value(W4); % W4 solution
W5_sol  = value(W5); % W5 solution
W6_sol  = value(W6); % W6 solution
W7_sol  = value(W7); % W7 solution
W8_sol  = value(W8); % W8 solution
W9_sol  = value(W9); % W9 solution
W10_sol = value(W10); % W10 solution
W11_sol = value(W11); % W11 solution
W12_sol = value(W12); % W12 solution
W13_sol = value(W13); % W13 solution
W14_sol = value(W14); % W14 solution
W15_sol = value(W15); % W15 solution
W16_sol = value(W16); % W16 solution
W17_sol = value(W17); % W17 solution
W18_sol = value(W18); % W18 solution
W19_sol = value(W19); % W19 solution
W20_sol = value(W20); % W20 solution
W21_sol = value(W21); % W21 solution
W22_sol = value(W22); % W22 solution
W23_sol = value(W23); % W23 solution
W24_sol = value(W24); % W24 solution
W25_sol = value(W25); % W25 solution
W26_sol = value(W26); % W26 solution
W27_sol = value(W27); % W27 solution
W28_sol = value(W28); % W28 solution
W29_sol = value(W29); % W29 solution
W30_sol = value(W30); % W30 solution
W31_sol = value(W31); % W31 solution
W32_sol = value(W32); % W32 solution
W33_sol = value(W33); % W33 solution
W34_sol = value(W34); % W34 solution
W35_sol = value(W35); % W35 solution
W36_sol = value(W36); % W36 solution
W37_sol = value(W37); % W37 solution
W38_sol = value(W38); % W38 solution
W39_sol = value(W39); % W39 solution
W40_sol = value(W40); % W40 solution
W41_sol = value(W41); % W41 solution
W42_sol = value(W42); % W42 solution
W43_sol = value(W43); % W43 solution
W44_sol = value(W44); % W44 solution
W45_sol = value(W45); % W45 solution
W46_sol = value(W46); % W46 solution
W47_sol = value(W47); % W47 solution
W48_sol = value(W48); % W48 solution
W49_sol = value(W49); % W49 solution
W50_sol = value(W50); % W50 solution
W51_sol = value(W51); % W51 solution
W52_sol = value(W52); % W52 solution
W53_sol = value(W53); % W53 solution
W54_sol = value(W54); % W54 solution
W55_sol = value(W55); % W55 solution
W56_sol = value(W56); % W56 solution
W57_sol = value(W57); % W57 solution
W58_sol = value(W58); % W58 solution
W59_sol = value(W59); % W59 solution
W60_sol = value(W60); % W60 solution
W61_sol = value(W61); % W61 solution
W62_sol = value(W62); % W62 solution
W63_sol = value(W63); % W63 solution
W64_sol = value(W64); % W64 solution



L1  = inv(Y_sol)*W1_sol.';
L2  = inv(Y_sol)*W2_sol.';
L3  = inv(Y_sol)*W3_sol.';
L4  = inv(Y_sol)*W4_sol.';
L5  = inv(Y_sol)*W5_sol.';
L6  = inv(Y_sol)*W6_sol.';
L7  = inv(Y_sol)*W7_sol.';
L8  = inv(Y_sol)*W8_sol.';
L9  = inv(Y_sol)*W9_sol.';
L10 = inv(Y_sol)*W10_sol.';
L11 = inv(Y_sol)*W11_sol.';
L12 = inv(Y_sol)*W12_sol.';
L13 = inv(Y_sol)*W13_sol.';
L14 = inv(Y_sol)*W14_sol.';
L15 = inv(Y_sol)*W15_sol.';
L16 = inv(Y_sol)*W16_sol.';
L17 = inv(Y_sol)*W17_sol.';
L18 = inv(Y_sol)*W18_sol.';
L19 = inv(Y_sol)*W19_sol.';
L20 = inv(Y_sol)*W20_sol.';
L21 = inv(Y_sol)*W21_sol.';
L22 = inv(Y_sol)*W22_sol.';
L23 = inv(Y_sol)*W23_sol.';
L24 = inv(Y_sol)*W24_sol.';
L25 = inv(Y_sol)*W25_sol.';
L26 = inv(Y_sol)*W26_sol.';
L27 = inv(Y_sol)*W27_sol.';
L28 = inv(Y_sol)*W28_sol.';
L29 = inv(Y_sol)*W29_sol.';
L30 = inv(Y_sol)*W30_sol.';
L31 = inv(Y_sol)*W31_sol.';
L32 = inv(Y_sol)*W32_sol.';
L33 = inv(Y_sol)*W33_sol.';
L34 = inv(Y_sol)*W34_sol.';
L35 = inv(Y_sol)*W35_sol.';
L36 = inv(Y_sol)*W36_sol.';
L37 = inv(Y_sol)*W37_sol.';
L38 = inv(Y_sol)*W38_sol.';
L39 = inv(Y_sol)*W39_sol.';
L40 = inv(Y_sol)*W40_sol.';
L41 = inv(Y_sol)*W41_sol.';
L42 = inv(Y_sol)*W42_sol.';
L43 = inv(Y_sol)*W43_sol.';
L44 = inv(Y_sol)*W44_sol.';
L45 = inv(Y_sol)*W45_sol.';
L46 = inv(Y_sol)*W46_sol.';
L47 = inv(Y_sol)*W47_sol.';
L48 = inv(Y_sol)*W48_sol.';
L49 = inv(Y_sol)*W49_sol.';
L50 = inv(Y_sol)*W50_sol.';
L51 = inv(Y_sol)*W51_sol.';
L52 = inv(Y_sol)*W52_sol.';
L53 = inv(Y_sol)*W53_sol.';
L54 = inv(Y_sol)*W54_sol.';
L55 = inv(Y_sol)*W55_sol.';
L56 = inv(Y_sol)*W56_sol.';
L57 = inv(Y_sol)*W57_sol.';
L58 = inv(Y_sol)*W58_sol.';
L59 = inv(Y_sol)*W59_sol.';
L60 = inv(Y_sol)*W60_sol.';
L61 = inv(Y_sol)*W61_sol.';
L62 = inv(Y_sol)*W62_sol.';
L63 = inv(Y_sol)*W63_sol.';
L64 = inv(Y_sol)*W64_sol.';



polos1_obs  = eig(Ad1 - L1*C);
polos2_obs  = eig(Ad2 - L2*C);
polos3_obs  = eig(Ad3 - L3*C);
polos4_obs  = eig(Ad4 - L4*C);
polos5_obs  = eig(Ad5 - L5*C);
polos6_obs  = eig(Ad6 - L6*C);
polos7_obs  = eig(Ad7 - L7*C);
polos8_obs  = eig(Ad8 - L8*C);
polos9_obs  = eig(Ad9 - L9*C);
polos10_obs = eig(Ad10 - L10*C);
polos11_obs = eig(Ad11 - L11*C);
polos12_obs = eig(Ad12 - L12*C);
polos13_obs = eig(Ad13 - L13*C);
polos14_obs = eig(Ad14 - L14*C);
polos15_obs = eig(Ad15 - L15*C);
polos16_obs = eig(Ad16 - L16*C);
polos17_obs  = eig(Ad17 - L17*C);
polos18_obs  = eig(Ad18 - L18*C);
polos19_obs  = eig(Ad19 - L19*C);
polos20_obs  = eig(Ad20 - L20*C);
polos21_obs  = eig(Ad21 - L21*C);
polos22_obs  = eig(Ad22 - L22*C);
polos23_obs  = eig(Ad23 - L23*C);
polos24_obs  = eig(Ad24 - L24*C);
polos25_obs  = eig(Ad25 - L25*C);
polos26_obs  = eig(Ad26 - L26*C);
polos27_obs  = eig(Ad27 - L27*C);
polos28_obs  = eig(Ad28 - L28*C);
polos29_obs  = eig(Ad29 - L29*C);
polos30_obs  = eig(Ad30 - L30*C);
polos31_obs  = eig(Ad31 - L31*C);
polos32_obs  = eig(Ad32 - L32*C);
polos33_obs  = eig(Ad33 - L33*C);
polos34_obs  = eig(Ad34 - L34*C);
polos35_obs  = eig(Ad35 - L35*C);
polos36_obs  = eig(Ad36 - L36*C);
polos37_obs  = eig(Ad37 - L37*C);
polos38_obs  = eig(Ad38 - L38*C);
polos39_obs  = eig(Ad39 - L39*C);
polos40_obs  = eig(Ad40 - L40*C);
polos41_obs  = eig(Ad41 - L41*C);
polos42_obs  = eig(Ad42 - L42*C);
polos43_obs  = eig(Ad43 - L43*C);
polos44_obs  = eig(Ad44 - L44*C);
polos45_obs  = eig(Ad45 - L45*C);
polos46_obs  = eig(Ad46 - L46*C);
polos47_obs  = eig(Ad47 - L47*C);
polos48_obs  = eig(Ad48 - L48*C);
polos49_obs  = eig(Ad49 - L49*C);
polos50_obs  = eig(Ad50 - L50*C);
polos51_obs  = eig(Ad51 - L51*C);
polos52_obs  = eig(Ad52 - L52*C);
polos53_obs  = eig(Ad53 - L53*C);
polos54_obs  = eig(Ad54 - L54*C);
polos55_obs  = eig(Ad55 - L55*C);
polos56_obs  = eig(Ad56 - L56*C);
polos57_obs  = eig(Ad57 - L57*C);
polos58_obs  = eig(Ad58 - L58*C);
polos59_obs  = eig(Ad59 - L59*C);
polos60_obs  = eig(Ad60 - L60*C);
polos61_obs  = eig(Ad61 - L61*C);
polos62_obs  = eig(Ad62 - L62*C);
polos63_obs  = eig(Ad63 - L63*C);
polos64_obs  = eig(Ad64 - L64*C);

tic

% Vector de valores propios (polos)
% Graficar el círculo unitario
theta = linspace(0, 2*pi, 100);
x_circulo = cos(theta); % Parte real
y_circulo = sin(theta); % Parte imaginaria

figure;
hold on;
plot(x_circulo, y_circulo, 'r--', 'LineWidth', 2); % Círculo unitario
for i = 1:64
    poloeval = eval(['polos' num2str(i) '_obs'])
    plot(real(poloeval), ...
         imag(poloeval), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Valores propios
end
hold off

% Configuración del gráfico
axis equal;
xlabel('Parte real');
ylabel('Parte imaginaria');
title('Valores propios y el círculo unitario');
legend('Círculo unitario', 'Valores propios');
grid on;

tic
n=1
enclave=0;
while(n<=iteraciones)
tic
    Qin = QinFiltrado(n,1);
    Qout = QoutFiltrado(n,1);
    Hin = HinFiltrado(n,1);
    Hout = HoutFiltrado(n,1);

    if (Qin-Qout>=0.8e-03)  && enclave ==0 && n>=20000
        enclave =1;
        tic
    end

    if enclave ==1
        Q1max=max(QinFiltrado);
        Q1min=min(QinFiltrado);

        H2max=17;
        H2min=11;


        Q2max=Q1max-(lam1*sqrt(H2max))+0.001;
        Q2min= Q1min-(lam1*sqrt(H2min))-0.002;

      %  Q2max = 0.025;
      %  Q2min = 0.005;

        H3max=17;
        H3min=5.5;

        Q3max = max(QoutFiltrado);
        Q3min = min(QoutFiltrado);

        z2max = 225;
        z2min = 0;

        Q1Sch = (Q1max - x1) / (Q1max - Q1min);
        H2Sch = (H2max - x2) / (H2max - H2min);
        Q2Sch = (Q2max - x3) / (Q2max - Q2min);
        H3Sch = (H3max - x4) / (H3max - H3min);
        Q3Sch = (Q3max - x5) / (Q3max - Q3min);
        z2Sch = (z2max - x6) / (z2max - z2min);

        Q1 = x1;
        H2 = x2;
        Q2 = x3;
        H3 = x4;
        Q3 = x5;
        z2 = x6;
        lam = x7;
        Hk = [1 0 0 0 0 0 0; 0 0 0 0 1 0 0];

        A = Ar;
        promediolongitud = L;
        dz1 = z1;
        dz2 = (z2-z1);
        dz3 = (L-z2);
        %% HAGO MIS Q ESTIMADAS (Qi gorro j+1)
        Q11=Q1+dt*((((-g*A)*(H2-Hin))/(dz1))-miu1*Q1^2);
        Q21=Q2+dt*((((-g*A)*(H3-H2))/(dz2))-miu2*Q2^2);
        Q31=Q3+dt*((((-g*A)*(Hout-H3))/(dz3))-miu3*Q3^2);

        %% HAGO MIS H ESTIMADAS (Hi gorro j+1)
        H21=H2+dt*((-b^2)*(Q2-Q1+lam1*sqrt((abs(H2))))/(g*A*dz1));
        H31=H3+dt*((-b^2)*(Q3-Q2+lam*sqrt((abs(H3))))/(g*A*dz2));

        phi1 = -miu1*Q11;
        phi2 = (g*Ar*(H21-1))/z1;
        phi3 = (-g*Ar*H21^2)/(z1*z2);
        phi4 = (b^2/(g*Ar*z1))*(1-((lam1*sqrt(abs(H21)))/(2*Q11)));
        phi5 = (-b^2/(g*Ar*z1))*(1+((lam1*sqrt(abs(H21)))/(2*Q21)));
        phi6 = (g*Ar)/(z2-z1);
        phi7 = -miu2*Q21;
        phi8 = (-g*Ar)/(z2-z1);
        phi9 = (b^2)/(g*Ar*(z2-z1));
        phi10 = (-b^2)/(g*Ar*(z2-z1));
        phi11 = (-b^2*sqrt(abs(H31)))/(g*Ar*(z2-z1));
        phi12 = (g*Ar)/(L-z2);
        phi13 = -miu3*Q31;

        Aact1 = [phi1 ,phi2, 0, 0 ,0 ,phi3, 0;
        phi4 ,0, phi5, 0, 0, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        0, 0, phi9 ,0 ,phi10, 0, phi11;
        0, 0, 0, phi12, phi13, 0, 0;
        zeros(2,7)];

        xact1 = [Q11, H21, Q21, H31, Q31, z2, lam]' * (dt / 2);
        Bact1=[(g*Ar)/z1, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z2);
        0, 0;
        0, 0];
        u = [Hin * (dt / 2); Hout * (dt / 2)];
        Ax1 = Aact1 * xact1;
        Bu1 = Bact1 * u;

        phi1 = -miu1*Q1;
        phi2 = (g*Ar*(H2-1))/z1;
        phi3 = (-g*Ar*H2^2)/(z1*z2);
        phi4 = (b^2/(g*Ar*z1))*(1-((lam1*sqrt(abs(H2)))/(2*Q1)));
        phi5 = (-b^2/(g*Ar*z1))*(1+((lam1*sqrt(abs(H2)))/(2*Q2)));
        phi6 = (g*Ar)/(z2-z1);
        phi7 = -miu2*Q2;
        phi8 = (-g*Ar)/(z2-z1);
        phi9 = (b^2)/(g*Ar*(z2-z1));
        phi10 = (-b^2)/(g*Ar*(z2-z1));
        phi11 = (-b^2*sqrt(abs(H3)))/(g*Ar*(z2-z1));
        phi12 = (g*Ar)/(L-z2);
        phi13 = -miu3*Q3;

        Aact2 = [phi1 ,phi2, 0, 0 ,0 ,phi3, 0;
        phi4 ,0, phi5, 0, 0, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        0, 0, phi9 ,0 ,phi10, 0, phi11;
        0, 0, 0, phi12, phi13, 0, 0;
        zeros(2,7)];

        xact2 = [Q1, H2, Q2, H3, Q3, z2, lam]' * (dt / 2);
        Bact2=[(g*Ar)/z1, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z2);
        0, 0;
        0, 0];
        u = [Hin * (dt / 2); Hout * (dt / 2)];
        Ax2 = Aact2 * xact2;
        Bu2 = Bact2 * u;

        xx = [Q1; H2; Q2; H3; Q3; z2; lam];
        dXmat = (Ax1 + Ax2) + (Bu1 + Bu2) + xx;
        dQ1 = dXmat(1);
        dH2 = dXmat(2);
        dQ2 = dXmat(3);
        dH3 = dXmat(4);
        dQ3 = dXmat(5);
        dXm = [dQ1; dH2; dQ2; dH3; dQ3; z2; lam];

        %% Cálculo de los coeficientes de pertenencia `mu`
        Q1SchMax= (1 - Q1Sch);
        Q1SchMin= Q1Sch;

        H2SchMax= (1 - H2Sch);
        H2SchMin= H2Sch;

        Q2SchMax= (1 - Q2Sch);
        Q2SchMin= Q2Sch;

        H3SchMax= (1 - H3Sch);
        H3SchMin= H3Sch;

        Q3SchMax= (1 - Q3Sch);
        Q3SchMin= Q3Sch;

        z2SchMax= (1-z2Sch);
        z2SchMin = z2Sch;

        mu_vals = obtenermusF2(Q1SchMax, Q1SchMin, H2SchMax, H2SchMin, Q2SchMax, Q2SchMin,H3SchMax,H3SchMin,Q3SchMax,Q3SchMin, z2SchMax, z2SchMin, Combinatoria, 6);
        tic

        mu(n) = sum(mu_vals);
        tic
        L_interp=0;
        L_interp = L1 * mu_vals(1) + ...
           L2 * mu_vals(2) + ...
           L3 * mu_vals(3) + ...
           L4 * mu_vals(4) + ...
           L5 * mu_vals(5) + ...
           L6 * mu_vals(6) + ...
           L7 * mu_vals(7) + ...
           L8 * mu_vals(8) + ...
           L9 * mu_vals(9) + ...
           L10 * mu_vals(10) + ...
           L11 * mu_vals(11) + ...
           L12 * mu_vals(12) + ...
           L13 * mu_vals(13) + ...
           L14 * mu_vals(14) + ...
           L15 * mu_vals(15) + ...
           L16 * mu_vals(16) + ...
           L17 * mu_vals(17) + ...
           L18 * mu_vals(18) + ...
           L19 * mu_vals(19) + ...
           L20 * mu_vals(20) + ...
           L21 * mu_vals(21) + ...
           L22 * mu_vals(22) + ...
           L23 * mu_vals(23) + ...
           L24 * mu_vals(24) + ...
           L25 * mu_vals(25) + ...
           L26 * mu_vals(26) + ...
           L27 * mu_vals(27) + ...
           L28 * mu_vals(28) + ...
           L29 * mu_vals(29) + ...
           L30 * mu_vals(30) + ...
           L31 * mu_vals(31) + ...
           L32 * mu_vals(32) + ...
           L33 * mu_vals(33) + ...
           L34 * mu_vals(34) + ...
           L35 * mu_vals(35) + ...
           L36 * mu_vals(36) + ...
           L37 * mu_vals(37) + ...
           L38 * mu_vals(38) + ...
           L39 * mu_vals(39) + ...
           L40 * mu_vals(40) + ...
           L41 * mu_vals(41) + ...
           L42 * mu_vals(42) + ...
           L43 * mu_vals(43) + ...
           L44 * mu_vals(44) + ...
           L45 * mu_vals(45) + ...
           L46 * mu_vals(46) + ...
           L47 * mu_vals(47) + ...
           L48 * mu_vals(48) + ...
           L49 * mu_vals(49) + ...
           L50 * mu_vals(50) + ...
           L51 * mu_vals(51) + ...
           L52 * mu_vals(52) + ...
           L53 * mu_vals(53) + ...
           L54 * mu_vals(54) + ...
           L55 * mu_vals(55) + ...
           L56 * mu_vals(56) + ...
           L57 * mu_vals(57) + ...
           L58 * mu_vals(58) + ...
           L59 * mu_vals(59) + ...
           L60 * mu_vals(60) + ...
           L61 * mu_vals(61) + ...
           L62 * mu_vals(62) + ...
           L63 * mu_vals(63) + ...
           L64 * mu_vals(64);

        Kk = L_interp;
        tic
        dX = dXm + Kk * ([Qin; Qout] - [dQ1; dQ3]);
        x1 = dX(1);
        x2 = dX(2);
        x3 = dX(3);
        x4 = dX(4);
        x5 = dX(5);
        x6 = dX(6);
        x7 = dX(7);

        tic
    else
        x1 = x1i;
        x2 = x2i;
        x3 = x3i;
        x4 = x4i;
        x5 = x5i;
        x6 = x6i;
        x7 = x7i;

        % Inicialización de vectores
        dX = [x1; x2; x3; x4; x5; x6; x7];
        mus = 0;
        mu_vals = zeros(64, 1);  % Inicializamos todos los `mu_i` a 0
    end


for l = 1:64
    mu_Leak2{l}(n) = mu_vals(l); % Almacena el valor
end

    tic

    Q1Fuga2Observador(n)=dX(1);
    H2Fuga2Observador(n)=dX(2);
    Q2Fuga2Observador(n)=dX(3);
    H3Fuga2Observador(n)=dX(4);
    Q3Fuga2Observador(n)=dX(5);
    z2Fuga2Observador(n)=dX(6);
    lamda2Fuga2Observador(n)=dX(7);

    n=n+1
end



tic
t=(out.tout(1000:end));
H2Real = out.H2Real.signals.values(1000:end)
H3Real = out.H3Real.signals.values(1000:end)
z1Real = out.z1Real.signals.values(1000:end)
z2Real = out.z2Real.signals.values(1000:end)
Lam1Real =  out.lam1Real.signals.values(1000:end)
Lam2Real =  out.lam2Real.signals.values(1000:end)


FugaH2 = H2Fuga1Observador_extendido((1:29016))
FugaH2Extend = H2Fuga2Observador(29017:end)

FugaH2Total = [FugaH2, FugaH2Extend];
 tic

output_folder = 'C:\Users\adria\Desktop\DOCUMENTOS TESIS\Resultados\';

%% GRÁFICA DE MUS, LEAK 1 Y LEAK2
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','SUMAMUS','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[min(t) max(t)]); 
ylim(axes1,[-0.1 1.1]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t, MUSL1Fuga1Observador_extendido, 'color', [0.5 0.5 0.5], 'linewidth', 1), hold on  % Un gris intermedio
plot(t, mu, 'color', [0.2 0.2 0.2], 'linewidth', 1), hold on  % Un gris más oscuro
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$\sum \psi_{Leak_1}(t)$';
leyenda2 = '$\sum \psi_{Leak_2}(t)$';
leyenda = legend(leyenda1,leyenda2,'Location','southeast');
print(gcf,fullfile(output_folder, 'MUS.eps'), '-depsc');
tic



%% GRAFICO LA EVOLUCION DE MU'S
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Hin&Hout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[280 max(t)]); 
ylim(axes1,[-0.01 0.35]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Genera colores RGB
colors = lines(64);  % Genera 64 colores distintos
grayColors = repmat(mean(colors, 2), 1, 3);  % Convierte cada color a escala de grises
% Plotea cada mu_Leak2 en un color diferente usando un ciclo for
for l = 1:64
    % Plotea cada mu_Leak2 con su color correspondiente
    plot(t, mu_Leak2{l}, 'Color', grayColors(l, :)); hold on;
    tic
end
print(gcf, fullfile(output_folder, 'musEvolution.eps'), '-depsc');

tic


%% GRÁFICA DE HIN Y HOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Hin&Hout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 max(t)]); 
ylim(axes1,[2 20]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,HinFiltrado, 'color', [0.5 0.5 0.5],'linewidth',0.5), hold on
plot(t,HoutFiltrado,'color', [0.2 0.2 0.2],'linewidth',0.5), hold on
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('[m]','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$H_{in}(t)$';
leyenda2 = '$H_{out}(t)$';
leyenda = legend(leyenda1,leyenda2,'Location','southeast');
print(gcf, fullfile(output_folder, 'presionesEntradaYSalidaLPV.eps'), '-depsc');

%% GRÁFICA DE QIN Y QOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Qin&Qout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 max(t)]); 
ylim(axes1,[0.0158 0.0185]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,QinFiltrado, 'color', [0.5 0.5 0.5],'linewidth',0.5), hold on
plot(t,QoutFiltrado,'color', [0.2 0.2 0.2],'linewidth',0.5), hold on
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('[m^3/s]','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$Q_{in}(t)$';
leyenda2 = '$Q_{out}(t)$';
leyenda = legend(leyenda1,leyenda2,'Location','southeast');
print(gcf,fullfile(output_folder, 'FlujoEntradaYSalidaLPV.eps'), '-depsc');
tic

%% GRÁFICA DE H2REAL,H2OBSERVADOR
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','DIFERENCIASH2','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 max(t)]); 
ylim(axes1,[10 18]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,H2Real,'color', [0.5 0.5 0.5],'linewidth',0.5), hold on
plot(t,FugaH2Total,'color', [0.2 0.2 0.2],'linewidth',0.5), hold on
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('[m]','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$H_{2}(t)$';
leyenda2 =  '$\hat{H}_{2}(t)$';
leyenda = legend(leyenda1,leyenda2,'Location','southeast');
print(gcf,fullfile(output_folder, 'H2EVOLUTION.eps'), '-depsc');
tic

%% GRÁFICA DE H3REAL,H3OBSERVADOR
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','DIFERENCIASH3','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 max(t)]); 
ylim(axes1,[4 16]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,H3Real,'color', [0.5 0.5 0.5],'linewidth',0.5), hold on
plot(t,H3Fuga2Observador,'color', [0.2 0.2 0.2],'linewidth',0.5), hold on
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('[m]','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$H_{3}(t)$';
leyenda2 =  '$\hat{H}_{3}(t)$';
leyenda = legend(leyenda1,leyenda2,'Location','southeast');
print(gcf,fullfile(output_folder, 'H3EVOLUTION.eps'), '-depsc');
tic


%% GRÁFICA DE z1,z2
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 max(t)]); 
ylim(axes1,[20 180]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
plot(t, z1Real, 'color',[0.5 0.5 0.5], 'linewidth', 1.5), hold on
plot(t, z1Fuga1Observador_extendido,'color', [0.5 0.5 0.5], 'linewidth', 0.5), hold on
plot(t, z2Real, '--', 'color',[0.2 0.2 0.2],'linewidth', 1.5), hold on
plot(t, z2Fuga2Observador, '--', 'color', [0.2 0.2 0.2], 'linewidth', 0.5), hold on
% plot(t,z2Real,'k--','linewidth',0.5), hold on
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('[m]','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$z_{1}(t)$';
leyenda2 =  '$\hat{z}_{1}(t)$';
leyenda3 =  '$z_{2}(t)$';
leyenda4 =  '$\hat{z}_{2}(t)$';
legend(leyenda1,leyenda2,leyenda3,leyenda4,'Location','southeast');
print(gcf,fullfile(output_folder, 'leakspositions.eps'), '-depsc');
    tic


%% GRÁFICA DE LAM1,LAM2
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','magnitudleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
factor = 20;
t_reducido = downsample(t, factor);
LAM1Fuga1Observador_reducido = downsample(LAM1Fuga1Observador_extendido, factor);
lamda2Fuga2Observador_reducido = downsample(lamda2Fuga2Observador, factor);

% plot(t_reducido, LAM1Fuga1Observador_reducido, 'color', [0.5 0.5 0.5], 'linewidth', 1.5), hold on
% plot(t_reducido, lamda2Fuga2Observador_reducido, '--', 'color', [0.2 0.2 0.2], 'linewidth', 1.5), hold on
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 max(t)]); 
ylim(axes1,[0 3.5e-04]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
plot(t, Lam1Real, 'color',[0.5 0.5 0.5], 'linewidth',1.5), hold on
%plot(t, LAM1Fuga1Observador_extendido, 'color', [0.5 0.5 0.5], 'linewidth', 0.01), hold on
plot(t_reducido, LAM1Fuga1Observador_reducido, 'color', [0.5 0.5 0.5], 'linewidth', 0.5), hold on
plot(t, Lam2Real, '--','color',[0.2 0.2 0.2], 'linewidth', 1.5), hold on
plot(t_reducido, lamda2Fuga2Observador_reducido, '--', 'color', [0.2 0.2 0.2], 'linewidth', 0.5), hold on
% plot(t, lamda2Fuga2Observador,'--', 'color', [0.2 0.2 0.2], 'linewidth', 0.01), hold on
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('[m^{5/2}/s]','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$\lambda_{1}(t)$';
leyenda2 =  '$\hat{\lambda}_{1}(t)$';
leyenda3 =  '$\lambda_{2}(t)$';
leyenda4 =  '$\hat{\lambda}_{2}(t)$';
leyenda = legend(leyenda1,leyenda2,leyenda3,leyenda4,'Location','southeast');
print(gcf,fullfile(output_folder, 'lambdasmagnitudes.eps'), '-depsc');
    tic


%% GRÁFICA DE QinFiltrado-QoutFiltrado
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','magnitudleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[min(t) max(t)]); 
ylim(axes1,[0 0.0020]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
plot(t, QinFiltrado-QoutFiltrado,'color', [0.5 0.5 0.5], 'linewidth', 0.5), hold on
grid on
xlabel('[s]','FontName','Times New Roman','FontSize', 20);
ylabel('[m^{5/2}/s]','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$|Q_{in}-Q_{out}|$';
leyenda = legend(leyenda1,'Location','southeast');
print(gcf,fullfile(output_folder, 'QRESTA.eps'), '-depsc');
tic

% 
% figure(1)
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
% plot(out.QinReal.signals.values((1000:end)),'b','linewidth',1);
% hold on
% plot(Q1Fuga2Observador,'r','linewidth',1)
% title('Qin (Q1) Flow rate ','FontName','Times New Roman','FontSize', 20);
% xlabel('Time [s]','FontName','Times New Roman','FontSize', 20);
% ylabel('Floe rate [m^3/s]','FontName','Times New Roman','FontSize', 20);
% legend('Qin-Real', 'Qin-Observador'); % Añade una leyenda para identificar las líneas
% grid on
% ylim([0.015 0.02]);  
% hold off
% saveas(gcf, fullfile(output_folder, 'Qin_Flow_rate.png'))  % Guardar en PNG
% 
% 
% figure(2)
% plot(out.H2Real.signals.values((1000:end)),'b','linewidth',1);
% hold on
% plot(H2Fuga2Observador,'r','linewidth',1)
% title('H2 Pressure head','FontName','Times New Roman','FontSize', 20);
% xlabel('Time [s]','FontName','Times New Roman','FontSize', 20);
% ylabel('Pressure head [m]','FontName','Times New Roman','FontSize', 20);
% legend('H2-Real', 'H2-Observador'); % Añade una leyenda para identificar las líneas
% grid on
% ylim([6 18]);  
% hold off
% saveas(gcf, fullfile(output_folder, 'H2.png'))  % Guardar en PNG
% 
% 
% figure(3)
% plot(out.Q2Real.signals.values((1000:end)),'b','linewidth',1);
% hold on
% plot(Q2Fuga2Observador,'r','linewidth',1)
% title('Q2 Flow rate','FontName','Times New Roman','FontSize', 20);
% xlabel('Time [s]','FontName','Times New Roman','FontSize', 20);
% ylabel('Floe rate [m^3/s]','FontName','Times New Roman','FontSize', 20);
% legend('Q2-Real', 'Q2-Observador'); % Añade una leyenda para identificar las líneas
% grid on
% ylim([0.015 0.02]);  
% hold off
% saveas(gcf, fullfile(output_folder, 'Q2_Flow_rate.png'))  % Guardar en PNG
% 
% 
% figure(4)
% plot(out.H3Real.signals.values((1000:end)),'b','linewidth',1);
% hold on
% plot(H3Fuga2Observador,'r','linewidth',1)
% title('H3 Pressure head','FontName','Times New Roman','FontSize', 20);
% xlabel('Time [s]','FontName','Times New Roman','FontSize', 20);
% ylabel('Pressure head [m]','FontName','Times New Roman','FontSize', 20);
% legend('H3-Real', 'H3-Observador'); % Añade una leyenda para identificar las líneas
% grid on
% ylim([4 18]);  
% hold off
% saveas(gcf, fullfile(output_folder, 'H3.png'))  % Guardar en PNG
% 
% 
% figure(5)
% plot(out.QoutReal.signals.values((1000:end)),'b','linewidth',1);
% hold on
% plot(Q3Fuga2Observador,'r','linewidth',1)
% title('Qout (Q3) Flow rate','FontName','Times New Roman','FontSize', 20);
% xlabel('Time [s]','FontName','Times New Roman','FontSize', 20);
% ylabel('Flow rate [m^3/s]','FontName','Times New Roman','FontSize', 20);
% legend('Qout-Real', 'Qout-Observador'); % Añade una leyenda para identificar las líneas
% grid on
% ylim([0.015 0.02]);  
% hold off
% saveas(gcf, fullfile(output_folder, 'Qout_Flow_rate.png'))  % Guardar en PNG
% 
% figure(6)
% plot(out.z2Real.signals.values((1000:end)),'b--','linewidth',1);
% hold on
% plot(z2Fuga2Observador,'r','linewidth',1)
% title('Leak 2 Position','FontName','Times New Roman','FontSize', 20);
% xlabel('Time [s]','FontName','Times New Roman','FontSize', 20);
% ylabel('Position [m]','FontName','Times New Roman','FontSize', 20);
% legend('z2-Real', 'z2-Observador'); % Añade una leyenda para identificar las líneas
% ylim([0 200]);  
% grid on
% hold off
% saveas(gcf, fullfile(output_folder, 'z2.png'))  % Guardar en PNG
% 
% 
% figure(7)
% plot(out.lam2Real.signals.values((1000:end)),'b','linewidth',1);
% hold on
% plot(lamda2Fuga2Observador,'r','linewidth',1)
% title('Magnitude of the Leak 2','FontName','Times New Roman','FontSize', 20);
% xlabel('Time [s]','FontName','Times New Roman','FontSize', 20);
% ylabel('Leak coefficient [m^{5/2}/s]','FontName','Times New Roman','FontSize', 20);
% legend('lam2-Real', 'lam2-Observador'); % Añade una leyenda para identificar las líneas
% grid on
% hold off
% saveas(gcf, fullfile(output_folder, 'LAM2.png'))  % Guardar en PNG














%%FUNCIÓN PARA CREAR COMBINATORIA (MAX,MIN)
function Combinatoria = matrizVerdad(NoParametrosAdaptar)
numerobits=NoParametrosAdaptar;
tamano=2^numerobits;
matriz=zeros(tamano,numerobits);
n=0;
acumulador=0;
inicializador=0;

for columna=numerobits:-1:1
inicializador=(2^(n))+1;
incrementos = 2^(n+1);
    if inicializador<2 
    inicializador = 2;
    end
    for filas=incrementos:incrementos:tamano
        matriz(inicializador:filas,columna)=1;
        inicializador=inicializador+incrementos;
    end
n=n+1;
end

    Combinatoria = matriz;  % Asignar la matriz al argumento de salida
end


%%% FUNCIONES PARA FUGA2
function [adis,bdis,cdis,ddis,matricesA,matricesB] = evaluarF2(z1,lam1,x1max,x1min,x2max,x2min,x3max,x3min,x4max,x4min,x5max,x5min,x6max,x6min,matrizcombinatoria,NoParametros,Cd,Dd,dt)


matricesA = cell(2^NoParametros,1);
matricesB = cell(2^NoParametros,1);
adis = cell(2^NoParametros, 1);
bdis = cell(2^NoParametros, 1);
cdis = cell(2^NoParametros, 1);
ddis = cell(2^NoParametros, 1);
SYS = cell(2^NoParametros, 1);
SYSd = cell(2^NoParametros, 1);

D=0.1048;
Ar=3.1416*(D/2)^2;
L=225;
b=721.6531;
f1=0.0276;
f2=0.0276;
f3=0.0276;
miu1=f1/(2*D*Ar);
miu2=f2/(2*D*Ar);
miu3=f3/(2*D*Ar);
g=9.81;

for i = 1:2^NoParametros

    
    % Leer el valor de cada bit en la fila i (todas las columnas de la fila i)
    fila = matrizcombinatoria(i, :);
    
    % Asignar x1 según el primer valor de la fila
    if fila(1) == 0
        x1 = x1min;
    else
        x1 = x1max;
    end
    
    % Asignar x2 según el segundo valor de la fila
    if fila(2) == 0
        x2 = x2min;
    else
        x2 = x2max;
    end
    
    % Asignar x3 según el tercer valor de la fila
    if fila(3) == 0
        x3 = x3min;
    else
        x3 = x3max;
    end
    
    % Asignar x4 según el cuarto valor de la fila
    if fila(4) == 0
        x4 = x4min;
    else
        x4 = x4max;
    end

    % Asignar x5 según el cuarto valor de la fila
    if fila(5) == 0
        x5 = x5min;
    else
        x5 = x5max;
    end

    % Asignar x6 según el cuarto valor de la fila
    if fila(6) == 0
        x6 = x6min;
    else
        x6 = x6max;
    end
    
phi1 = -miu1*x1;
phi2 = (g*Ar*(x2-1))/z1;
phi3 = (-g*Ar*x2^2)/(z1*x6);
phi4 = (b^2/(g*Ar*z1))*(1-((lam1*sqrt(abs(x2)))/(2*x1)));
phi5 = (-b^2/(g*Ar*z1))*(1+((lam1*sqrt(abs(x2)))/(2*x3)));
phi6 = (g*Ar)/(x6-z1);
phi7 = -miu2*x3;
phi8 = (-g*Ar)/(x6-z1);
phi9 = (b^2)/(g*Ar*(x6-z1));
phi10 = (-b^2)/(g*Ar*(x6-z1));
phi11 = (-b^2*sqrt(abs(x4)))/(g*Ar*(x6-z1));
phi12 = (g*Ar)/(L-x6);
phi13 = -miu3*x5;

matricesA{i} = [phi1 ,phi2, 0, 0 ,0 ,phi3, 0;
        phi4 ,0, phi5, 0, 0, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        0, 0, phi9 ,0 ,phi10, 0, phi11;
        0, 0, 0, phi12, phi13, 0, 0;
    zeros(2,7)];


        matricesB{i}=[(g*Ar)/z1, 0;
    0, 0;
    0, 0;
    0, 0;
    0, -(g*Ar)/(L-x6);
    0, 0;
    0, 0];
    
    dt = dt;
    % Crear el sistema continuo en espacio de estados
    SYS{i} = ss(matricesA{i}, matricesB{i}, Cd, Dd);

    % Discretizar el sistema
    SYSd{i} = c2d(SYS{i,1}, dt);

    % Extraer las matrices discretas
    [adis{i}, bdis{i}, cdis{i}, ddis{i}] = ssdata(SYSd{i,1});
    tic
end
    
end

function mu_vals = obtenermusF2(Q1SchMax, Q1SchMin, H2SchMax, H2SchMin, Q2SchMax, Q2SchMin,H3SchMax,H3SchMin,Q3SchMax,Q3SchMin, z2SchMax, z2SchMin, MatrizCombinatoria, NoParametros)

    % Inicializamos el vector mu_vals con ceros
    mu_vals = zeros(2^NoParametros, 1);  % Se crea un vector de ceros con el tamaño adecuado

    for i = 1:(2^NoParametros)
        % Obtener la fila actual de la Matriz Combinatoria
        fila = MatrizCombinatoria(i, :);

        % Asignar los valores dependiendo del bit
        if fila(1) == 1
            Q1Sch = Q1SchMax;
        else
            Q1Sch = Q1SchMin;
        end

        if fila(2) == 1
            H2Sch = H2SchMax;
        else
            H2Sch = H2SchMin;
        end

        if fila(3) == 1
            Q2Sch = Q2SchMax;
        else
            Q2Sch = Q2SchMin;
        end

        if fila(4) == 1
            H3Sch = H3SchMax;
        else
            H3Sch = H3SchMin;
        end

        if fila(5) == 1
            Q3Sch = Q3SchMax;
        else
            Q3Sch = Q3SchMin;
        end

        if fila(6) == 1
            z2Sch = z2SchMax;
        else
            z2Sch = z2SchMin;
        end

        % Calcular y almacenar el valor de mu_vals(i)
        mu_vals(i) = Q1Sch * H2Sch * Q2Sch * H3Sch * Q3Sch * z2Sch;
    end

end

