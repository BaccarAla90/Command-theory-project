%Les données
Ks=5.468;jl=0.003882;jeq=0.00035;Rm=15.5;Km=0.0089;Kt=0.0089;Beq=0.001;Kg=1/36;
Bl=-0.01;
Const=(((Km^2)*(Kg^2))/(Rm*jeq))-(Bl+Beq)/jeq;
%Les matrices du modèles d'état
A=[0 1 0 0;-Ks/jl Bl/jl Ks/jl -Bl/jl; 0 0 0 1; -Ks/jeq Bl/jeq Ks/jeq Const ];
B=[0; 0; (Kg*Km)/(Rm*jeq);-(Kg*Km)/(Rm*jeq)];
C=[1 1 0 0];
[D,poles]=eig(A);
%Etude de commandabilité et d'observabilité du système
Co = ctrb(A,B);%Matrice de commandabilité
dCo=det(Co);
oB=obsv(A,C); %Matrice d'observabilité
doB=det(oB);
%Controle par retour d'état
%u=-Kx
desired_poles=[-17.5+17.8i,-17.5-17.8i,-35+35.7i,-35-35.7i];
K=place(A,B,desired_poles);
%[D,poles]=eig(A-B*K);
%u=-Kx+Nv
%v??
X=inv(A-B*K);
invN=(-C*X)*B;
N=1/invN;
D=0;
[b,a]=ss2tf(A,B,C,D); %Détermination de la fonction de transfert à partir du répresentation
  %d'état
sysc=tf(b,a); %Fonction de transfert continu
sysd = c2d(sysc,0.1); %Fonction de transfert discrete
%Commande RST
phi=[1         0             0             0          0          0        0       0;
    -4.85e05   1             0             0          1041      0         0        0 ;
    1.011e06   -4.85e05      1             0          -179.7    1041       0      0;
    -5.1129e05 1.011e06     -4.85e05      1           -2049     -179.7     1041    0;
    10.11      -5.1129e05   1.011e06     -4.85e05      1230      -2049    -179.7   1041;
    0         10.11        -5.1129e05    1.011e06       0         1230    -2049    -179.7;
    0         0             10.11       -5.1129e05     0         0        1230    -2049; 
    0         0             0             10.11         0         0         0       1230];
P=[1; -1.3; 0.48; -1.3; 0; 0; 0; 0];
RST=linsolve(phi,P);

s=filt(1,[1 -0.1725 -1.9682 1.1814 ]);
