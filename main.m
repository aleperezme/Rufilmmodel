function [] = main() 
close all
% --- universal constants
F = 96485.3329; R = 8.314472;  %C/mol, J/(molK)
% --- system constants
%solution reaction
n=1; T= 298; p.alp = 0.5; p.E0=0.21; p.ks=8; % #e-, K, [],V, 1/s

%film characteristics g0 = total load in 1 monolayer (nmol/cm2)
Nav = 6.022*10^23; amol=71.36; p.gs = 10^9/(amol*Nav*10^-16); nl=32;
%molec/mol,        a2/molec        nmol/mol       a=cm2/a2
p.Eado = (-6900-1200*(nl-1))/nl; p.Eadr = (-5000-6000*(nl-1))/nl; p.g = -1.25; 
%J/mol,                         J/mol,                         Ei=RTg/g0

%electron difussion
p.D=1e-4; p.L = 26e-8*nl; %cm/s, cm; 

% CV experimental setup scan rate v, E1 and E2 initial and final pot
p.v=0.05; p.E1=0; p.E2=0.62; tcv = 2*abs(p.E1-p.E2)/p.v; %V s-1, V, V 

% --- parameters calculated no depends on experiment
p.ec=n*F/(R*T);  p.te = 1./(R*T); % C/J, mol/J

% --- saving parameters
fid=fopen('paramsim.dat','w'); fprintf(fid,'%.6g \n',[p.Eado p.Eadr p.g p.D p.L tcv]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- simulation
opt = odeset('refine',4,'RelTol',1e-8); %options ode
%create mesh
nz = nl; zl = 1; [z,zw,dz]=fvMesh1D(nz,zl);
% 1. equilibration (nc=1)
nc =1; p.tf = 15; p.kn = p.D*p.tf/p.L; %experiment time and dif - parameter
IC = zeros(length(z),1);IC=IC+3/5;
[t, c] = ode15s(@(t,c)echemfilm(t,c,z,dz,nc,p),[0 1],IC,opt); plot(t,c);
%cyclic voltammetry (nc=2)
nc=2; df=p.E1-p.E2; p.tf=2*abs(df)/p.v; p.kn = p.D*p.tf/p.L;
IC=c(end,:)'; [t, c] = ode15s(@(t,c)echemfilm(t,c,z,dz,nc,p),[0 1],IC,opt);
V=genV(t,nc,p); figure(2); plot(V,c); hold on; plot(V,1-c); 
i=cur(t,c,nc,p); figure(3);plot(V,i); hold on;
%cyclic voltammetry (nc=2)
nc=2; df=p.E1-p.E2; p.tf=2*abs(df)/p.v; p.kn = p.D*p.tf/p.L;
IC=c(end,:)'; [t, c] = ode15s(@(t,c)echemfilm(t,c,z,dz,nc,p),[0 1],IC,opt);
%plot and saving files
V=genV(t,nc,p); figure(2); plot(V,c); hold on; plot(V,1-c); 
i=cur(t,c,nc,p); figure(3);plot(V,i);
fid=fopen('i-cv.dat','w'); fprintf(fid,'%.6g %.6g\n',[V';i']);
for i=1:length(c(1,:))
    fid=fopen(strcat('g-cv',num2str(i),'.dat'),'w'); 
    fprintf(fid,'%.6g %.6g \n',[V';c(:,i)']); 
end

%verifying kinetic model leads to same resut than Nernst
 E=linspace(0,0.62,26)'; gam = zeros(length(E),2);
 for i=1:length(E)
     Ea =E(i); %Ea is a choosen E
     if i== 1 
         x0 = [1/3*p.gs 2/3*p.gs]; %initial guess for the gr and go
     else
         x0 = gam(i-1,:); %guess for other value should be close to the last 
     end
     %fmincon solve (1) go + gr = g0 (2) -go + gr*exp(ec*eta+adin/RT) = 0
     %with go=x(1),gr=x(2), so Eq. (1) is [1 1]*x = g0 --> A*x = b                                
     % adin=f(x(1),x(2),E), so Eq. (2) is nonlinear and is constructed in 
     % the function equil2. The values of x are bounded between 0 and g0.
     A=[1 1]; b=1; low = [0 0]; up = [1 1];  
     gam(i,:) = fmincon(@(x)equil2(x,Ea,p),x0,[],[],A,b,low,up);
end
%plot and saving
 figure(2);plot(E,gam,'g'); hold on;
 fid=fopen('gamma-eq.dat','w'); fprintf(fid,'%.6g %.6g %.6g\n',[E';gam']); 
return

function f = echemfilm (t,c,z,dz,nc,p)
% Generate odes for transient diffusion 1D Eq wit reaction in z=0. 
% with spatial variables discretized with central differences 
%Equation
% d(gr/Dx)/dt= D*d2gr/dx2  
% at x = 0, dGR/dt = Rgen-JRout  (Rgen=kr*go-ko*gr) (JRout=-D*d(gr)/dx) 
% with CR=gr/Dx R=gr/gs z= x/L t'=t/tf dimensionless.
% dR/dt'= tf*Dz*D/L d2R/dz2   
% Rigth hand discrete:
% dR/dt = (tf*D/(Dz*L))*(R_i-1 + -2 R_i + R_i+1) 
% dR/dt = [aE aP aW 0 ...; 0 aE aP aW...]*C(i) i=1,...nx
% In:    t  = time variable (as tspan, or t vector).
%        z =  nodes position vector, including frontiers (x(1),x(m)).
%        dz  = cell width 
%        p = system parameter structure (tf,D,L,krxn,eta)
% Out: f = dR/dt' 

% --- variables R = gr/gs 
m = length(z); R = c; 

% --- Vector ctes
b = zeros (m,1);

% --- Calcula coeficientes de DifCen
aW = zeros(m,1); aP=aW;
aW = aW + p.kn/dz; aP = aP - 2*p.kn/dz; aE= aW;

% --- Boundary Conditions
pot = genV(t,nc,p);eta = pot-p.E0;
[ko, kr] = calck(eta,1-R(1),R(1),p);
aP(1)=-(p.kn/dz+(kr+ko)*p.tf); b(1) = kr*p.tf;   %  Boundary in 0

if m == 1
    B = [aW aP aE];
    A = spdiags(B,-1:1,m,2); 
    f = sum(A*R)+b;
else
    aP(m) = -p.kn/dz; b(m) = 0;  %  Boundary in L.
    B = [aW aP aE]; A = spdiags(B,-1:1,m,m);
    f = A*R+b;
end
return

function pot = genV(t,n,p)
switch n
    case 1
    pot = p.E1; %pot signal
    case 2
    df= p.E1-p.E2; pot= p.E2+df*abs(abs(df)-p.v*p.tf*t)/abs(df); % pot signal
    otherwise
         warning('Unexpected element. No circuit created.')
end
return

function [ko, kr] = calck(eta,o,r,p)
% calculate reaction rates depends on Ead(J/s), eta(V), go,gr(dimensionless)
% and parameters te=1/RT, ec=nF/RT and g interaction (aR-R+aO-O-2aR-O)  
% reaction rate parameters
ko=p.ks*exp(p.te*p.Eadr)*exp((1-p.alp)*p.ec*eta).*exp(-p.g*r); %1/s 
kr=p.ks*exp(p.te*p.Eado)*exp(-(p.alp)*p.ec*eta).*exp(-p.g*o);  %1/s
return

function f = cur(t,c,nc,p)
pot = genV(t,nc,p); eta= pot-p.E0; %overpotential V
%balance
% superficial concentration mol/cm2
gr=c(:,1); go=1-gr; %X=g (mol/cm2) 
% reaction rate parameters
[ko, kr] = calck(eta,go,gr,p);
% reaction rates
rox = ko.*gr*p.gs;  %mol/(s cm2)
rre = kr.*go*p.gs;  %mol/(s cm2) 
%balance total
f(:,1) = (rox-rre)*p.ec/p.te; 
return

function f = equil2(x,e1,p)
%constants
go=x(1); gr=x(2); %variables
% parameters calculated
eta= e1-p.E0; % C/J, overpotential V
adin= p.Eadr-p.Eado-(p.g*(gr-go)/p.te); 
%    go + gr - p.g0 = 0
%    -go + gr*exp(ec*eta) = 0
f(1) = (-go + gr*exp(eta*p.ec+adin*p.te))^2;
return