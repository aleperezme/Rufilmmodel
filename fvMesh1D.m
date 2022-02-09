function [x,xw,dx] = fvMesh1D(nx,xlen)
% fvMesh1D  Crea malla 1D para volumenes finito. 
%
% [x,xw,dx] = fvMesh1D(nx,xlen,r)
%
% Input: nx   = Numero de nodos. Default: nx = 10
%        xlen = longitud de la malla.  Default:  xlen = 1
% Output:  x = vector con la localizacion de los nodos, incluyendo las fronteras (0, xlen) 
%          xw = vector con las fronteras de cada volumen de control 
%               xw(1) = 0,  xw(nx+2) = xlen
%          dx = ancho de los volumenes de control

if nargin<1,  nx = 10;   end
if nargin<2,  xlen = 1;  end

% --- Crear malla
dx = xlen/(nx);          
% --- Ubicar nodos
x = zeros(nx,1); x(1) = dx/2;
for i=2:nx
  x(i) = x(i-1) + dx;
end
% --- Ubicar limites de los volumenes
xw = zeros(size(x));
for i=2:nx
  xw(i) = (x(i-1)+x(i))/2;
end
xw(1)=0; xw(nx+1) = xlen;
%plotMesh(x,xw) %grafica la malla

