function plotMesh(x,xw)
% plotMesh  dibuja la malla
dy = (x(2)-x(1))/2;
  
  %Graficar 
  figure(4)
  %lineas horizontales
  xh = repmat([xw(1) xw(length(xw))]',1,2); dyh=repmat([-dy dy],2,1); plot(xh,dyh,'-k');
  hold('on');
  %lineas frontera entre volumenes
  xv = [xw xw]'; dyv = repmat([-dy dy]',1,length(xw)); plot(xv,dyv,'--k');
  % puntos (nodos)
  plot(x(1),0,'*r');plot(x(length(x)),0,'*r'); plot(x(2:(length(x)-1)),zeros(1,length(x)-2),'*b');
  hold('off')
  axis('equal'); 