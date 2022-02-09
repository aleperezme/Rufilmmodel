clear all; close all
area=0.07;

%read data from file
[f,fp] = uigetfile('MultiSelect', 'on'); fp = strcat(fp,f); %readnamefiles
 
for k = 1:length(fp) 
    data = dlmread(cell2mat(fp(k)),'',1, 1); %read datafrinfile k
    [m,n] = size(data); data = data(1:m,1:n-1); %chooseusefuldata
    cur(:,k) = data(:,2); %save current
end
pot = data(3:m/2-1,1); vel=str2double(f)/1000; %potential and scan rate (potential equal in all files)
cre= cur(3:m/2-1,:)/area  ; cox = cur(end-1:-1:m/2+3,:)/area; %current ox and re

%regresion lineal (kdfv^0.5+kcpv.5=i)
%{[v^-0.5 1]*[k2;k1]=[i/v]}={Xa*ka=ya} or {[1 v^0.5]*[k1;k2]=[i/v^0.5]}={Xb*kb=yb}
Xb=[ones(k,1), (vel.^0.5)']; %matrix X
%1. average approx. (cox-cre) 
ya=((cox-cre)./vel.^0.5)'; km=Xb\ya; % y=cur(pot,vel)

%errors
df=(k-size(Xb,2)); sse_a =sum((ya-Xb*km).^2); sst_a=sum((ya-mean(ya)).^2); 
rmsea=(sse_a/df).^.5; Rsqa = 1- sse_a./sst_a; se_ka = (diag(inv(Xb'*Xb))*rmsea.^2).^0.5; 
tska = km./se_ka; pval_ka = 2*(1-tcdf(abs(tska),df)); gofm=[rmsea; Rsqa; se_ka; pval_ka];

%2. biased cre and cox
ya=((cox-cur(m/2+1,:)/area)./vel.^0.5)'; yb=((cre-cur(1,:)/area)./vel.^0.5)';%cox-coxi
%ka = (inv(Xa'*Xa))*Xa'*ya; kb = (inv(Xb'*Xb))*Xb'*yb; %parametros k1 k2 for each pot
kox=Xb\ya; kre=Xb\yb; %parametros kre kox for each pot 

%errors
df=(k-size(Xb,2)); sse_a=sum((ya-Xb*kox).^2); sst_a=sum((ya-mean(ya)).^2); 
rmsea=(sse_a/df).^.5; Rsqa = 1- sse_a./sst_a; se_ka = (diag(inv(Xb'*Xb))*rmsea.^2).^0.5; 
tska = kox./se_ka; pval_ka = 2*(1-tcdf(abs(tska),df)); gofa=[rmsea; Rsqa; se_ka; pval_ka];

df=(k-size(Xb,2)); sse_b=sum((yb-Xb*kre).^2); sst_b=sum((yb-mean(yb)).^2); 
rmseb=(sse_b/df).^.5; Rsqb = 1- sse_b./sst_b; se_kb = (diag(inv(Xb'*Xb))*rmseb.^2).^0.5; 
tskb = kre./se_kb; pval_kb = 2*(1-tcdf(abs(tskb),df)); gofb=[rmseb; Rsqb; se_kb; pval_kb];

% data correct (neg/pos kdf when is close to zero)
ind = find(km(1,:)<0); km(1,ind)= (km(1,ind)-min(km(1,:)))*0.1;
ind = find(kre(1,:)>0); kre(1,ind)= (kre(1,ind)-max(kre(1,:)))*0.1;
ind = find(kox(1,:)<0); kox(1,ind)= (kox(1,ind)-min(kox(1,:)))*0.1;

%extracting % 
kdfm= mean([abs(kre(1,:));kox(1,:)]); kcpm = mean([abs(kre(2,:));kox(2,:)]);
idf = km(1,:)'*vel.^0.5 ; icp=km(2,:)'*vel; percp=icp./(icp+idf);
idfm = kdfm'*vel.^0.5 ; icpm=kcpm'*vel; percpm=icpm./(icpm+idfm);
idfo = kox(1,:)'*vel.^0.5 ; icpo=kox(2,:)'*vel; percpo=icpo./(icpo+idfo);
idfr = kre(1,:)'*vel.^0.5 ; icpr=kre(2,:)'*vel; percpr=abs(icpr)./(abs(icpr)+abs(idfr));

%selecting 
selV=[pot(1:7:15); pot(21:5:31); pot(36:4:48); pot(53:5:end)]; 
selA=[cre(1:7:15,:); cre(21:5:31,:); cre(36:4:48,:); cre(53:5:end,:);... 
   cox(1:7:15,:);cox(21:5:31,:); cox(36:4:48,:); cox(53:5:end,:)]'; 
selP= [[selV; selV]'; selA];  

%saving 
writetable(table([[0;vel'] selP]),'sp.dat') 
fid=fopen('kivsv5.dat','w'); fprintf(fid,'%.6g %.6g\n',km); 
fid=fopen('gofivsv5.dat','w'); fprintf(fid,'%.6g %.6g %.6g %.6g %.6g %.6g\n',gofm);
fid=fopen('korc.dat','w'); fprintf(fid,'%%.6g %.6g %.6g %.6g\n',[kox; kre]); 
fid=fopen('gofoxcor.dat','w'); fprintf(fid,'%.6g %.6g %.6g %.6g %.6g %.6g\n',gofa);
fid=fopen('gofrecor.dat','w'); fprintf(fid,'%%.6g %.6g %.6g %.6g %.6g %.6g\n',gofb);
fid=fopen('ks.dat','w'); fprintf(fid,'%.6g %.6g\n',[kdfm;kcpm]);
fid=fopen('percp.dat','w'); fprintf(fid,'%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g \n',percp');
fid=fopen('percpm.dat','w'); fprintf(fid,'%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g \n',percpm');

%plot 
figure(1); plot(pot,percpo); figure(2); plot(pot,percpr); figure(3); plot(pot,percp); legend(num2str(vel'*1000)); % 
figure(4); plot(pot,percpm);

%regresion lineal (kdfv^0.5+kcpv.5=i) = {[1v^0.5]*[kdf;kcp]=[i/v^0.5]}={Xb*kb=yb} % 
%y data capacitive current for ploting 
%ya=((cox-cur(m/2+1,:)/area)./vel.^0.5)'; yb=((cre-cur(1,:)/area)./vel.^0.5)'; %yo=cur at pot x at diff
ya=((cox)./vel.^0.5)'; yb=((cre)./vel.^0.5)'; %yo=cur at pot x at diff
%in= round(length(pot)/20,0); ya = ya(1:in:end,:); yb = yb(1:inc:end,:); %select less data (20 default) 
% kox =(inv(Xb'*Xb))*Xb'*ya; kre = (inv(Xb'*Xb))*Xb'*yb; %parametros kdf kcp %

%calculo parametros
kox= Xb\ya; ioxd=kox(1,:)'*vel.^0.5; ioxc=kox(2,:)'*vel; 
kre=Xb\yb; ired=kre(1,:)'*vel.^0.5; irec=kre(2,:)'*vel; % red and ox currents  

%errors
sse_a=sum((ya-Xb*kox).^2); sst_a=sum((ya-mean(ya)).^2); rmsea=(sse_a/df).^.5; Rsqa = 1- sse_a./sst_a;
sse_b=sum((yb-Xb*kre).^2); sst_b=sum((yb-mean(yb)).^2); rmseb=(sse_b/df).^.5; Rsqb = 1- sse_b./sst_b;
se_ka =(diag(inv(Xb'*Xb))*rmsea.^2).^0.5 ; tska = kox./se_ka; pval_ka = 2*(1-tcdf(abs(tska),df));
se_kb =(diag(inv(Xb'*Xb))*rmseb.^2).^0.5; tskb = kre./se_kb; pval_kb = 2*(1-tcdf(abs(tskb),df));
gofa = [rmsea; Rsqa; se_ka; pval_ka]; gofb=[rmseb; Rsqb; se_kb; pval_kb];

%saving
fid=fopen('koxred.dat','w'); fprintf(fid,'%.6g %.6g %.6g %.6g\n',[kox; kre]); 
fid=fopen('gofox.dat','w');fprintf(fid,'%.6g %.6g %.6g %.6g %.6g %.6g\n',gofa); 
fid=fopen('gofre.dat','w');fprintf(fid,'%.6g %.6g %.6g %.6g %.6g %.6g\n',gofb); 

% plot
figure(5);plot([pot;pot(end:-1:1)],[irec; ioxc(end:-1:1,:)]); legend(num2str(vel'*1000)); 
figure(6);plot([pot;pot(end:-1:1);pot;pot(end:-1:1)],[cre(:,1:5:end);cox(end:-1:1,1:5:end);...
     irec(:,1:5:end);ioxc(end:-1:1,1:5:end)]); legend(num2str(vel(1:5:end)'*1000));