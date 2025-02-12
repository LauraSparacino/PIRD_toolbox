% PROGRAMMA PER FILTRO AUTOREGRESSIVO, PASSA ALTO E PASSA BASSO

% sig è un file multivariato con can canali
% can è il canale che si vuole filtrare
% p è l'ordine del filtro AR

function [fia,fib]=AR_filter(sig,can,p)

s=sig(:,can);

s1=size(s,1);%dimensione dato

usc = zeros(s1,1);
fib = zeros(s1,1);

% "andata" 
usc(1)= s(1);
for i = 2 : s1
   usc(i)= p*usc(i-1)+(1-p)*s(i);
end 

% "ritorno"
fib(s1)=usc(s1);
for i = (s1-1):-1:1
     fib(i)= p*fib(i+1)+(1-p)*usc(i);
end
   
%passa alto
fia = s-fib+mean(s);


% % crea un vettore contenente segnale filtrato e da filtrare 
% % e li visualizza su un unico grafico
% figure(2);
% subplot(2,1,1);
% plot ([s fib]),title('originale+passa basso');
% zoom xon;
% subplot(2,1,2);
% plot ([s fia]),title('originale+passa alto');
% zoom xon;