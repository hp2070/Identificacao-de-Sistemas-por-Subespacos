load appl1.mat % CARREGANDO MATRIZES 

u = dtrend(u); % REMOVENDO O LINEAR FIT
y = dtrend(y);

m = min(size(u)); 			% Numero de Entradas
l = min(size(y)); 			% Numero de Saidas
ax_id = [1:900];n_id = length(ax_id); 	% Setando tamanho da matriz para identifica��o(1:900)
ax_val = [901:1247];n_val = length(ax_val); % Setando tamanho da matriz de valida��o (901:1247)
u_id = u(ax_id,:);y_id = y(ax_id,:); 	% Matrizes para identifica��o
u_val = u(ax_val,:);y_val = y(ax_val,:); % Matrizes para valida��o
i = 10; 				% Numero de blocos linha

disp(['     Numero de Entradas:             ',num2str(m)]);
disp(['     Numero de Saidas:            ',num2str(l)]);
disp(['     Amostras para identifica��o:     ',num2str(n_id)]);
disp(['     Amostras para valida��o    ',num2str(n_val)]);
disp(['     Numero de blocos Linha:         ',num2str(i)]);


figure(1);hold off;subplot;clf;
subplot(221);plot(u_id(:,1));title('Input 1');
subplot(222);plot(u_id(:,2));title('Input 2');
subplot(223);plot(y_id(:,1));title('Output 1');
subplot(224);plot(y_id(:,2));title('Output 2');

pause


[l,ny] = size(y_id);if (ny < l);y_id = y_id';[l,ny] = size(y_id);end
[m,nu] = size(u_id);if (nu < m);u_id = u_id';[m,nu] = size(u_id);end


j = ny-2*i+1;
AUXin = [];
Uaux = u(1,1);


Y=zeros(l*2*i,j);                               % MATRIZ DE HENKEL DAS SAIDAS
y_idnormal = y_id/sqrt(j);
for k=1:2*i
	Y((k-1)*l+1:k*l,:)=y_idnormal(:,k:k+j-1);
end

U=zeros(l*2*i,j);                               % MATRIZ DE HENKEL DAS SAIDAS
u_idnormal = u_id/sqrt(j);
for k=1:2*i
	U((k-1)*l+1:k*l,:)=u_idnormal(:,k:k+j-1);
end

                               	
disp('      Computando Fator R');
R = triu(qr([U;Y]'))';            % C�LCULO DO FATOR DE DECOMPOSI��O R
clear U Y
R = R(1:2*i*(m+l),1:2*i*(m+l)); 	


mi2  = 2*m*i; % DOBRO DAS LINHAS EXISTENTES


Rf = R((2*m+l)*i+1:2*(m+l)*i,:); 	% MATRIZ DE HENKEL SAIDAS FUTURAS(DECOMPOSTA)
Rp = [R(1:m*i,:);R(2*m*i+1:(2*m+l)*i,:)]; % MATRZ DE HENNKEL, ENTRADAS E SAIDAS PASSADAS(DECOMPOSTA)
Ru  = R(m*i+1:2*m*i,1:mi2); 	% MATRIZ DE HENKEL ENTRADAS FUTURAS(DECOMPOSTA)
Rfp = [Rf(:,1:mi2) - (Rf(:,1:mi2)/Ru)*Ru,Rf(:,mi2+1:2*(m+l)*i)]; %PERPENDICULARES FUTURAS(DECOMPOSTA)
Rpp = [Rp(:,1:mi2) - (Rp(:,1:mi2)/Ru)*Ru,Rp(:,mi2+1:2*(m+l)*i)];%PERPENDICULARES PASSADAS(DECOMPOSTA)




 if (norm(Rpp(:,(2*m+l)*i-2*l:(2*m+l)*i),'fro')) < 1e-10 % EVITAR DEFICI�NCIA DE POSI��ES DE MATRIZ
 Ob  = (Rfp*pinv(Rpp')')*Rp; 	% CALCULO DA PROJE��O OBLIQUA Ob = Yf/Wp(Uf)
 else
 Ob = (Rfp/Rpp)*Rp;
 end
 
 
 [U,S,V] = svd(Ob);
 ss = diag(S);
 clear V S WOW
 
 
 n = 8;
  
%ENTRAR AQUI A ORDEM DO SISTEMA
U1 = U(:,1:n); 				%DETERMINAMOS U1 A PARTIR DA ORDEM ESCOLHIDA


gam  = U1*diag(sqrt(ss(1:n))); % GAMA COM TRA�O EM CIMA
gamm = U1(1:l*(i-1),:)*diag(sqrt(ss(1:n)));%DETERMINANDO GAMA i
gam_per  = U(:,n+1:l*i)'; 		% COMPLEMENTO ORTOGONAL GAMA i
gamm_inv = pinv(gamm); 			% PSEUDO INVERSA GAMA i

A = gamm_inv*gam(l+1:l*i,:); % DETERMINAMOS A DIRETAMENTE A E C
C = gam(1:l,:);% DETERMINAR C 
M = gam_per*(R((2*m+l)*i+1:2*(m+l)*i,:)/R(m*i+1:2*m*i,:)); % DETERMINANDO MATRIZ M 
L = gam_per; % DETERMINANDO MATRIZ L

% RESOLVENDO SISTEMA DE EQUA��O
Lhs = zeros(i*(l*i-n),m);
Rhs = zeros(i*(l*i-n),l*i);
aa = 0;
for k=1:i
  Lhs((k-1)*(l*i-n)+1:k*(l*i-n),:) = M(:,(k-1)*m+1:k*m);
  Rhs((k-1)*(l*i-n)+1:k*(l*i-n),1:(i-k+1)*l) = L(:,(k-1)*l+1:l*i);
end
Rhs = Rhs*[eye(l),zeros(l,n);zeros(l*(i-1),l),gamm];

sol = Rhs\Lhs;
B = sol(l+1:l+n,:);
D = sol(1:l,:);



[ys,ers] = simul(y_val,u_val,A,B,C,D);
ers
figure(1)
subplot(211);plot([y_val(:,1),ys(:,1)]);
title('Real Azul e simulado laranja')
subplot(212);plot([y_val(:,2),ys(:,2)]);

pause
