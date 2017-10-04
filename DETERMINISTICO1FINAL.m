load appl1.mat % CARREGANDO MATRIZES 

u = dtrend(u); % REMOVENDO O LINEAR FIT
y = dtrend(y);

m = min(size(u)); 			% Numero de Entradas
l = min(size(y)); 			% Numero de Saidas
ax_id = [1:900];n_id = length(ax_id); 	% Setando tamanho da matriz para identificação(1:900)
ax_val = [901:1247];n_val = length(ax_val); % Setando tamanho da matriz de validação (901:1247)
u_id = u(ax_id,:);y_id = y(ax_id,:); 	% Matrizes para identificação
u_val = u(ax_val,:);y_val = y(ax_val,:); % Matrizes para validação
i = 10; 				% Numero de blocos linha

disp(['     Numero de Entradas:             ',num2str(m)]);
disp(['     Numero de Saidas:            ',num2str(l)]);
disp(['     Amostras para identificação:     ',num2str(n_id)]);
disp(['     Amostras para validação    ',num2str(n_val)]);
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

                               	

R = triu(qr([U;Y]'))';            % CÁLCULO DO FATOR DE DECOMPOSIÇÃO R
clear U Y
R = R(1:2*i*(m+l),1:2*i*(m+l)); 	


mi2  = 2*m*i; % DOBRO DAS LINHAS EXISTENTES


Rf = R((2*m+l)*i+1:2*(m+l)*i,:); 	% MATRIZ DE HENKEL SAIDAS FUTURAS(DECOMPOSTA)
Rp = [R(1:m*i,:);R(2*m*i+1:(2*m+l)*i,:)]; % MATRZ DE HENNKEL, ENTRADAS E SAIDAS PASSADAS(DECOMPOSTA)
Ru  = R(m*i+1:2*m*i,1:mi2); 	% MATRIZ DE HENKEL ENTRADAS FUTURAS(DECOMPOSTA)
Rfp = [Rf(:,1:mi2) - (Rf(:,1:mi2)/Ru)*Ru,Rf(:,mi2+1:2*(m+l)*i)]; %PERPENDICULARES FUTURAS(DECOMPOSTA)
Rpp = [Rp(:,1:mi2) - (Rp(:,1:mi2)/Ru)*Ru,Rp(:,mi2+1:2*(m+l)*i)];%PERPENDICULARES PASSADAS(DECOMPOSTA)




 if (norm(Rpp(:,(2*m+l)*i-2*l:(2*m+l)*i),'fro')) < 1e-10 % EVITAR DEFICIÊNCIA DE POSIÇÕES DE MATRIZ
 Ob  = (Rfp*pinv(Rpp')')*Rp; 	% CALCULO DA PROJEÇÃO OBLIQUA Ob = Yf/Wp(Uf)
 else
 Ob = (Rfp/Rpp)*Rp;
 end
 
 
 [U,S,V] = svd(Ob); % VALORES SINGULARES
 ss = diag(S);
 clear V S WOW
 
 n = 8; % ORDEM DO SISTEMA
 U1 = U(:,1:n); 				% DETERMINANDO U1 A PARTIR DA ORDEM DO SISTEMA
 
gam  = U1*diag(sqrt(ss(1:n))); % GAMA
gamm = gam(1:l*(i-1),:); % GAMA i-1
gam_inv  = pinv(gam); % PSEUDO INVERSA DE GAMA
gamm_inv = pinv(gamm); % PSEUDO INVERSA DE GAMA i-1

% CALCULO DA SEGUNDA PROJEÇÃO OBLIQUA
Rf = R((2*m+l)*i+l+1:2*(m+l)*i,:); 	% SAIDAS FUTURAS
Rp = [R(1:m*(i+1),:);R(2*m*i+1:(2*m+l)*i+l,:)]; % ENTRADAS E SAIDAS PASSADAS
Ru  = R(m*i+m+1:2*m*i,1:mi2); 		% ENTRADAS FUTURAS
% PERPENDICULAR FUTURA 
Rfp = [Rf(:,1:mi2) - (Rf(:,1:mi2)/Ru)*Ru,Rf(:,mi2+1:2*(m+l)*i)]; 
% PERPENDICULAR PASSADA
Rpp = [Rp(:,1:mi2) - (Rp(:,1:mi2)/Ru)*Ru,Rp(:,mi2+1:2*(m+l)*i)]; 
if (norm(Rpp(:,(2*m+l)*i-2*l:(2*m+l)*i),'fro')) < 1e-10
  Obm  = (Rfp*pinv(Rpp')')*Rp; 		% PROJEÇÃO OBLIQUA
else
  Obm = (Rfp/Rpp)*Rp;
end

% DETERMINANDO OS ESTADOS Xi e Xi+1
Xi  = gam_inv  * Ob;
Xip = gamm_inv * Obm;


% RESOLVENDO O SISTEMA DE EQUAÇÕES 
Rhs = [       Xi   ;  R(m*i+1:m*(i+1),:)]; 
Lhs = [      Xip   ;  R((2*m+l)*i+1:(2*m+l)*i+l,:)];
sol = Lhs/Rhs;

% EXTRAIR AS MATRIZES DA SOLUÇÃO
A = sol(1:n,1:n);
B = sol(1:n,n+1:n+m);
C = sol(n+1:n+l,1:n);
D = sol(n+1:n+l,n+1:n+m);


[ys,ers] = simul(y_val,u_val,A,B,C,D);
ers
figure(1)
subplot(211);plot([y_val(:,1),ys(:,1)]);
title('Real Azul e simulado laranja')
subplot(212);plot([y_val(:,2),ys(:,2)]);

pause



