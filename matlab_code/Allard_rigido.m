
function [Zs alpha] = Allard_rigido(L,fmax,sigma,phi,alpha_inf,Lambda_material,Lambda_l_material)

%% Modelo de material poroso de estrutura rigido
% C�lcula densidade efetiva din�mica e compressibilidade efetiva para
% modelos de flu�do equivalente de Johnson-Lafarge. 
% [Zs alpha] = Lafarge(L,freq,sigma,phi,alpha_inf,Lambda_material,Lambda_l_material)
% Par�mentros de entrada:
% L = comprimento da amostra;
% freq = frequ�ncia m�xima;
% sigma = resistividade ao fluxo [Ns/m^2];
% phi = porosidade [%]
% alpha_inf = tortuosidade (normalmente varia de 1 � 4)
% Lambda_material = Comprimento caracter�stico viscoso
% Lambda_l_material = Comprimento caracter�stico t�rmico

% Como sa�da tem-se:
% Zs = impedancia de superf�cie
% Alpha = Coeficiente de absor��o
% 
L_material = L;     % espessura do material poroso 
rho_0 = 1.204;      % Densidade [kg/m^3]
T = 20;               % Temperatura
c_0 = (331.2+0.6*T);% velocidade de propaga��o no meio [m/s]
eta = 1.84e-5;      % Coeficiente de viscosidade, viscosidade absoluta ou viscosidade din�mica
gamma = 1.4;        % Raz�o de calores espec�ficos (gamma=c_p/c_v)
k_f = 0.026;        % Condutibilidade t�rmica do fluido [W/mK]
P_0 = 101320;       % Press�o est�tica do meio [Pa]
c_p = 1.0035e3;     % Coeficiente de Calor espec�fico a press�o constante (c_p), (c_v � coef. calor esp. a VOLUME cte)
Pr = eta*c_p/k_f;   % N�mero de Prandtl (???)--- c_p calor espec�fico, k_f � 

% Par�metros Macro
sigma_material = sigma;      % Resistividade ao fluxo [kN/s]
phi_material =  phi;         % Porosidade 60%
alpha_inf_material = alpha_inf;    % Tortuosidade
Lambda = Lambda_material;
Lambda_l = Lambda_l_material;
q_0 = eta/sigma;
% q_l_0 = q_0*Lambda_l/Lambda;
q_l_0 = phi_material*Lambda_l^2/8;

% M = rho_1+phi_material*rho_0;

f = 0:1:fmax;

w = 2.*pi.*f;%% MODELO DE Johnson - Lafarge

D = (eta*(phi_material^2)*(Lambda^2));
C = (4*w*rho_0*(q_0^2)*(alpha_inf_material^2));
B = ((phi_material*eta)./(1i*w*rho_0*alpha_inf_material*q_0));
A = 1 + B .* ( 1 + 1i*C ./ D).^(1/2);
rho_rigido = rho_0*alpha_inf_material * A; % densidade din�mica efetiva

disp('rho_rigido');
mean(rho_rigido)

% rho_limp = (rho_ef*M - rho_0^2)./(M+rho_ef-2*rho_0);
F = eta*phi_material^2*Lambda_l^2;
E = 4*rho_0*Pr*q_l_0^2.*w;
D =  1 + (1i*E ./ F);
C = (phi_material*eta)./ (1i*w*rho_0*Pr*q_l_0);
B = 1 + C .* D.^(1/2);
A = gamma-(gamma-1) ./ B;
K_ef = gamma*P_0 ./ A; % m�dulo de compressibilidade
disp('compressibilidade');
mean(K_ef)

Zc = sqrt(rho_rigido.*K_ef);
kc_L = w.*((rho_rigido./K_ef).^(1/2));
Zs = -1i.*Zc./phi_material .*cot(kc_L.*L_material);
alpha = 1 - (abs((Zs - rho_0*c_0)./(Zs + rho_0*c_0))).^2;