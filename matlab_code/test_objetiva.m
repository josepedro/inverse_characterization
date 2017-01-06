X = [24930 0.91 1 4.1000e-05 1.4400e-04 64];

[Zs1 alpha_simples] = Allard_rigido(0.025,2000,X(1),X(2),X(3),X(4),X(5));
[Zs2 alpha_duplo] = Allard_rigido(2*0.025,2000,X(1),X(2),X(3),X(4),X(5));

%load('A_exp')
load('Z_exp')
Z_referencia_esp1 = Z_3.';
F_obj_1 = sum(abs((Z_referencia_esp1(500:1650)/413)-(Zs1(500:1650)/413)).^2);

disp(F_obj_1);

Z_referencia_esp2 = Z_6.';
F_obj_2 = sum(abs((Z_referencia_esp2(500:1650)/413)-(Zs2(500:1650)/413)).^2);

disp(F_obj_2);

F_obj = 3*F_obj_1 + F_obj_2
