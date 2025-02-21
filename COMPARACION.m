%Limpieza de pantalla
clear all
close all
clc

%Declaración de variables simbólicas
syms th1(t) th2(t) t l1 l2 q1(t) q2(t) q3(t) l3 

%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP=[0 0];

%Creamos el vector de coordenadas articulares
Q= [th1, th2];
disp('Coordenadas generalizadas');
pretty (Q);

%Creamos el vector de velocidades generalizadas
Qp= diff(Q, t);
disp('Velocidades generalizadas');
pretty (Qp);
%Número de grado de libertad del robot
GDL= size(RP,2);
GDL_str= num2str(GDL);

%Junta 1 
%Posición de la junta 1 respecto a 0
P(:,:,1)= [l1*cos(th1); l1*sin(th1);0];
%Matriz de rotación de la junta 1 respecto a 0
R(:,:,1)= [cos(th1) -sin(th1)  0;
           sin(th1)  cos(th1)  0;
           0         0         1];


%Junta 2
%Posición de la junta 1 respecto a 0
P(:,:,2)= [l2*cos(th2); l2*sin(th2);0];
%Matriz de rotación de la junta 1 respecto a 0
R(:,:,2)= [cos(th2) -sin(th2)  0;
           sin(th2)  cos(th2)  0;
           0         0         1];

%Creamos un vector de ceros
Vector_Zeros= zeros(1, 3);

%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las matrices de transformación Homogénea globales
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las posiciones vistas desde el marco de referencia inercial
PO(:,:,GDL)= P(:,:,GDL); 
%Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO(:,:,GDL)= R(:,:,GDL); 
%Inicializamos las INVERSAS de las matrices de rotación vistas desde el marco de referencia inercial
RO_inv(:,:,GDL)= R(:,:,GDL); 

for i = 1:GDL
    i_str= num2str(i);
    disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
    pretty (A(:,:,i));

   %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end
    disp(strcat('Matriz de Transformación global T', i_str));
    T(:,:,i)= simplify(T(:,:,i));
    pretty(T(:,:,i))

    RO(:,:,i)= T(1:3,1:3,i);
    RO_inv(:,:,i)= transpose(RO(:,:,i));
    PO(:,:,i)= T(1:3,4,i);
    %pretty(RO(:,:,i));
    %pretty(RO_inv(:,:,i));
    %pretty(PO(:,:,i));
end

%Calculamos la matriz de transformación del marco de referencia inercial
%visto desde el actuador final
% disp(strcat('Matriz de Transformación T', GDL_str,'_O calculada de forma manual'));
% RF_O=RO_inv(:,:,GDL);
% PF_O=-RF_O*PO(:,:,GDL);
% TF_O= simplify([RF_O PF_O; Vector_Zeros 1]);
% pretty(TF_O);   

%disp(strcat('Matriz de Transformación T', GDL_str,'_O calculada de forma auntomática'));
%pretty(simplify(inv(T(:,:,GDL))));


%Calculamos el jacobiano lineal de forma diferencial
disp('Jacobiano lineal obtenido de forma diferencial');
%Derivadas parciales de x respecto a th1 y th2
Jv11= functionalDerivative(PO(1,1,GDL), th1);
Jv12= functionalDerivative(PO(1,1,GDL), th2);
%Derivadas parciales de y respecto a th1 y th2
Jv21= functionalDerivative(PO(2,1,GDL), th1);
Jv22= functionalDerivative(PO(2,1,GDL), th2);
%Derivadas parciales de z respecto a th1 y th2
Jv31= functionalDerivative(PO(3,1,GDL), th1);
Jv32= functionalDerivative(PO(3,1,GDL), th2);

%Creamos la matríz del Jacobiano lineal
jv_d=simplify([Jv11 Jv12;
              Jv21 Jv22;
              Jv31 Jv32]);
pretty(jv_d);


%Calculamos el jacobiano lineal de forma analítica
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);

for k= 1:GDL
    if RP(k)==0 %Casos: articulación rotacional
       %Para las juntas de revolución
        try
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
            Jw_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));
            Jw_a(:,k)=[0,0,1];
        end
        
        %Para las juntas prismáticas
     elseif RP(k)==1 %Casos: articulación prismática
%         
        try
            Jv_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)=[0,0,1];
        end
            Jw_a(:,k)=[0,0,0];
     end
 end    

Jv_a= simplify (Jv_a);
Jw_a= simplify (Jw_a);
V=simplify (Jv_a*Qp');
W=simplify (Jw_a*Qp');

% Configuración del robot antropomórfico, 0 para junta rotacional, 1 para junta prismática
joint_types_A=[0 0 0]; % Robot antropomórfico

% Creamos el vector de coordenadas articulares para el robot antropomórfico
Q_A= [q1 q2 q3];
disp('Coordenadas articulares (Robot Antropomórfico)');
pretty (Q_A);

% Creamos el vector de velocidades articulares para el robot antropomórfico
Qp_A= diff(Q_A, t); % Utilizo diff para derivadas cuya variable de referencia no depende de otra: ejemplo el tiempo
disp('Velocidades articulares (Robot Antropomórfico)');
pretty (Qp_A);

% Número de grado de libertad del robot antropomórfico
GDL_A= size(joint_types_A,2); % ***Siempre se coloca 3, ya que indica la dimensión de las columnas
GDL_str_A= num2str(GDL_A); % Convertimos el valor numérico a una cadena de carácteres tipo string

% Articulación 1 para el robot antropomórfico
% Posición de la junta 1 respecto a 0
P_A(:,:,1)= [0;
             0;
             l1]; % *** Vector de posición indexado por página

% Articulación 2 para el robot antropomórfico
P_A(:,:,2)=[l2*cos(q2); 
          l2*sin(q2);
                   0];

% Articulación 3 para el robot antropomórfico
P_A(:,:,3)=[l3*cos(q3); 
          l3*sin(q3);
                   0];

% Matriz de rotación de la articulación 1 respecto a 0
R_A(:,:,1)= [cos(q1)  0  sin(q1);
           sin(q1)  0 -cos(q1);
           0         1        0];

R_A(:,:,2)= [cos(q2) -sin(q2) 0;
           sin(q2)  cos(q2) 0;
           0         0        1];

R_A(:,:,3)= [cos(q3) -sin(q3) 0;
           sin(q3)  cos(q3) 0;
           0         0        1];

% Inicializamos las matrices de transformación Homogénea locales
A_A(:,:,GDL_A)=simplify([R_A(:,:,GDL_A) P_A(:,:,GDL_A); Vector_Zeros 1]);

% Inicializamos las matrices de transformación Homogénea globales
T_A(:,:,GDL_A)=simplify([R_A(:,:,GDL_A) P_A(:,:,GDL_A); Vector_Zeros 1]);

% Inicializamos los vectores de posición vistos desde el marco de referencia inercial
PO_A(:,:,GDL_A)= P_A(:,:,GDL_A); 

% Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO_A(:,:,GDL_A)= R_A(:,:,GDL_A);

for i = 1:GDL_A
    i_str_A= num2str(i);
    % Locales
    A_A(:,:,i)=simplify([R_A(:,:,i) P_A(:,:,i); Vector_Zeros 1]);

    % Globales
    try
       T_A(:,:,i)= T_A(:,:,i-1)*A_A(:,:,i);
    catch
       T_A(:,:,i)= A_A(:,:,i);  % Caso específico cuando i=1 nos marcaría error en try
    end
    disp(strcat('Matriz de Transformación global T', i_str_A));
    T_A(:,:,i)= simplify(T_A(:,:,i));
    pretty(T_A(:,:,i));

    % Obtenemos la matriz de rotación "RO "y el vector de translación PO de la
    % matriz de transformación Homogénea global T_A(:,:,GDL_A)
    RO_A(:,:,i)= T_A(1:3,1:3,i);
    PO_A(:,:,i)= T_A(1:3,4,i);
    pretty(RO_A(:,:,i));
    pretty(PO_A(:,:,i));
end

% Calculamos el jacobiano lineal de forma diferencial
disp('Jacobiano lineal obtenido de forma diferencial (Robot Antropomórfico)');
% Derivadas parciales de x respecto a q1 
Jv_A11= functionalDerivative(PO_A(1,1,GDL_A), q1);
Jv_A12= functionalDerivative(PO_A(1,1,GDL_A), q2);
Jv_A13= functionalDerivative(PO_A(1,1,GDL_A), q3);
% Derivadas parciales de y respecto a q1 
Jv_A21= functionalDerivative(PO_A(2,1,GDL_A), q1);
Jv_A22= functionalDerivative(PO_A(2,1,GDL_A), q2);
Jv_A23= functionalDerivative(PO_A(2,1,GDL_A), q3);
% Derivadas parciales de z respecto a q1, q2, q3
Jv_A31= functionalDerivative(PO_A(3,1,GDL_A), q1);
Jv_A32= functionalDerivative(PO_A(3,1,GDL_A), q2);
Jv_A33= functionalDerivative(PO_A(3,1,GDL_A), q3);
% Creamos la matriz del Jacobiano lineal
jv_d_A=simplify([Jv_A11 Jv_A12 Jv_A13;
               Jv_A21 Jv_A22 Jv_A23;
               Jv_A31 Jv_A32 Jv_A33]); 
pretty(jv_d_A);

% Calculamos el jacobiano lineal de forma analítica
% Inicializamos jacobianos analíticos (lineal y angular)
Jv_a_A(:,GDL_A)=PO_A(:,:,GDL_A);
Jw_a_A(:,GDL_A)=PO_A(:,:,GDL_A);

for k= 1:GDL_A
    if ((joint_types_A(k)==0)|(joint_types_A(k)==1)) % Casos: articulación rotacional y prismática

       % Para las articulaciones rotacionales
        try
            Jv_a_A(:,k)= cross(RO_A(:,3,k-1), PO_A(:,:,GDL_A)-PO_A(:,:,k-1)); % *****
            Jw_a_A(:,k)= RO_A(:,3,k-1);
        catch
            Jv_a_A(:,k)= cross([0,0,1], PO_A(:,:,GDL_A)); % Matriz de rotación de 0 con respecto a 0 es la Matriz Identidad, la posición previa también será 0
            Jw_a_A(:,k)=[0,0,1]; % Si no hay matriz de rotación previa se obtiene la Matriz identidad
         end
     else
        % Para las articulaciones prismáticas
        try
            Jv_a_A(:,k)= RO_A(:,3,k-1);
        catch
            Jv_a_A(:,k)=[0,0,1]; % Si no hay matriz de rotación previa se obtiene la Matriz identidad
        end
            Jw_a_A(:,k)=[0,0,0];
     end
end    

Jv_a_A= simplify (Jv_a_A);
Jw_a_A= simplify (Jw_a_A);
V_A=simplify (Jv_a_A*Qp_A');
W_A=simplify (Jw_a_A*Qp_A');
% Mostrar jacobianos y velocidades para el robot convencional
disp('----------- Robot Planar -----------');
disp('Jacobiano lineal obtenido de forma analítica');
pretty(Jv_a);
disp('Jacobiano angular obtenido de forma analítica');
pretty(Jw_a);
disp('Velocidad lineal obtenida mediante el Jacobiano lineal');
pretty(V);
disp('Velocidad angular obtenida mediante el Jacobiano angular');
pretty(W);

% Separación para mayor claridad
disp('------------------------------------------');

% Mostrar jacobianos y velocidades para el robot antropomórfico
disp('----------- Robot Antropomórfico -----------');
disp('Jacobiano lineal obtenido de forma analítica (Robot Antropomórfico)');
pretty(Jv_a_A);
disp('Jacobiano angular obtenido de forma analítica (Robot Antropomórfico)');
pretty(Jw_a_A);
disp('Velocidad lineal obtenida mediante el Jacobiano lineal (Robot Antropomórfico)');
pretty(V_A);
disp('Velocidad angular obtenida mediante el Jacobiano angular (Robot Antropomórfico)');
pretty(W_A);
