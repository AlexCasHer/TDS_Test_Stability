%% PROGRAMA PARA SISTEMAS CON RETARDOS: SISTEMAS DISCTRIBUIDOS VEHÍCULOS CONECTADOS 
%% LIMPIEZA DEL ÁREA DE TRABAJO:
clc
clear
hold on
format long
tic

%% PARÁMETROS GENERALES DEL PROGRAMA:
r = 5;

val_exp=2;

%% Vectores de almacenamiento:
k1_vector_Th1 = zeros();
k2_vector_Th1 = zeros();
k1_vector_Th2 = zeros();
k2_vector_Th2 = zeros();
indice_Th1 = 1;
indice_Th2 = 1;

%% Retardos del sistema:
h1 = 0.05;          %Originalmente el retardo era de 0.1 y lo dejé en 0.1
h2 = 3*h1;


for PUNTOS=1:1   %Definición de intervalos de graficación:
   
        %Definición del número de puntos.
        puntos  = 120;
        epsilon = 0.00005;

        %Paso para a:
        k1_min = 0;
        k1_max =  20;
        paso_k1=(k1_max-k1_min)/puntos;

        %Paso para b:
        k2_min = 0;
        k2_max =  20; %Inicialmente el valor es 30
        paso_k2=(k2_max-k2_min)/puntos;
        
end


%% PROGRAMA PRINCIPAL.
for k1 = k1_min:paso_k1:k1_max
    for k2= k2_min:paso_k2:k2_max
        
        if (h1==0 && h2==0)
        else

            for SISTEMA=1:1 % DEFINICIÓN DEL SISTEMA

                % Definición de los retardos del sistema
                h0 = 0;

                %Máximo Común Divisor:
                %h_base=gcd(100000*h1,100000*h2)/100000;
                h_base = h1;

                %Definición de las matrices del sistema con retardos:
                A1 = [  0  -2.5;
                    -2.5   0];

                A0 = zeros(length(A1));

                A2 = zeros(length(A1));

                C0 = zeros(length(A1));

                C1 = [-k1/(2*h_base)          0;
                    0          -k2/(2*h_base)];

                C2 = C1;

                G = C1;
            end

            for PARAMETROS_INICIALES=1:1 %PARÁMETROS INICIALES REQUERIDOS
                % Dimensiones de la matrices del sistema
                n = length(A0);

                % Matriz Identidad de dimensión n x n y de n^2 x n^2:
                In = eye(n);
                Inn = eye(n^2);

                % Matriz Cero de dimensión n x n;
                Ceros_n = zeros(n);


                %Cálculo del Número de Dinámicas:
                H=max(h1,h2);
                m=round(H/h_base);
                N_Dinamicas=round(4*m-1);
            end

            for MATRIZ_A=1:1 % MATRIZ A.
                %Generación del bloque 1 de la Matriz A:
                A_Main_B1=zeros(n^2,(m+1)*n^2);

                A_delta0=round(h0/h_base);
                A_delta1=round(h1/h_base);
                A_delta2=round(h2/h_base);

                A_Main_B1(:,A_delta0*n^2+1:(A_delta0+1)*n^2)=A_Main_B1(:,A_delta0*n^2+1:(A_delta0+1)*n^2) + kron(A0',In);
                A_Main_B1(:,A_delta1*n^2+1:(A_delta1+1)*n^2)=A_Main_B1(:,A_delta1*n^2+1:(A_delta1+1)*n^2) + kron(A1',In);
                A_Main_B1(:,A_delta2*n^2+1:(A_delta2+1)*n^2)=A_Main_B1(:,A_delta2*n^2+1:(A_delta2+1)*n^2) + kron(A2',In);

                %Generación del bloque 2 de la Matriz A:
                A_Main_B2=zeros(n^2,(m+1)*n^2);

                A_delta0=round((H-h0)/h_base);
                A_delta1=round((H-h1)/h_base);
                A_delta2=round((H-h2)/h_base);

                A_Main_B2(:,A_delta0*n^2+1:(A_delta0+1)*n^2)=A_Main_B2(:,A_delta0*n^2+1:(A_delta0+1)*n^2) - kron(In',A0');
                A_Main_B2(:,A_delta1*n^2+1:(A_delta1+1)*n^2)=A_Main_B2(:,A_delta1*n^2+1:(A_delta1+1)*n^2) - kron(In',A1');
                A_Main_B2(:,A_delta2*n^2+1:(A_delta2+1)*n^2)=A_Main_B2(:,A_delta2*n^2+1:(A_delta2+1)*n^2) - kron(In',A2');

                %Matriz de ceros para A:
                A=zeros(2*m*n^2);

                %Primer Bloque de la Matriz A:
                for i=1:m
                    A((i-1)*n^2+1:i*n^2,(i-1)*n^2+1:(m+i)*n^2)=A_Main_B1;
                end
                %Segundo Bloque de la Matriz A:
                for i=m:2*m-1
                    A(i*n^2+1:(i+1)*n^2,(i-m)*n^2+1:(i+1)*n^2)=A_Main_B2;
                end
                    
            end

            for MATRIZ_C=1:1 % MATRIZ C.
                % Generación del bloque 1 de la Matriz C:
                C_Main_B1=zeros(n^2,m*n^2);

                C_delta0=round(h0/h_base);
                C_delta2=round((h2-2*h_base)/h_base);
                C_delta3=round((h2-h_base)/h_base);

                C_Main_B1(:,C_delta0*n^2+1:(C_delta0+1)*n^2)=C_Main_B1(:,C_delta0*n^2+1:(C_delta0+1)*n^2) + kron(C0',In);
                C_Main_B1(:,C_delta2*n^2+1:(C_delta2+1)*n^2)=C_Main_B1(:,C_delta2*n^2+1:(C_delta2+1)*n^2) + kron(C1',In);
                C_Main_B1(:,C_delta3*n^2+1:(C_delta3+1)*n^2)=C_Main_B1(:,C_delta3*n^2+1:(C_delta3+1)*n^2) + kron(C2',In);

                % Generación del bloque 2 de la Matriz C:
                C_Main_B2=zeros(n^2,m*n^2);

                C_delta0=round((H-h0)/h_base);
                C_delta2=round((H-(h2-2*h_base))/h_base);
                C_delta3=round((H-(h2-h_base))/h_base);

                C_Main_B2(:,(C_delta0-1)*n^2+1:C_delta0*n^2)=C_Main_B2(:,(C_delta0-1)*n^2+1:C_delta0*n^2) - kron(In',C0');
                C_Main_B2(:,(C_delta2-1)*n^2+1:C_delta2*n^2)=C_Main_B2(:,(C_delta2-1)*n^2+1:C_delta2*n^2) - kron(In',C1');
                C_Main_B2(:,(C_delta3-1)*n^2+1:C_delta3*n^2)=C_Main_B2(:,(C_delta3-1)*n^2+1:C_delta3*n^2) - kron(In',C2');


                % Matriz de ceros para C:
                C=zeros(2*m*n^2,(2*m-1) * n^2);

                % Primer Bloque de la Matriz C:
                for i=1:m
                    C((i-1)*n^2+1:i*n^2,(i-1)*n^2+1:(m-1+i)*n^2)=C_Main_B1;
                end
                % Segundo Bloque de la Matriz L:
                for i=m:2*m-1
                    C(i*n^2+1:(i+1)*n^2,(i-m)*n^2+1:i*n^2)=C_Main_B2;
                end

            end

            for MATRIZ_D=1:1 % MATRIZ D.
                % Generación de la matriz D:
                D_main_B1 = [Inn -Inn];
                D = zeros((2*m-1)*n^2,2*m*n^2);
                for i=1:2*m-1
                    D( (i-1)*n^2+1:i*n^2, (i-1)*n^2+1:(i+1)*n^2 ) = D_main_B1;
                end
            end

            for MATRIZ_L=1:1 % MATRIZ L.
                % Generación de la matriz L.
                L = [ A              C;
                    D   zeros( (2*m-1)*n^2 )];
            end
            
            for MATRIZ_M=1:1 % MATRIZ M.
                % Generación de la matriz M:
                M = zeros((4*m-1)*n^2);
                M_main_B1 = eye((2*m-1)*n^2);
                M_main_B2 = eye((2*m-2)*n^2);

                M(1:(2*m-1)*n^2,1:(2*m-1)*n^2)                                           = M_main_B1;
                M((2*m-1)*n^2+1:(N_Dinamicas-2)*n^2,2*m*n^2+1:(N_Dinamicas-1)*n^2)       = M_main_B2;
                M((N_Dinamicas-2)*n^2+1:(N_Dinamicas-1)*n^2,(3*m-1)*n^2+1:3*m*n^2)       = Inn;
                M((N_Dinamicas-1)*n^2+1: N_Dinamicas*n^2,(m-1)*n^2+1:2*m*n^2 )           = A_Main_B1;
                M((N_Dinamicas-1)*n^2+1: N_Dinamicas*n^2,(3*m-1)*n^2+1:N_Dinamicas*n^2)  = C_Main_B1;
            end
            
            for MATRIZ_N=1:1 % MATRIZ N.
                % Generación de la matriz N:
                N = zeros((4*m-1)*n^2);
                N_main_B1 = M_main_B1;
                N_main_B2 = M_main_B2;
                N(1:(2*m-1)*n^2,n^2+1:2*m*n^2)                                           = -N_main_B1;
                N((2*m-1)*n^2+1:(N_Dinamicas-2)*n^2,(2*m+1)*n^2+1:N_Dinamicas*n^2)       = -N_main_B2;
                N((N_Dinamicas-1)*n^2+1: (N_Dinamicas)*n^2,1:(m+1)*n^2)                  = -A_Main_B2;
                N((N_Dinamicas-1)*n^2+1: (N_Dinamicas)*n^2,2*m*n^2+1:3*m*n^2)            = -C_Main_B2;
            end

            for MATRIZ_E=1:1 % MATRICES E1 Y E2
                % MATRIZ E1 posición 4m-2:
                E1 = zeros(n^2,(N_Dinamicas)*n^2);
                E1(:,(4*m-3)*n^2+1:(4*m-2)*n^2 ) = Inn;

                % MATRIZ E2 posición m+1:
                E2 = zeros(n^2,(N_Dinamicas)*n^2);
                E2(:,m*n^2+1:(m+1)*n^2 ) = Inn;
            end

            for MATRIZ_X=1:1 % MATRIZ X
                % Generación de la matriz X
                X1  = [zeros(n^2,n^2)             E2;
                    zeros(N_Dinamicas*n^2,n^2) L  ];
                I_0 = [zeros(n^2,N_Dinamicas*n^2); eye(N_Dinamicas*n^2)];
                X   = M + [-E1' N]*expm(X1*h_base)*I_0;
            end

            
            for COND_FRONT=1:1 %CONDICIONES DE FRONTERA.
                I_Ceros=zeros(1,(N_Dinamicas-1)*n^2);
                W=eye(n);
                x0 = -X\[vec(I_Ceros); vec(W)];
            end

            %disp('Si');
            
            for MAT_KHARITONOV=1:1 %MATRIZ DE KHARITONOV
                
                %Definir parametros iniciales P1 y P2:
                p1=r-1;
                p2=m;

                %Número de elementos repetidos y con el número base.
                h_k=gcd(p1,p2); 

                %Cuantos numeros diferentes existen:
                num_difere=p1/h_k;

                %Generar el vector de argumentos:
                for i=1:r-1
                argumentos(i)=mod(m*(i-1),(r-1));
                end
                                
                %R inicia con las condiciones iniciales;
                R_ac=x0;
                h_ac=0;
                R=expm(L*h_k/(r-1)*h_base);
                
                cont_exp=1;
                aux=1;

                %disp('Si 2')
                
                for i=1:num_difere
                    contador=1;
                        for v=1:r-1
                            if h_ac==argumentos(v)
                                %Calcular bloque en donde se almacena.
                                bloque=m-floor((v-1)*m/(r-1));

                                %Extraer de R la información y coloca en Kr de acuerdo al
                                %índice detectado:
                                V1(:,(v-1)*n+1:(v-1+1)*n)=mat(R_ac((bloque-1)*n^2+1:bloque*n^2)); 

                                %Al encontrar todos los valores, romper el for más cercano.
                                if contador==h_k
                                    break
                                end
                                contador=contador+1;

                            end

                        end
    
                        if cont_exp==val_exp
                            R_ac=expm(L*h_k/(r-1)*h_base*val_exp*aux)*x0;
                            aux=aux+1;
                            cont_exp=0;
                        else
                            R_ac=R*R_ac;
                        end
                        
                        cont_exp=cont_exp+1;
                        h_ac=h_k+h_ac;

  
                end
               
                %Agregar la U(H)
                V1(:,(r-1)*n+1:r*n)=mat(R_ac(1:n^2));

                %Generación de la matriz de Kharitonov:
                Kr=zeros(n*r); 
                for i=1:r
                    Kr( n*(i-1)+1:n*i , n*(i-1)+1:r*n ) = V1(:,1:(r-i+1)*n);
                end
                %disp('Si 3')

                for i=2:r
                Kr(n*(i-1)+1:n*r,n*(i-2)+1:n*(i-1))=(V1(:,n+1:n*(r-i+2)))';
                end
                
                
            end

            for CALCULO_r_ESTRELLA = 1:1  %CALCULO DE r ESTRELLA
                % Obtención de la matriz P:
                sigma   = 2;
                W_barra = [W - W/sigma     zeros(n);
                             zeros(n)         W];
                A_barra = [A0+A0'-H*eye(n) A1;
                           A1'           zeros(n)];
        
                P = inv(W_barra)*A_barra;
        
                % Cálculo de alpha0_1:
                alpha0_1 = -1/( (m+1)*min( eig(P) ));
        
                %Cálculo de P_barra:
                P_barra = inv(W)*(G'*G);
                
                % Cálculo de alpha0_2
                alpha0_2 = 1/( sigma*H*(m+1)*max(eig(P_barra))  );
        
                % Verificación de resultados: R1 y R2:
                alpha0 = min(alpha0_1,alpha0_2)-min(alpha0_1,alpha0_2)/1000;
            end

            for MATRIZ_FUNDAMENTAL=1:1 % CÁLCULO DE LA MATRIZ FUNDAMENTAL
                Pr = In;
                for i=1:r-1
                    t=i*H/(r-1);

                    % Método paso a paso:

                    if(t<=h1)
                    Pr(:,n*i+1:(i+1)*n) = In;
                    end

                    if (t>h1 && t<=2*h1)
                        
                        Pr(:,n*i+1:(i+1)*n) = In +(A1 - h1*G)*( t - h1 ) + 1/2 * G * (t^2 - h1^2);
                    end

                    if (t>2*h1 && t<=3*h1)

                        q0 = A1 - 2*h1*A1^2 + 2*h1^2*A1*G - h1*G + 2*h1^2*G*A1 - 4/3*h1^3*G^2;
                        q1 = A1^2 - 2*h1*A1*G + G - 2*h1*G*A1 + 2*h1^2*G^2;
                        q2 = 1/2*A1*G + 1/2*G*A1 - h1*G^2;
                        q3 = 1/6*G^2;

                        Pr(:,n*i+1:(i+1)*n) =  In + (A1-h1*G)*h1 + 3/2*h1^2*G + q0*(t-2*h1) + 1/2*q1*(t^2 -4*h1^2)+1/3*q2*(t^3 -8*h1^3) + 1/4 * q3*(t^4 -16*h1^4);
                    end

                end
            end
              

                        
            for ESTABILIDAD=1:1 %TEOREMA DE ESTABILIDAD
                %Empleo de criterio de Estabilidad 
                        
                % % Primera forma de evaluar Kr -alpha0*PP':
                %Kr_eig=eig(Kr);
                %if(min(Kr_eig)>0)
                %   k1_vector(indice)=k1;
                %   k2_vector(indice)=k2;
                %   indice=indice+1;                       
                %end

                % Teorema Kr - alpha0*(P'P)>0
                [Rx_Th1, Px_Th1] = chol( Kr );

                
                if( Px_Th1 == 0)
                   % disp('Estable')
                   k1_vector_Th1(indice_Th1) = k1;
                   k2_vector_Th1(indice_Th1) = k2;
                   indice_Th1=indice_Th1+1;
                else
                   % disp('Inestable')
                end



                % Teorema Kr - alpha0*(P'P)>0
                [Rx_Th2, Px_Th2] = chol( Kr - alpha0*(Pr'* Pr) );

                
                if( Px_Th2 == 0)
                   % disp('Estable')
                   k1_vector_Th2(indice_Th2) = k1;
                   k2_vector_Th2(indice_Th2) = k2;
                   indice_Th2 = indice_Th2 + 1;
                else
                   % disp('Inestable')
                end

                
            
            end
                
        end
    end
end

%% RESULTADO GRÁFICO OBTENIDO DEL TEOREMA 1 Kr > 0:
figure(1)
subplot(1,2,1)
plot(k1_vector_Th1,k2_vector_Th1,... 
    '.',...
    'color',[1 0 0],...
    'MarkerSize',6,...
    'Marker','.')
grid on

title('Kr>0')
xlabel('k_{1}')
ylabel('k_{2}','Rotation',0)
grid on
axis([-10 30 -10 30])

%% RESULTADO GRÁFICO OBTENIDO DEL TEOREMA 2 Kr - a0 (P' * P) > 0:
subplot(1,2,2)
plot(k1_vector_Th2,k2_vector_Th2,... 
    '.',...
    'color',[0 0 0],...
    'MarkerSize',6,...
    'Marker','.')
grid on

title('Kr-\alpha_{0}P_{r}^{T}P_{r}>0')
xlabel('k_{1}')
ylabel('k_{2}','Rotation',0)
grid on
axis([-10 30 -10 30])


toc