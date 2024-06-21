%% PROGRAM: Finite Lyapunov Stability Test for a class of Linear Systems with Pointwise and Distributed Delays 
% CASE: Distribuited delay with Exponential Kernel.
% PDh Student: Alejandro Castaño Hernández
% Advisors: Dr. Carlos Cuvas / Dra. Sabine Mondié / Dr. Alexey Egorov

%% Clean the workspace 
clear, clc, hold on

% Parameter ranges:
for PARAMETROS=1:1
    puntos = 180;

    % Parameter p1:
    p1_min  = -500;
    p1_max  =  500;
    paso_p1 = (p1_max-p1_min)/puntos;

    % Parameter p2:
    p2_min  = -500;
    p2_max  =  500;
    paso_p2 = (p2_max-p2_min)/puntos;
end

%Storage vectors
for ALMACENAMIENTO=1:1
    p1_vector = zeros();
    p2_vector = zeros();
    indice    = 1;
end

% Definition of the order r:
r=5;            % Degree of the approximation.
val_exp = 1;


for p1 = p1_min:paso_p1:p1_max
    for p2 = p2_min:paso_p2:p2_max
  
        % System parameters:
        % x'(t) = A0x(t)+A1x(t-h)+B1*int_{-h}^{0} exp(\gamma \tau) x(t-\tau) d\tau
        h = 0.1;
        A0 = eye(2);
        A1 = -[0 2.5; 2.5 0];
        B1 =  [p1 0; 0 p2];
        g  = -2.5;


        % ====================== DELAY LYAPUNOV MATRIX: EXPONENTIAL KERNEL ======================

        % Required variables:
        A = g;
        B = B1;
        W = eye(size(A0,1));
        e0 = expm(A*0)*1;
        eh = expm(-A*h)*e0;
        
        % Initial arrays:
        n   = size(A0,1);
        m   = size(A,1);
        In  = eye(n);
        ni  = 1:n;
        n2i = 1:n^2;
        d   = (2 * m + 2)*n^2;
        
        % Kronecker products:
        xA0 = kron(A0',eye(n)); 
        A0x = kron(eye(n),A0'); 
        xA1 = kron(A1',eye(n));
        A1x = kron(eye(n),A1');
        xB  = zeros(n*fliplr(size(B)));  % El comando fliplr devuelve B con sus columnas invertidas en la dirección izquierda-derecha.
        Bx  = zeros(size(xB));
        
        for i=0:m-1
            s = i * n;
            t = s * n;
            xB(:, t + n2i) = kron(B(s + ni,:)', eye(n));
            Bx(:, t + n2i) = kron(eye(n), B(s+ni,:)');
        end
        
        % Matrix L:
        L = [  xA0,                xA1,                 xB,                       zeros(n^2, m * n^2);
              -A1x,               -A0x,                 zeros(n^2, m * n^2),     -Bx;
              kron(e0, eye(n^2)), -kron(eh, eye(n^2)), -kron(A, eye(n^2)),        zeros(m * n^2, m * n^2);
              kron(eh, eye(n^2)), -kron(e0, eye(n^2)),  zeros(m * n^2, m * n^2),  kron(A,eye(n^2))];

        if rank(L)<4
            disp('Pierde rango')
        end
        
        % Vector E:
        E = [zeros(n^2), eye(n^2), zeros(n^2, m * n^2), zeros(n^2, m * n^2) ];
        
        % Matrices I1, I2:
        I1 = [eye(m*n^2), zeros(m*n^2, m*d)] * expm(-[zeros(m*n^2), -kron(eye(m), E); zeros(m*d, m*n^2), kron(A, eye(d)) - kron(eye(m), L)]*h) * [zeros(m*n^2, d); kron(e0, eye(d))];
        I2 = [eye(m*n^2), zeros(m*n^2, m*d)] * expm( [zeros(m*n^2),  kron(eye(m), E); zeros(m*d, m*n^2), kron(A, eye(d)) + kron(eye(m), L)]*h) * [zeros(m*n^2, d); kron(eh, eye(d))];
        
        % Matrix M:
        M = [eye(n^2), zeros(n^2), zeros(n^2, m*n^2),   zeros(n^2, m * n^2);
            [zeros(m * n^2, n^2),  zeros(m * n^2, n^2), zeros(m * n^2), eye(m * n^2)]   - I1;
            [zeros(m * n^2, n^2),  zeros(m * n^2, n^2), eye(m * n^2),   zeros(m * n^2)] - I2;
            xA0,                   xA1,                 xB,             zeros(n^2, m * n^2)];
        
        % Matrix N:
        N = [zeros(n^2),     -eye(n^2),    zeros(n^2, m * n^2),     zeros(n^2, m * n^2);
             zeros(m * n^2, d);
             zeros(m * n^2, d);
            A1x,    A0x,   zeros(n^2, m * n^2),   Bx];
        
        % Matrix X
        X = M + N * expm(L*h);
        w = -[zeros((2 * m + 1 )*n^2,1); W(:)];
        
        % Initial Conditions of the Semi-Analytical Method
        X0 = linsolve(X,w);

        % =====================  STABILITY TEST =====================

        for MAT_KHARITONOV=1:1 %MATRIZ DE KHARITONOV
                
                % Define initial parameters pp1 and pp2:
                H      = h;
                h_base = h;
                m_m    = round(H/h_base);
                pp1    = r-1;
                pp2    = m_m;

                % Number of elements repeated and with the base number.
                h_k = gcd(pp1,pp2); 

                % How many different numbers are there:
                num_difere = pp1/h_k;

                % Generate the argument vector:
                argumentos = zeros();
                for i=1:r-1
                argumentos(i)=mod(m_m*(i-1),(r-1));
                end
                                
                % R starts with the initial conditions:
                R_ac=X0;
                h_ac=0;
                R=expm(L*h_k/(r-1)*h_base);
                
                cont_exp=1;
                aux=1;
                
                for i=1:num_difere
                    contador=1;
                        for v=1:r-1
                            if h_ac==argumentos(v)

                                % Calculate block where it is stored.
                                bloque=m_m-floor((v-1)*m_m/(r-1));

                                % Extract the information from R and place it in Kr according to the detected index
                                V1(:,(v-1)*n+1:(v-1+1)*n)=reshape(R_ac((bloque-1)*n^2+1:bloque*n^2),n,n); 

                                % When you find all the values, break the nearest for.
                                if contador==h_k
                                    break
                                end
                                contador=contador+1;

                            end

                        end
    
                        if cont_exp==val_exp
                            R_ac=expm(L*h_k/(r-1)*h_base*val_exp*aux)*X0;
                            aux=aux+1;
                            cont_exp=0;
                        else
                            R_ac=R*R_ac;
                        end
                        
                        cont_exp=cont_exp+1;
                        h_ac=h_k+h_ac;

  
                end
               
                % Add the U(H)
                V1(:,(r-1)*n+1:r*n)=reshape(R_ac(1:n^2),n,n);

                % Construction of the Kharitonov matrix:
                Kr=zeros(n*r); 
                for i=1:r
                    Kr( n*(i-1)+1:n*i , n*(i-1)+1:r*n ) = V1(:,1:(r-i+1)*n);
                end

                for i=2:r
                Kr(n*(i-1)+1:n*r,n*(i-2)+1:n*(i-1))=(V1(:,n+1:n*(r-i+2)))';
                end
                   
        end


        % Stability Condition, Positivity test: Kr>0
        [Rx, Px] = chol(Kr);

        if(Px == 0)
            p1_vector(indice) = p1;
            p2_vector(indice) = p2;
            indice            = indice + 1;

        end


    end
end

plot(p1_vector,p2_vector,'.b')
axis([p1_min p1_max p2_min p2_max])
grid on
box on
axis([p1_min p1_max p2_min p2_max])