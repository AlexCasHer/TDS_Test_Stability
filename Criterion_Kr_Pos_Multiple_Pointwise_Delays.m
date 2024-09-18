%% PROGRAM: Finite Lyapunov Stability Test for a class of Linear Systems with Pointwise Delays
% Author/Developer: Alejandro Castaño Hernández

%% Clean the workspace 
clc, clear, hold on, format long
tic

% Parameter ranges:
for PARAMETERS=1:1
    
        % Initial parameters:
        points  = 120;
        epsilon = 5e-4;

        % Parameter p1:
        p1_min  = 0;
        p1_max  = 30;
        step_p1 = (p1_max - p1_min)/points;

        % Parameter p2:
        p2_min  = 0;
        p2_max  = 30; 
        step_p2 = (p2_max - p2_min)/points;
end

% Storage vectors
for ALMACENAMIENTO=1:1
    p1_vector = zeros();
    p2_vector = zeros();
    indice    = 1;
end

% Degree of the approximation: order r
r       = 10;            
val_exp = 10;


%% MAIN PROGRAM FOR STABILITY TEST.
for p1 = p1_min:step_p1:p1_max
    for p2 = p2_min:step_p2:p2_max

        % System parameters: x'(t) = A0x(t-h0)+A1x(t-h1)+A2x(t-h2)
        h0 = 0;
        h1 = p1;
        h2 = p2;


        A0 = -1.3;
        A1 = -1;
        A2 = -0.5;

        
        if (h1==0 && h2==0)
        else

            % ====================== DELAY LYAPUNOV MATRIX: MULTIPLE POINTWISE DELAYS ======================

            % Required variables:
            n  = size(A0,1);
            hb = gcd(100000*h1,100000*h2)/100000;
            H  = max(h1,h2);
            m  = round(H/hb);
            N_dynamics = round(2*m); 

            % Initial arrays:
            In = eye(n);
            W  = eye(n);
   
            % ========== Matrix L ==========

            % Generation of block 1 of Matrix L:
            L_Main_B1=zeros(n^2,(m+1)*n^2);

            delta0 = round(h0/hb);
            delta1 = round(h1/hb);
            delta2 = round(h2/hb);

            L_Main_B1(:,delta0*n^2+1:(delta0+1)*n^2) = L_Main_B1(:,delta0*n^2+1:(delta0+1)*n^2) + kron(A0',In);
            L_Main_B1(:,delta1*n^2+1:(delta1+1)*n^2) = L_Main_B1(:,delta1*n^2+1:(delta1+1)*n^2) + kron(A1',In);
            L_Main_B1(:,delta2*n^2+1:(delta2+1)*n^2) = L_Main_B1(:,delta2*n^2+1:(delta2+1)*n^2) + kron(A2',In);

            % Generation of block 2 of Matrix L:
            L_Main_B2 = zeros(n^2,(m+1)*n^2);

            delta0 = round((H-h0)/hb);
            delta1 = round((H-h1)/hb);
            delta2 = round((H-h2)/hb);

            L_Main_B2(:,delta0*n^2+1:(delta0+1)*n^2) = L_Main_B2(:,delta0*n^2+1:(delta0+1)*n^2) - kron(In',A0');
            L_Main_B2(:,delta1*n^2+1:(delta1+1)*n^2) = L_Main_B2(:,delta1*n^2+1:(delta1+1)*n^2) - kron(In',A1');
            L_Main_B2(:,delta2*n^2+1:(delta2+1)*n^2) = L_Main_B2(:,delta2*n^2+1:(delta2+1)*n^2) - kron(In',A2');

            % Initialization of matrix L:
            L = zeros(N_dynamics*n^2);

            % First Block of Matrix L:
            for i=1:m
                L((i-1)*n^2+1:i*n^2,(i-1)*n^2+1:(m+i)*n^2) = L_Main_B1;
            end

            % Second Block of Matrix L:
            for i=m:2*m-1
                L(i*n^2+1:(i+1)*n^2,(i-m)*n^2+1:(i+1)*n^2) = L_Main_B2;
            end

            % ========== Matrix M ==========

            % Initialization of matrix M:
            M=eye(N_dynamics*n^2);

            % Construction of the matrix M:
            M(n^2*(2*m-1)+1:2*m*n^2,:)=L(n^2*(m-1)+1:m*n^2,:);

            % ========== Matrix N ==========
            % Initialization of matrix M:
            N = zeros(N_dynamics*n^2);

            % Construction of the matrix M:
            N(1:n^2*(2*m-1),n^2+1:2*m*n^2) = -eye((N_dynamics-1)*n^2);
            N(n^2*(2*m-1)+1:2*m*n^2,:)     = -L(n^2*m+1:n^2*(m+1),:);


            % ========== Initial Conditions of the Semi-Analytical Method ==========  
            X  = M + N *expm(L*hb);
            w  = -[vec(zeros(1,(N_dynamics-1)*n^2)); vec(W)];
            X0 = linsolve(X,w);


            % =====================  EFFICIENT ALGORITHM FOR THE CONSTRUCTION OF MATRIX Kr =====================

            for MAT_KHARITONOV=1:1 %MATRIZ DE KHARITONOV
                
                % Define initial parameters pp1 and pp2:
                pp1 = r-1;
                pp2 = m;

                % Number of elements repeated and with the base number.
                h_k = gcd(pp1,pp2); 

                % Number of different elements existing:
                num_difere = pp1/h_k;

                % Generate the argument vector:
                argumentos = zeros();
                for i=1:r-1
                argumentos(i) = mod(m*(i-1),(r-1));
                end
                                
                % R starts with the initial conditions:
                R_ac = X0;
                h_ac = 0;
                R    = expm(L*h_k/(r-1)*hb);
                
                cont_exp = 1;
                aux      = 1;
                
                for i=1:num_difere
                    contador=1;
                        for v=1:r-1
                            if h_ac==argumentos(v)

                                % Calculate block where it is stored.
                                bloque = m-floor((v-1)*m/(r-1));

                                % Extract the information from R and place it in Kr according to the detected index
                                V1(:,(v-1)*n+1:(v-1+1)*n) = reshape(R_ac((bloque-1)*n^2+1:bloque*n^2),n,n); 

                                % When you find all the values, break the nearest for.
                                if contador == h_k
                                    break
                                end
                                contador = contador+1;

                            end

                        end
    
                        if cont_exp == val_exp
                            R_ac     = expm(L*h_k/(r-1)*hb*val_exp*aux)*X0;
                            aux      = aux+1;
                            cont_exp = 0;
                        else
                            R_ac = R*R_ac;
                        end
                        
                        cont_exp = cont_exp + 1;
                        h_ac     = h_k + h_ac;

  
                end
               
                % Add the U(H):
                V1(:,(r-1)*n+1:r*n) = reshape(R_ac(1:n^2),n,n);

                % Construction of the Kharitonov matrix:
                Kr = zeros(n*r); 
                for i=1:r
                    Kr( n*(i-1)+1:n*i , n*(i-1)+1:r*n ) = V1(:,1:(r-i+1)*n);
                end

                for i=2:r
                    Kr(n*(i-1)+1:n*r,n*(i-2)+1:n*(i-1)) = (V1(:,n+1:n*(r-i+2)))';
                end
                   
            end

            % =====================  STABILITY CONDITION =====================
            % Positivity test: Kr>0
            [Rx, Px] = chol(Kr);
    
            if(Px == 0)
                p1_vector(indice) = p1;
                p2_vector(indice) = p2;
                indice            = indice + 1;
    
            end

           
                
        end
    end
end


% Visualization of the stability map:
plot(p1_vector,p2_vector,... 
    '.',...
    'color',[0 0 1],...
    'MarkerSize',6,...
    'Marker','.')
xlabel('$p_{1}$','fontsize',14,'interpreter','latex')
ylabel('$p_{2}$','fontsize',14,'interpreter','latex','Rotation',0)

grid on
box on
axis([p1_min p1_max p2_min p2_max])

toc