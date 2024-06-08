%% PROGRAM: Legendre approximation-based stability test for distributed delay systems 
% CASE: Distribuited delay with Exponential Kernel.
% PDh Student: Alejandro Castaño Hernández
% Advisors: Dr. Mathieu Bajodek / Dra. Sabine Mondié

%% Clean the workspace
clear, clc, hold on, format long
tic

% Parameter ranges:
for PARAMETERS=1:1
    puntos = 80;
    
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
for STORAGE_VECTORS=1:1
    p1_vector = zeros();
    p2_vector = zeros();
    indice    = 1;
end

digits(32);

% Definition of the order n:
n  = 4;      % Degree of the Legendre polynomial (polynomial approximation).  

%% PROGRAMA PRINCIPAL:
for p1 = p1_min:paso_p1:p1_max
    for p2 = p2_min:paso_p2:p2_max

        % System parameters:
        % x'(t) = A0x(t)+A1x(t-h)+B1*int_{-h}^{0} exp(\gamma \tau) x(t-\tau) d\tau
        h = 0.1;                   % Delay
        A0 = eye(2);               % Matrix A0
        A1 = -[0 2.5; 2.5 0];      % Matrix A1
        B1 =  [p1 0; 0 p2];        % Matrix B1
        g  = -2.5;                 % Parameter \gamma
        
        % ====================== DELAY LYAPUNOV MATRIX: EXPONENTIAL KERNEL ======================

        % Required variables
        A = g;
        B = B1;
        W = eye(size(A0,1));
        e0 = expm(A*0)*1;
        eh = expm(-A*h)*e0;
        
        % Initial arrays:
        nx  = size(A0,1);
        m   = size(A,1);
        In  = eye(nx);
        ni  = 1:nx;
        n2i = 1:nx^2;
        d   = (2 * m + 2)*nx^2;
        
        % Kronecker products
        xA0 = kron(A0',eye(nx)); 
        A0x = kron(eye(nx),A0'); 
        xA1 = kron(A1',eye(nx));
        A1x = kron(eye(nx),A1');
        xB  = zeros(nx*fliplr(size(B)));  %The fliplr command returns B with its columns reversed in the left-right direction.
        Bx  = zeros(size(xB));
        
        for i=0:m-1
            s = i * nx;
            t = s * nx;
            xB(:, t + n2i) = kron(B(s + ni,:)', eye(nx));
            Bx(:, t + n2i) = kron(eye(nx), B(s+ni,:)');
        end
        
        % Matrix L:
        L = [  xA0,                xA1,                 xB,                       zeros(nx^2, m * nx^2);
              -A1x,               -A0x,                 zeros(nx^2, m * nx^2),     -Bx;
              kron(e0, eye(nx^2)), -kron(eh, eye(nx^2)), -kron(A, eye(nx^2)),        zeros(m * nx^2, m * nx^2);
              kron(eh, eye(nx^2)), -kron(e0, eye(nx^2)),  zeros(m * nx^2, m * nx^2),  kron(A,eye(nx^2))];
        
        % Vector E:
        E = [zeros(nx^2), eye(nx^2), zeros(nx^2, m * nx^2), zeros(nx^2, m * nx^2) ];
        
        % Matrices I1, I2:
        I1 = [eye(m*nx^2), zeros(m*nx^2, m*d)] * expm(-[zeros(m*nx^2), -kron(eye(m), E); zeros(m*d, m*nx^2), kron(A, eye(d)) - kron(eye(m), L)]*h) * [zeros(m*nx^2, d); kron(e0, eye(d))];
        I2 = [eye(m*nx^2), zeros(m*nx^2, m*d)] * expm( [zeros(m*nx^2),  kron(eye(m), E); zeros(m*d, m*nx^2), kron(A, eye(d)) + kron(eye(m), L)]*h) * [zeros(m*nx^2, d); kron(eh, eye(d))];
        
        % Matrix M:
        M = [eye(nx^2), zeros(nx^2), zeros(nx^2, m*nx^2),   zeros(nx^2, m * nx^2);
            [zeros(m * nx^2, nx^2),  zeros(m * nx^2, nx^2), zeros(m * nx^2), eye(m * nx^2)]   - I1;
            [zeros(m * nx^2, nx^2),  zeros(m * nx^2, nx^2), eye(m * nx^2),   zeros(m * nx^2)] - I2;
            xA0,                   xA1,                 xB,             zeros(nx^2, m * nx^2)];
        
        % Matrix N:
        N = [zeros(nx^2),     -eye(nx^2),    zeros(nx^2, m * nx^2),     zeros(nx^2, m * nx^2);
             zeros(m * nx^2, d);
             zeros(m * nx^2, d);
            A1x,    A0x,   zeros(nx^2, m * nx^2),   Bx];
        
        % Matrix X
        X = M + N * expm(L*h);
        w = -[zeros((2 * m + 1 )*nx^2,1); W(:)];
        
        % Initial Conditions of the Semi-Analytical Method
        Uvec = linsolve(X,w);

        % =====================  STABILITY TEST =====================

        % Auxiliary variables:
        E1 = [eye(nx^2), zeros(nx^2, (2 * m + 1)*nx^2)];
        U0 = reshape(E1*Uvec,nx,nx);

        % Calculation of required inverses:
        Li = L\eye(length(L));               
        gi = 1/g;                            
        IL = eye(length(L));
        gIL_Minus_L_in = (g*IL-L)\eye(length(L));
        gIL_Plus_L_in  = (g*IL+L)\eye(length(L));
        
        %% Qn and Rn terms:
      
        % Initialization for Gamma:
        Gamma    = cell(n,1);
        Gamma{1} = Li*( expm(h*L) - IL );                    %% Gamma[0]
        Gamma{2} = Li*( expm(h*L) + IL ) - 2/h*Li*Gamma{1};  %% Gamma[1]

        % Initialization for Gammab_k:
        Gammab    = cell(n,1);
        Gammab{1} = gi*IL*( IL - expm(-g*h*IL) );                        %% Gammab[0]
        Gammab{2} = 2*gi*IL*expm(-g*h*IL) - (2/h*gi*IL - IL)*Gammab{1};  %% Gammab[1]  

        % Initialization for Lambdak
        Lambda    = cell(n,1);
        Lambda{1} = gIL_Minus_L_in*(Gammab{1} - Gamma{1}*expm(-g*h*IL));   %% Lambda[0]
        Lambda{2} = gIL_Minus_L_in*(Gammab{2} - Gamma{2}*expm(-g*h*IL));   %% Lambda[1]

        % Recursive relations for Gamma, Lambdab and Lambda:
        for k=2:n-1
            Gamma{k+1}  = Gamma{k-1}  - 2*(2*k-1)/h*Li*Gamma{k};

            Gammab{k+1} = Gammab{k-1} - 2*(2*k-1)/h*gi*Gammab{k};
            Lambda{k+1} = gIL_Minus_L_in*( Gammab{k+1} - Gamma{k+1}*expm(-g*h*IL) );
        end

        % Construction of the terms Qn and Rn:
        Qn = zeros(nx,n*nx);
        Rn = zeros(nx,n*nx);

        for k=0:n-1
            Qk = mat(E1*Gamma{k+1}*Uvec);
            Qn(:,k*nx+1:(k+1)*nx) = Qk'*A1;

            Rk = mat(E1*Lambda{k+1}*Uvec);
            Rn(:,k*nx+1:(k+1)*nx) = Rk'*B1;
        end

        %% Term Tn:
        % Initialization for GammaBar
        GammaBar = cell(n);
        GammaBar{1,1} =  Li*( Gamma{1} - h*IL );                        %% GammaBar[0,0]
        GammaBar{1,2} = -Li*Gamma{2};                                   %% GammaBar[0,1]
        GammaBar{2,2} =  Li*( ( 2/h*Li - IL )*Gamma{2} - h/3*IL );      %% GammaBar[1,1]

        % Recursive relations for GammaBar
        for j=1:n
            for k=1:n
                if k >= max(3,j)
                    GammaBar{j,k} = GammaBar{j,k-2} + 2*(2*k-3)/h*Li*GammaBar{j,k-1};
                    % Addition of Kronecker Delta terms:
                    if k == j
                        GammaBar{j,k} = GammaBar{j,k} - h/(2*j-1)*Li;
                    end
                    % Addition of Kronecker Delta terms:
                    if k == j+2
                        GammaBar{j,k} = GammaBar{j,k} + h/(2*j-1)*Li;
                    end
                else
                    if k < j
                        GammaBar{j,k} = (-1)^(j+k)*GammaBar{k,j};
                    end
                end
            end
        end

      


        %% Terms Tn, Sn and Dn:
        Tn = zeros(n*nx);
        Sn = zeros(n*nx);
        Dn = zeros(n*nx);

        %  Recursive relations and construction of the terms Tn, Sn and Dn:
        for j = 0:n-1
            for k = 0:n-1
                Tjk = reshape(E1*GammaBar{j+1,k+1}*Uvec,nx,nx);
                Tn(j*nx+1:(j+1)*nx,k*nx+1:(k+1)*nx) = A1'*(Tjk+(-1)^(j+k)*Tjk')*A1;

                Sjk  = reshape(E1*gIL_Plus_L_in*( Gamma{j+1} - (-1)^(j) * Gammab{j+1})*Gammab{k+1}*Uvec,nx,nx);
                SjkB = reshape(E1*gIL_Minus_L_in*(-1)^(j)*( Gammab{j+1}*Gammab{k+1} - Gamma{j+1}*Gamma{k+1}*expm(-h*(g*IL + L) ) ) *Uvec,nx,nx);
                Sn(j*nx+1:(j+1)*nx,k*nx+1:(k+1)*nx) = A1'*(Sjk + SjkB')*B1;

                Djk  = reshape( E1*gIL_Plus_L_in*(  1/2*(-1)^(j)*gi*Gammab{j+1}*exp(-g*h) - gIL_Minus_L_in*Gamma{j+1}*exp(-g*h) + gIL_Minus_L_in*Gammab{j+1} - 1/2*gi*Gammab{j+1}   ) * Gammab{k+1}*Uvec,nx,nx );
                DjkB = reshape( E1*gIL_Minus_L_in*( gIL_Plus_L_in*(-1)^(j)*Gamma{j+1}*Gamma{k+1}*expm(-h*L)*exp(-2*h*g) - gIL_Plus_L_in*Gammab{j+1}*Gamma{k+1}*exp(-g*h) - 1/2*gi*(-1)^(j)*Gammab{j+1}*Gammab{k+1}*exp(-g*h) + 1/2*gi*Gammab{j+1}*Gammab{k+1} )*Uvec,nx,nx);
                Dn(j*nx+1:(j+1)*nx,k*nx+1:(k+1)*nx) = B1'*(Djk+DjkB')*B1;
            end
        end

                


        %% Term IIn:
        IIn = 1/h*kron(diag(2*(0:n-1)+1),eye(nx));

        % Matrix Pn: Stability Condition
        Pn = [U0 Qn+Rn; Qn'+Rn' Tn+Sn+Sn'+Dn+inv(IIn)];

        % Positivity test: Pn>0
        [Rx, Px] = chol(Pn);

        if(Px == 0)
            p1_vector(indice) = p1;
            p2_vector(indice) = p2;
            indice            = indice + 1;
        end

    end
end
toc
plot(p1_vector,p2_vector,'.r')
axis([p1_min p1_max p2_min p2_max])
grid on
box on
hold on

