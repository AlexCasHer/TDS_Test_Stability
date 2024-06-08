%% PROGRAM: Legendre approximation-based stability test for distributed delay systems 
% CASE: Distribuited delay with Constant Kernel.
% PDh Student: Alejandro Castaño Hernández
% Advisors: Dr. Mathieu Bajodek / Dra. Sabine Mondié

%% Clean the workspace
clear, clc, hold on, format long
tic

% Parameter ranges:
for PUNTOS=1:1  
   
        puntos = 160;
        epsilon = 1e-5;

        % Parameter p1: Example 1 [-100, 200] x [-100, 200]
        p1_min  = 0.01;
        p1_max  = 5;
        paso_p1 = (p1_max-p1_min)/puntos;

        % Parameter p2: Example 2 [0.1, 7] x [0.1, 7]
        p2_min  = 0.01;
        p2_max  = 5; 
        paso_p2 = (p2_max-p2_min)/puntos;
        
end

%Storage vectors
indice  = 1;
Estable = zeros(1,2);

digits(32);


% Definition of the order n:
n  = 1;  % Degree of the Legendre polynomial (polynomial approximation).
nD = 3;  % Number of dynamics of the semi-analytical method: U(tau),U(tau-h),int_{0}^{h}U(tau+theta)dtheta;

for p1 = p1_min:paso_p1:p1_max 
    for p2 = p2_min:paso_p2:p2_max

        % System parameters: 
        % x'(t) = A0x(t)+A1x(t-h)+G*int_{-h}^{0}x(t-\tau) d\tau
% %         h  = 0.1;
% %         A0 = p1;
% %         A1 = p2;
% %         G  = 1;
        h     = p1;
        gamma = p2;

        A0 = [ gamma-1  0   0       0;
                0       0   0      -5;
               -0.5556  0  -0.5556  0.5556;
                0       1  -1       0];

        A1 = zeros(length(A0));


        G = [-gamma/h  0  1/h 0;
                    0  0  0   0;
                    0  0  0   0;
                    0  0  0   0];

        % Initial arrays:
        nx = length(A0);
        In = eye(nx);
        I0 = zeros(nx);


        % ====================== DELAY LYAPUNOV MATRIX: CONSTANT KERNEL ======================
        % Matrix L:
        L = [ kron(A0',In)  kron(A1',In)  kron(G',In);
             -kron(In,A1') -kron(In,A0') -kron(In,G');
              kron(In,In)  -kron(In,In)   kron(I0,I0)];

        % Decomposition in Jordan form:
        [P,D] = eig(L);

        % Matrix M:
        M  = [ kron(In,In)   kron(I0,I0)   kron(I0,I0); 
               kron(I0,I0)   kron(I0,I0)   kron(In,In);
               kron(A0',In)  kron(A1',In)  kron(G',In)];

        % Matrix N:
        N  = [ kron(I0,I0)  -kron(In,In)   kron(I0,I0);
               kron(I0,I0)   kron(I0,I0)   kron(I0,I0);
               kron(In,A1')  kron(In,A0')  kron(In,G')];

        % Matrix E:
        E  = [kron(I0,I0) kron(In,In) kron(I0,I0)];

        % Matrix X:
        X  = M + [-E' N]*expm([kron(I0,I0) E; zeros(nD*nx^2,nx^2) L]*h)*[zeros(nx^2,nD*nx^2); eye(nD*nx^2)];
        
        % Initial Conditions of the Semi-Analytical Method:
        W  = eye(nx);
        x0 = X\[vec(I0); vec(I0); -vec(W)];

        % Calculo de U0:
        tau = 0;
        U   = [kron(In,In) kron(I0,I0) kron(I0,I0)]*expm(tau*L)*x0;
        U0  = double(reshape(U,nx,nx));


        % % Calculation of required inverses:
        Di = eye(length(L));
        for i=1:length(L)
            if  abs(D(i,i))  <= epsilon
                Di(i,i) = 0;
            elseif D(i,i) ~= 0
                Di(i,i) = 1/D(i,i);
            end
        end

        % =========================  STABILITY TEST =========================

        %% Qn and Rn terms:

        % % Initialization for Gamma:
        Gamma    = cell(n,1);
        Gamma{1} = Di*( expm(h*D) - eye(nD*nx^2) );                    %% Gamma[0]
        Gamma{2} = Di*( expm(h*D) + eye(nD*nx^2) ) - 2/h*Di*Gamma{1};  %% Gamma[1]

        % Initialization for Lambda:
        Lambda    = cell(n,1);
        Lambda{1} = Di*( Gamma{1} - h*eye(nD*nx^2) );                  %% Lambda[0]
        Lambda{2} = Di*( Gamma{2} );                                   %% Lambda[1]

        % Recursive relations for Gamma and Lambda: Nilpotent block
        for i = 1:length(L)
            if abs( D(i,i) ) <= epsilon
                % Initialization for Gamma0
                Gamma_aux      = zeros(nD*nx^2);
                Gamma_aux(i,i) = h;
                Gamma{1}       = Gamma{1} + Gamma_aux;

                % Initialization for Lambda0
                Lambda_aux      = zeros(nD*nx^2);
                Lambda_aux(i,i) = h^2/2;
                Lambda{1}       = Lambda{1} + Lambda_aux;

                % Initialization for Lambda1
                Lambda_aux      = zeros(nD*nx^2);
                Lambda_aux(i,i) = h^2/6;
                Lambda{2}       = Lambda{2} + Lambda_aux;
            end
        end
        
        % Recursive relations for Gamma and Lambda:
        for k=2:n-1
            Gamma{k+1}  = Gamma{k-1} - 2*(2*k-1)/h*Di*Gamma{k};
            Lambda{k+1} = Di*Gamma{k+1};
        end

        % Construction of the terms Qn and Rn:
        Qn = zeros(nx,n*nx);
        Rn = zeros(nx,n*nx);

        for k=0:n-1
            Qk = reshape([kron(In,In) kron(I0,I0) kron(I0,I0)]*P*Gamma{k+1}*(eye(3*nx^2)/P)*x0,nx,nx);
            Qn(:,k*nx+1:(k+1)*nx) = Qk'*A1;

            Rk = reshape([kron(In,In) kron(I0,I0) kron(I0,I0)]*P*Lambda{k+1}*(eye(3*nx^2)/P)*x0,nx,nx);
            Rn(:,k*nx+1:(k+1)*nx) = Rk'*G;
        end

        %% Term Tn:
        % Initialization for GammaBar
        GammaBar = cell(n);
        GammaBar{1,1} =  Di*( Gamma{1} - h*eye(nD*nx^2) );
        GammaBar{1,2} = -Di*Gamma{2};
        GammaBar{2,2} =  Di*( ( 2/h*Di - eye(nD*nx^2) )*Gamma{2} - h/3*eye(nD*nx^2) );

        % Recursive relations for GammaBar:
        for j=1:n
            for k=1:n
                if k >= max(3,j)
                    GammaBar{j,k} = GammaBar{j,k-2} + 2*(2*k-3)/h*Di*GammaBar{j,k-1};
                    % Addition of Kronecker Delta terms:
                    if k == j
                        GammaBar{j,k} = GammaBar{j,k} - h/(2*j-1)*Di;
                    end
                    % Addition of Kronecker Delta terms:
                    if k == j+2
                        GammaBar{j,k} = GammaBar{j,k} + h/(2*j-1)*Di;
                    end
                else
                    if k < j
                        GammaBar{j,k} = (-1)^(j+k)*GammaBar{k,j};
                    end
                end
            end
        end

        % Recursive relations for Gamma and Lambda: Nilpotent block
        for i=1:length(L)
            if abs( real(D(i,i)) ) <= epsilon
                GammaBar_aux      = zeros(nD*nx^2);
                GammaBar_aux(i,i) = h^2/2;
                GammaBar{1,1}     = GammaBar{1,1} + GammaBar_aux;
            end
        end

        %% Terms Sn:

        % Recursive relations for LambdaB=Zeta

        Zeta = cell(n);
        for j=1:n
            for k=1:n
                if k==1
                   Zeta{j,k} = h*Di*Gamma{j};

                   if j==1
                       Zeta{j,k} = Zeta{j,k} - h^2*Di*eye(nD*nx^2);
                   end

                else
                   Zeta{j,k} = zeros(nD*nx^2);
                end

            end
        end 

        % Recursive relations for  LambdaBF=ZetaBar
        ZetaBar = cell(n);
        for j=1:n
            for k=1:n
                ZetaBar{j,k} = ((-1)^(j-1))*Di*Gamma{k}*Gamma{j}*expm(-h*D);
                if j==1 && k==1
                    ZetaBar{j,k} = ZetaBar{j,k} - h^2*Di*eye(nD*nx^2);
                end
            end
        end

        % Recursive relations for  LambdaBF  and LambdaBF: Nilpotent block
        for i = 1:length(L)
            if abs( real(D(i,i)) ) <= epsilon
                % Inialización de Zeta{1,1}
                Zeta_aux      = zeros(nD*nx^2);
                Zeta_aux(i,i) = h^3/2;
                Zeta{1,1}     = Zeta{1,1} + Zeta_aux;

                if n>=2
                    % Initialization for Zeta{2,1}
                    Zeta_aux      = zeros(nD*nx^2);
                    Zeta_aux(i,i) = h^3/6;
                    Zeta{2,1}     = Zeta{2,1} + Zeta_aux;

                    % Initialization for ZetaBar{1,2}
                    ZetaBar_aux      = zeros(nD*nx^2);
                    ZetaBar_aux(i,i) = h^3/6;
                    ZetaBar{1,2}     = ZetaBar{1,2} + ZetaBar_aux;

                    % Initialization for{2,1}
                    ZetaBar_aux      = zeros(nD*nx^2);
                    ZetaBar_aux(i,i) = -h^3/6;
                    ZetaBar{2,1}     = ZetaBar{2,1} + ZetaBar_aux;
                end

            end
        end



        %% Terms D

        % Recursive relations for DD
        DD = cell(n);
        for j=1:n
            for k=1:n
                if k==1
                    DD{j,k} = h*Di^2*Gamma{j};
                    
                    if j==1
                        DD{j,k} = -h^3/2*Di - h^2*Di^2 + h*Di^2*Gamma{j};
                    end

                    if j==2
                        DD{j,k} = -h^3/6*Di^2+ h*Di^2*Gamma{j};
                    end

                else
                    DD{j,k} = zeros(nD*nx^2);
                end


            end
        end

        % Recursive relations for DBar:
        DBar = cell(n);

        for j=1:n
            for k=1:n
                DBar{j,k} = -(-1)^(j-1)*Di^2*Gamma{k}*Gamma{j}*expm(-h*D);

                if j==1
                   DBar{j,k} = DBar{j,k} + Di^2*h*Gamma{k};
                end

                if j==1 && k==1
                   DBar{j,k} = DBar{j,k} - h^3/2*Di;
                end

                if j==1 && k==2
                   DBar{j,k} = DBar{j,k} - h^3/6*Di;
                end
                
            end
        end

        % % Recursive relations for D and DBar: Nilpotent block
        for i = 1:length(L)
            if abs( real(D(i,i)) ) <= epsilon
                % Inialización de D{1,1}
                D_aux       = zeros(nD*nx^2);
                D_aux(i,i)  = h^4/6;
                DD{1,1}     = DD{1,1} + D_aux;

                if n>=2
                    % Initialization for D{2,1}
                    D_aux      = zeros(nD*nx^2);
                    D_aux(i,i) = h^4/12;
                    DD{2,1}    = DD{2,1} + D_aux;
                end

                if n>=3
                    % Initialization for D{3,1}
                    D_aux      = zeros(nD*nx^2);
                    D_aux(i,i) = h^4/60;
                    DD{3,1}     = DD{3,1} + D_aux;
                end


                % Initialization for DBar{1,1}
                DBar_aux      = zeros(nD*nx^2);
                DBar_aux(i,i) = h^4/12;
                DBar{1,1}     = DBar{1,1} + DBar_aux;

                if n>=2
                    % Initialization for DBar{1,2}
                    DBar_aux      = zeros(nD*nx^2);
                    DBar_aux(i,i) = h^4/12;
                    DBar{1,2}     = DBar{1,2} + DBar_aux;

                    % Initialization for DBar{2,2}
                    DBar_aux      = zeros(nD*nx^2);
                    DBar_aux(i,i) = h^4/36;
                    DBar{2,2}     = DBar{2,2} + DBar_aux;
                end

                if n>=3
                    % Initialization for DBar{3,1}
                    DBar_aux      = zeros(nD*nx^2);
                    DBar_aux(i,i) = -h^4/60;
                    DBar{3,1}     = DBar{3,1} + DBar_aux;
                end




            end
        end

        %% Terms Tn, Sn and Dn:
        Sn = zeros(n*nx);
        Tn = zeros(n*nx);
        Dn = zeros(n*nx);
        for j = 0:n-1
            for k = 0:n-1
                Tjk = reshape([kron(In,In) kron(I0,I0) kron(I0,I0)]*P*GammaBar{j+1,k+1}*(eye(3*nx^2)/P)*x0,nx,nx);
                Tn(j*nx+1:(j+1)*nx,k*nx+1:(k+1)*nx) = A1'*(Tjk+(-1)^(j+k)*Tjk')*A1;

                Sjk    = reshape([kron(In,In) kron(I0,I0) kron(I0,I0)]*P*Zeta{j+1,k+1}*(eye(3*nx^2)/P)*x0,nx,nx);
                SjkBar = reshape([kron(In,In) kron(I0,I0) kron(I0,I0)]*P*ZetaBar{j+1,k+1}*(eye(3*nx^2)/P)*x0,nx,nx);
                Sn(j*nx+1:(j+1)*nx,k*nx+1:(k+1)*nx) = A1'*(Sjk+SjkBar')*G;

                Djk    = reshape([kron(In,In) kron(I0,I0) kron(I0,I0)]*P*DD{j+1,k+1}*(eye(3*nx^2)/P)*x0,nx,nx);
                DjkBar = reshape([kron(In,In) kron(I0,I0) kron(I0,I0)]*P*DBar{j+1,k+1}*(eye(3*nx^2)/P)*x0,nx,nx);
                Dn(j*nx+1:(j+1)*nx,k*nx+1:(k+1)*nx) = G'*(Djk+DjkBar')*G;

            end
        end


        %% Termino IIn
        IIn = 1/h*kron(diag(2*(0:n-1)+1),eye(nx));

        % Matrix Pn: Stability Condition
        Pn = [U0 Qn+Rn; Qn'+Rn' Tn+Sn+Sn'+1/2*(Dn+Dn')+inv(IIn)];
        
        % Positivity test: Pn>0
        [Cholx,Px] = chol( double(Pn) );

        if( Px == 0)
            Estable(indice,:) = [p1;p2];
            indice = indice + 1;
        end


  
    end
end
plot(Estable(:,1),Estable(:,2),...
    '.',...
    'color',[0.9 0.9 0.1],...
    'MarkerSize',6,...
    'Marker','.')
axis square
grid on


toc