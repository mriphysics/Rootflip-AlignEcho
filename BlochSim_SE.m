function [Mxy,Mz] = BlochSim_SE(B1,dt,G,r,M0,T2,T1,Anime)
%BLOCHSIM Simulates the change of Magnetization in presence of Gradients
%and B1 pulse according to the Bloch Equation.
% 29/09/16 BlochSim_SE is specifically equipped for spin-echo simulations
% with T1/T2 decay.

%   B1 (Ntx1 Vector) holding the B1-pulse. Can be complex-valued.
%   G  (Ntx3 Matrix) holding the time-varying Gradient (x,y,z) values
%   r  (Nrx3 Matrix) the spatial coordinates where the simulation runs.
%       Currently only Nz works. Extension for the other two dimensions is
%       TBD
%   M0 (1x3 vector) Initial Magnetisation state [Mx My Mz]. Typically [0 0 1]
%   Anime (Boolean) Set to 1 if you want to store the magnetisation profile
%       over time. 
gyro=2*pi*42.58*1e3; %in rad/s/mT
B1x=real(B1);
B1y=imag(B1);
% M=M0*ones(1,length(r));
Nt=length(B1);
Nr=length(r);
if Anime 
    Mout = zeros(Nt,Nr,3);
else
    Mout = zeros(Nr,3);
end

% 12/02/2015
%For each time-step i, find for each location z, the change in
%magnetisation over a tiny time-instance dt. The general problem is given
%by M'= BM * M. This can be solved by a method of solving a system of
%differential equations. I used: http://tutorial.math.lamar.edu/Classes/DE/ComplexEigenvalues.aspx

% M = C * V * exp(D*dt) where:
% C = Initial conditions constants, which only apply for a dt time span
% V = 3x3 matrix where the columns are eigenvectors of BM
% D = eigenvalues of BM. 

% Parfor implementation:

parfor j =1:Nr
    M = M0*ones(1,Nt);
    BM= [-1/T2 gyro*G(1,:)*r(j,:)' gyro*-B1y(1);
            -gyro*G(1,:)*r(j,:)' -1/T2 gyro*B1x(1);
            gyro*B1y(1) -gyro*B1x(1) 0];
    [V,D]=eig(BM);
    C = M(:,1)'/(V');
    M(:,1)= C*diag(exp(diag(D)*dt))*V';
    for i=2:Nt

        BM0 = BM;
        BM= [-1/T2 gyro*G(i,:)*r(j,:)' gyro*-B1y(i);
            -gyro*G(i,:)*r(j,:)' -1/T2 gyro*B1x(i);
            gyro*B1y(i) -gyro*B1x(i) -1/T1];
        
        if ~isequal(BM,BM0)
            [V,D]=eig(BM);
        end
        C = M(:,i-1)'/(V');
        M(:,i)= C*diag(exp(diag(D)*dt))*V' + dt*[0 0 M0(3)/T1];
   
    end           
    Mout(:,j,:)=M(:,:)';

end

if Anime
    Mxy = Mout(:,:,1)+1i*Mout(:,:,2);
    Mz = Mout(:,:,3);
else
    Mxy = Mout(end,:,1) + 1i*Mout(end,:,2);
    Mz = Mout(end,:,3);
end
end