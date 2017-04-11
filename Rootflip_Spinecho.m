clear all;
close all;


%%%---------- User defined setting here  ----------------------------%%%

gamma_mT = 267522.1; %<--- radians/mT/s
slthick = 2e-3; %<-- slice thickness in mm

% Set excitation type.
% For Align-all, set exc_type to 'alignte'.

% exc_type = 'minduration';
% exc_type = 'alignedecho';
exc_type = 'alignte';

CS = 2; %<-- 0 for Time-symmetric (Sharma MRM 2016). Set to 1 for AM-only. Set to 2 for "all-flip" (No symmetry).

% Relaxation parameters:
T2 = 0.080;   %<--- T2 relaxation in s
T2s = 0.045; %<--- T2* relaxation in s
T1 = inf;

nt = 256;%<-- number of time-points
mb = 3; %<-- Multiband factor

% Note: I've empirically found that TB 2.5-4 leads to 2 roots per passband.
% This work was developed and tested for TB4 pulses. But other TB values
% should work as long as there are a even number of roots per passbands.

% For 4 roots per passband, I found TB 7.2 to 7.5 to work. 
% Configurations for 6 roots and above have not been attempted.

tb = 4; %<-- time-bandwidth product
bs = 3; %<-- slice sep in units of slice thickness
maxb1 = 0.020; %<-- mT

M = round(tb/2); %<-- this will need some manipulation when is not 4.

%%%------------------------------------------------------------------%%%
fprintf('Designing Align-Echo Rootflipped pulse with M=%.d\n',M);
[ ~,rf180,tb,bout] = rootflipAE(nt,mb,tb,bs,CS,M);

figure;
cscatter(roots(bout));
title('Roots plotted');
ylabel('Imaginary axis');
xlabel('Real axis');
RF.rf180 = rf180;

%% Find exictation pulse from the refocusing pulse for matched-excitation

[a90,b90di] = rf2ab_sas(RF.rf180,(-nt/2:1/2:nt/2-1/2)',0); % get beta of 180 pulse
b90d = (b90di.^2)/sqrt(2); % target 90-deg beta profile
b90d= -conj(b90d);
bx=fftshift( fft(ifftshift(b90d))/length(b90d) );
[~,ax]=b2amp(bx);
rf90 = -1i*conj(islr(ax,bx));
if CS==1
    rf90=real(rf90);
end

RF.rf90 = rf90;

RF.dt180 = max(abs(RF.rf180))/(gamma_mT*maxb1);
RF.T180 = length(RF.rf180)* RF.dt180;


switch exc_type
    case 'alignedecho'
        RF.dt90 = RF.dt180;
    case 'alignte'
        RF.dt90 = RF.dt180/2;
    case 'minduration'
        RF.dt90 = max(abs(RF.rf90))/(gamma_mT*maxb1);
    otherwise
        warning('Excitation method not specified! Minimum duration excitation chosen as default')
        RF.dt90 = max(abs(RF.rf90))/(gamma_mT*maxb1);
end
RF.T90 = length(RF.rf90)*RF.dt90;

% This step choses the shortest dwell-time for interpolation and simulation.
% If desired, you could set your own dwell-time. But make sure it's short
% enough to avoid sampling artefacts

dt_str={'dt90','dt180'};
[dt,ind]=min([RF.dt90 RF.dt180 ]); %Avoid Aliasising in simulation.

RF.tb = tb;

%    calculate slice bandwidths
RF.bw180 = RF.tb/RF.T180;
RF.bw90 = 2*RF.tb/RF.T90;

%% Generate Time-axis

Tsim = 150*1e-3; %<-- Total time of simulation
nominal_TE = 90*1e-3; %<-- Time in simulation when spin-echoes should form.
t1 = 0;
t2 = RF.T90; %<-- Excitation pulse duration
t4 = RF.T180;%<-- Refocusing pulse duration
t3 = (nominal_TE -(t2+t4))/2; %< -- time between end of Excitation and beginning of refocusing
t5 = Tsim - sum([t1 t2 t3 t4]); 
if t5<0
    error('Tsim too low');
end
tt = ceil([t1 t2 t3 t4 t5]/dt)*dt;
nn = ceil(tt/dt);
nnc = cumsum(nn);
%% Interpolate RF pulses and generate gradient waveforms.
intmethod = 'linear';
RF.rf90i =(interp1(linspace(0,1,length(RF.rf90)) ,RF.rf90 ,linspace(0,1,nn(2)),intmethod)/(gamma_mT*RF.dt90))';
RF.rf180i=(interp1(linspace(0,1,length(RF.rf180)),RF.rf180,linspace(0,1,nn(4)),intmethod)/(gamma_mT*RF.dt180))';

RF.g90  = 2*pi*RF.bw90 /(gamma_mT*slthick); %mT/m
RF.g180 = 2*pi*RF.bw180/(gamma_mT*slthick); %mT/m

RF.Gz90 = RF.g90*ones(length(RF.rf90),1);
RF.Gz180= RF.g180*ones(length(RF.rf180),1);

RF.G90 = [0*RF.Gz90 0*RF.Gz90 RF.Gz90];
RF.G180 = [0*RF.Gz180 0*RF.Gz180 RF.Gz180];

%% Evaluate optimal rewind gradient
% Define the space-grid on which the rewind
% gradient-moment is evaluated.
Nz = 4000; %<-- This needs to be high enough to correct in-phase slice errors in the optimal rewind search
FOVz=1*(mb*bs*slthick)+3*slthick;
z=linspace(-FOVz/2,FOVz/2,Nz)';
pos = [0*z(:) 0*z(:) z(:)];% in meters

% simulate Excitation pulse
[~,~,~,~,a90f,b90f] = blochsim_CK(RF.rf90(:)/(gamma_mT*RF.dt90),RF.G90,pos,ones(Nz,1),zeros(Nz,1),'dt',RF.dt90);
% simulate refocusing pulse
[~,~,~,~,~  ,b180f] = blochsim_CK(RF.rf180(:)/(gamma_mT*RF.dt180),RF.G180,pos,ones(Nz,1),zeros(Nz,1),'dt',RF.dt180);

% Evaluate the un-wound slice profile
Mxy_uw = 2*a90f(:,end).*conj(b90f(:,end)).*b180f(:,end).^2;

% optimal MB rewind code:
[idz,~] = idmxy(z,slthick,mb,bs);
inds = find(abs(diff(idz)) > 0.5);
inds(1:2:end) = inds(1:2:end) + 2;
inds(2:2:end) = inds(2:2:end) -1;
[phir,maxerr] = minimaxsmsphi(mb,inds,Mxy_uw,z);
fprintf('Rewind area:%.3f rad/m leads to error %.2f rad\n',phir,maxerr);
gphs = phir*z;
Mxy_rw = 2*a90f(:,end).*conj(b90f(:,end)).*conj(exp(1i*gphs)).*b180f(:,end).^2;

RF.RF = [
    zeros(nn(1),1);
    RF.rf90i;
    zeros(nn(3),1);
    RF.rf180i; 
    zeros(nn(5),1)];

%% Now that the rewind area is found, regenerate gradient pulses and re-evaluate slice profile 
grew = -phir/(gamma_mT*dt*nn(3)); %<-- find gradient rewind value
RF.Gz = [
    zeros(nn(1),1);
    RF.g90*ones(nn(2),1);
    grew*ones(nn(3),1);
    RF.g180*ones(nn(4),1);
    zeros(nn(5),1)];

RF.Gx = 0*RF.Gz;
RF.G = [RF.Gx(:),0*RF.Gz(:),RF.Gz(:)];


Nt = sum(nn);
t = (1:Nt)*dt*1000;
pos = [0*z(:) 0*z(:) z(:)];% in meters

exidc=1:nnc(2); %Excitation Indices
rwidc=nnc(2)+1:nnc(3); %Rewinding Indices
rfidc=(nnc(3)+1) : nnc(4); %Refocusing Indices

%         Combine Excitation and Rewinding
[~,~,~,~,RF.a90t,RF.b90t] =...
    blochsim_CK(RF.RF([exidc rwidc]),RF.G([exidc rwidc],:),pos,ones(Nz,1),zeros(Nz,1),'dt',dt);

[~,~,~,~,~,RF.b180t] = ...
    blochsim_CK(RF.RF(rfidc),RF.G(rfidc,:),pos,ones(Nz,1),zeros(Nz,1),'dt',dt);

RF.Mxy=2*RF.a90t.*conj(RF.b90t);
RF.Mxyexc = RF.Mxy;

RF.Mxy =[RF.Mxy,repmat(RF.Mxy(:,end),1,length(rfidc))];
RF.Mxy(:,nnc(3)+1:nnc(4)) = ...
    RF.Mxy(:,nnc(3)+1:nnc(4)).*RF.b180t.^2;

fh = figure;
set(fh, 'Position', [1, 1, 1397,791]);
subplot(221),plot(t,abs(RF.RF));axis([-inf inf -inf inf]);
title(sprintf('tb=%.2f, bs=%.2f',tb,bs));
subplot(222),imagesc(z,t,abs(RF.Mxy'));
colormap gray
title(sprintf('Exc Dur: %.4f ms Ref Dur: %.4f ms',RF.T90*1e3,RF.T180*1e3));
subplot(223),plotMxy(z,RF.Mxy(:,end));
subplot(224),cplot(RF.Mxy(:,end)); axis([-inf inf -inf inf]); 


%% Lorentzian-weighted dephasing simulation

% Simulate N_sc points per slice:
z_sim = [];
Nsc = 7; %<--- number of spatial points per slice.
spos = (1:mb)-(mb+1)/2; 
for i = 1:mb
    cp = spos(i)*slthick*bs;%<-- center point.
    z_sim = [z_sim linspace(cp-0.3*slthick,cp+0.3*slthick,Nsc)];
end
Nz_sim = length(z_sim);

% Make FOVx 1m
FOVx = 1; 
Nx = 200; %<-- number of spinors per spatial point. If there are not enough spinors, spin-echo "Aliases" will appear in time.
x = linspace(-FOVx/2,FOVx/2,Nx);

% Create dephasing gradient to simulate \Delta B0
% Note: I picked 50*1/T2s because I empirically found that it gives a
% sufficiently good exponential decay.
RF.gx = 2*pi*50*(1/T2s)/(gamma_mT*FOVx);
fx = gamma_mT*RF.gx*x/(2*pi);
theta = 1/(T2s);
mask = theta./(theta^2 + (fx).^2);   
% Scale mask such that integrated energy equals number of spatial points
mask = Nx*mask./sum(mask);

% plot(fx,mask); %<-- for interest.

% Keep off-resonance on all the time!
RF.Gx = RF.gx*ones(Nt,1);

RF.G = [RF.Gx(:),0*RF.Gz(:),RF.Gz(:)];

[xm,zm] = meshgrid(x,z_sim(:));
posm = [xm(:) 0*xm(:) zm(:)];

%% The heavy bit - uses parpool for acceleration
delete(gcp('nocreate')); %deletes current local pool if it exists.
myCluster=parcluster('local'); 
Workers=4;
myCluster.NumWorkers=Workers; 
parpool(myCluster,Workers);

tic;
[mxyt,~] = BlochSim_SE([RF.RF],dt,RF.G,posm,[0;0;1],T2,T1,1);
Evaltime = toc;
fprintf('Simulation took %.f seconds\n',Evaltime);
Mxydp = reshape(mxyt',[Nz_sim Nx Nt]); 
%%
Mxydp_L = repmat(mask,[Nz_sim 1 Nt]).*Mxydp;

RF.Mxydp = squeeze(sum(Mxydp_L,2));

dum = 1;
echo_str = zeros(mb,Nt);
for i=1:mb
    idc = ((i-1)*Nsc+1):i*Nsc;                            
    tmp_echo = abs(sum(RF.Mxydp(idc,:))/Nx/Nsc);                                
    echo_str(i,:) = abs(sum(RF.Mxydp(idc,:))/Nx/Nsc);
end
%% Display
t_rf = ((1:Nt)*dt - t2/2)*1e3;
figure;
subplot(2,1,1);
plot(t_rf,abs(RF.RF)*1e3);
xlabel('Time [ms]');
ylabel('|B_1| [\muT]');
subplot(2,1,2);
plot(t_rf,echo_str);
xlabel('Time [ms]');
ylabel('|M_{xy}(t)|');
