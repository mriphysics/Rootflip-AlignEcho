function plotMxy(z,Mxy)
% Plot Mxy excitation profile
%   Plot the abs(~) and unwrap(angle(~)) of a Excitation profile over
%   space/frequency axis z.

% 21/05/2015 sas
Nz=length(z);
% If z-axis not specified.. assign Mxy to the z-input and create a standard
% z-axis

if nargin==1
    Mxy=z;
    Nz=length(z);
    z=-Nz/2:Nz/(Nz-1):Nz/2;
    if not(mod(Nz,2))
        z=z-z(Nz/2 + 1);
    end
end

% 20/8/2015 If Mxy part of a multidimensional matrix, squeeze it.
if length(size(Mxy))>2
    Mxy=squeeze(Mxy);
end

if length(z)~=length(Mxy)
    error('PlotMxy: Axis and Mxy not of equal length');
end
% plotyy(z,abs(Mxy),z,unwrap(angle(Mxy)));
% Without unwrap..
[ax,h1,h2]=plotyy(z,abs(Mxy),z,angle(Mxy));
% set(h1,'linewidth',2)% to change the first line
% set(h2,'linewidth',2) % to change the second line
hold on;
plot(z,ones(Nz,1)*0.98,'-.k');
plot(z,ones(Nz,1)*0.99,'-.k');
plot(z,ones(Nz,1)*1,'-.k');
% Mag only
% plot(z,abs(Mxy),'b','linewidth',2)
% hold on;
% plot(z,angle(Mxy),'color',[0 0.5 0]);

hold off;
end