function [ rf90,rf180,tb,bout] = rootflipAE(n,nb,tblin,bandsep,CS,M)
%   Detailed explanation goes here
% 31/10/2015 sas Function that returns a pair of root-flipped pulses for a
% spin-echo sequence. The returned pulses are scaled in radians and defined
% in normalized time given by n-time points for the refocusing pulse and
% 2*n-time points for the excitation pulse


% n  : number of time points
% nb : number of bands/slices
% tb : Time-bandwidth product
% bandsep : number of slices of separation between bands (measured from
% passband centre to passband centre)
% CS :% Boolean for conjugate Symmetry. 1 for AM-only pulses
% M  : Number of passband roots per passband. Set to 2 for tb4. Otherwise see documentation

% Outputs:
% rf90: Excitation pulse in radians
% rf180: Refocusing pulse in radians
% tb : time-bandwidth product of the minimum-phase Multiband pulse. This is
% typically 30% less than the time-bandwidth of the linear-phase pulse.
%
bandsep = bandsep*tblin;

d1 = 0.01;             % combined Mxy ripple, passband
d2 = 0.01;             % combined Mxy ripple, stopband
                        % 1 Gauss = 0.1mT, or 1e-4T
rftype = 'matched';     % 'matched' or '2xrefd' (twice-refocused)

osfact = 10; % oversampling factor

N = 2*(n-1); % length of linear-phase filter we will factor to get min-phase filter

% directly design a MB min-phase filter, whose roots we can flip
% sas - adjust the ripples depending on 'matched' or 'twice-refocused'
if strcmp(rftype,'matched')
  d1 = d1/4; % target beta passband ripple, considering full 90-180 pair
  d2 = (d2/sqrt(2))^0.25; % target beta stopband ripple, considering full 90-180 pair
else
  d1 = d1/8; % target beta passband ripple, considering twice-refocused
  d2 = d2.^(1/4); % target beta stopband ripple, considering twice-refocused
end

nn = (0:N/2*osfact)'/(N/2*osfact);  % 0 to 1 - profile indices
d = zeros(N/2*osfact+1,1);          % passband mask
s = zeros(N/2*osfact+1,1);          % stopband mask
wts = zeros(N/2*osfact+1,1);        % ripple taper weights

dinfmin = 1/2*dinf(2*d1,d2^2/2); % d-infinity for a min-phase pulse with beta ripples (d1,d2)
dinflin = dinf(d1,d2);      % d-infinity for a linear phase pulse with the same ripples


tb = tblin/dinflin*dinfmin; % scale TBW product so as to get the same transition 
                            % width as linear phase pulse with same ripples, 
                            % after scaling back to desired slice thickness. This 
                            % makes comparison to other MB excitations more 
                            % meaningful, since all will have same slice characteristics.
w = dinfmin/tb; % transition width

if rem(nb,2) % if Nb odd
    % start out the f, m and w vectors with the DC band
    f = [0 (1-w)*(tb/2) (1+w)*(tb/2)];%*di/dilp;
    d = nn <= f(2)/(n/2); % target pattern
    wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
else
    f = 0;
end   

% add non-DC bands to the profiles
for ii = 1:floor(nb/2)
    cent = (ii - (rem(nb,2) == 0)/2)*(bandsep)*dinfmin/dinflin;
%     Old version of code used:
%     cent = (ii - (rem(nb,2) == 0)/2)*(bandsep-2/osfact)*dinfmin/dinflin
    f = [f (cent-(1+w)*(tb/2)) (cent-(1-w)*(tb/2)) (cent+(1-w)*(tb/2)) (cent+(1+w)*(tb/2))];
    d = d | (nn >= f(end-2)/(n/2) & nn <= f(end-1)/(n/2));
    s = s | (nn >= f(end-4)/(n/2) & nn <= f(end-3)/(n/2));
    nnc = nn - (f(end-1)+f(end-2))/2/(n/2); % indices centered with passband for weight calcs
    wts = max(wts,1./abs(nnc).^2); % quadratically-decaying ripple weights
end
% append the last stopband
s = s | (nn >= f(end)/(n/2));
wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

% build system matrix for cvx design
A = 2*cos(2*pi*(0:N/2*osfact)'*(-N/2:0)/(N*osfact));A(:,end) = 1/2*A(:,end);

% mask everything to get rid of transition bands
% Note how Ad has fewer rows than A, because the transition bands are not
% included.
Ad = A(s | d,:);
dd = double(d(s | d));
ss = wts(s | d).*double(s(s | d));

% use cvx to do the constrained optimization
cvx_begin
  variable delta(1) 
  variable x(N/2+1)
  minimize( delta )
  subject to
    -delta*dd <= Ad*x - dd <= delta*dd + delta*d2^2/(2*d1)*ss
cvx_end

% stack two halves together to get full linear-phase filter
x = [x;x(end-1:-1:1)]';
blin=x;

% factor the linear phase filter to get a min-phase filter b
b = real(fmp2(x));
b = b(end:-1:1);

% sas- Note How TBW scaling is used in the selection of passband locations
bp = 2*((1:nb)-1/2-nb/2)*bandsep/n*pi*dinfmin/dinflin; % centers of passbands
% sas - find bp as coordinates
bpc=1.*exp(1i*bp);

flip=pi;
N = length(b); % number of time points in pulse
b = b./max(abs(fft(b))); % normalized beta coefficients
b = b*sin(flip/2 + atan(d1*2)/2); % scale to target flip angle
[~,a] = b2amp(b);

rfinit_max = max(abs(islr(a,b))); % get initial peak RF amplitude

r = roots(b); % calculate roots of beta polynomial

% sort the roots by phase angle
[~,idx] = sort(angle(r)); 
r = r(idx);

rn = length(r);

% sas 02/06/2016
% What's different in my code from Sharma's code is that idxPass
% is of length M * NB and finds the M closest roots (i.e. smallest
% difference in radians from passband centre) where-as Sharma's code
% includes all roots which are 3 slice-thicknesses in radians near-to the
% passband center. 

% idxPass is a vector containing the indices of roots made elligible for
% flipping. For each passband, find the radial difference between each
% other root. Then find the M roots with the smallest difference and append
% them to idxPass.

idxPass = []; 
for i = 1:nb
    radial_diff = abs ( bp(i) - angle(r));

    [srd,sorted_radial_diff] = sort(radial_diff);

    %     idxPass = [idxPass sorted_radial_diff(1:M)'];

    % Added find(..) term to deal with DC-roots. DC roots are not considered
    % for flipping. This line adds the first M root-indices for which the
    % radial difference is not zero.
	diff_idc = (1:M)+find(srd>0,1)-1;
	
	% 04/10/16 Added code to ignore roots which appear in
	% the radial range of the passband roots, but have a
	% magnitude larger than 1.03
	if any(abs(r(sorted_radial_diff(diff_idc)))>1.03)
		k = 1;
    else
        k = 0;
	end
	idxPass = [idxPass sorted_radial_diff(diff_idc+k)'];    
end
idxPass = sort(idxPass);

% Depending on M, define a matrix called lpc (linear-phase combination,
% although this is not true in a DSP-sense!). Then iteratively add random
% combinations of the available rows into pFlip.

if M == 2
    lpc = [0 1;1 0];
    if and(CS == 1,mod(nb,2))
        warning('Not possible to create a aligned-echo AM Root-flipped pulse with a DC-slice');
    end
elseif M==3
warning('Echoes will not align in this case, as not possible to have equal number of roots inside and outside the unit when there are an uneven amount of roots per passband')
    lpc = [1 0 1;0 1 0];
elseif M==4
 % sas   11/04/2017 - I have not experimented a lot with flipping four
 % roots  per passband, but invite others to try.
 
%     lpc = [1 0 0 1;0 1 1 0];
%     Second possibility:
%     lpc = [0 1 0 1;1 0 1 0];
%     Third possibility:
%     lpc = [0 0 1 1;1 1 0 0];
%     lpc = [0 0 0 0;0 0 0 0];
%     lpc = [0 1 1 0];
    lpc = [0 1 1 0;1 0 0 1;0 1 0 1;1 0 1 0;0 0 1 1;1 1 0 0];
    if CS == 1
        % For AM-only, constrain to these two cases.
        lpc = [0 1 1 0;1 0 0 1];
    end
else
    error('Rootflip error: No combinations for M>4')
end
[nps,~] = size(lpc);


%%
ntrials = 50;
minpeakrf = Inf;
for ii = 1:ntrials;
    flip = zeros(1,rn);
    pFlip = [];
    for i=1:nb
        pFlip = [pFlip lpc(randi(nps),:)];
    end

    nrPb = length(pFlip);

    if CS==0
        symtype='TS';
    elseif CS==1
        symtype='CS';
    elseif CS==2
        symtype='AF';
    end

    % Manipulate pFlip depending on the desired symmetry.

    if strcmp(symtype,'TS')==1 %|| strcmp(symtype,'CS')==1
        %         pFlip = [pFlip(1:nrPb/2) not(pFlip(1:nrPb/2))];
        % 02/07/16 pFlip(nrPb/2 +1) is the
        % root on the bottom half of the uc. Keep that as it is, and conjugate the
        % other passband configs.
        %         pFlip = [pFlip(1:nrPb/2) pFlip(nrPb/2 + 1) not(pFlip(1:nrPb/2 -1))];
        if mod(nrPb,2) % If number of passband roots are odd
            pFlip = [pFlip(1:ceil(nrPb/2)) (not(pFlip(1:floor(nrPb/2))))];
        else
        % 12/08/2016 sas
        % If there are odd number of passband root on each half of the
        % unit circle (i.e. nb is odd) then keep the DC passband
        % unchanged.
            if mod(nb,2)
                pFlip = [pFlip(1:nrPb/2 + floor(M/2)) fliplr(not(pFlip(1:nrPb/2 - floor(M/2))))];
            else
                % If even on each side, keep the first half of pFlip and
                % conjugate the second half.
                pFlip = [pFlip(1:nrPb/2) fliplr(not(pFlip(1:nrPb/2)))];
            end
        end
    elseif strcmp(symtype,'CS')==1 
            if mod(nb,2)
                pFlip = [pFlip(1:nrPb/2 + floor(M/2)) fliplr((pFlip(1:nrPb/2 - floor(M/2))))];
            else
                pFlip = [pFlip(1:nrPb/2) fliplr((pFlip(1:nrPb/2)))];
            end
    elseif strcmp(symtype,'AF')==1
        pFlip = pFlip;
    end

    % Flip the roots!
    flip(idxPass) = pFlip;
    
    if strcmp(symtype,'TS')==1
        % 26/10 Find the element related to final bp(end)+wp
        fnz = floor( (N/2 + 1) + (bp(1)-3*tb*pi/N)*N/2/pi );
        % Find all roots on the bottom-halve of the complex plane that are beyond
        % the final flipped- passband root. This is to try and aim each stop-band
        tmp = and(angle(r)<0 , angle(r)< angle((r(fnz))) );

        flip = or(flip,tmp');
        
    elseif or(strcmp(symtype,'CS')==1,strcmp(symtype,'AF')==1)
        % 26/10 Find the element related to final bp(end)+wp
        fnz = floor( (N/2 + 1) + (bp(1)-3*tb*pi/N)*N/2/pi );

        tmp = and(angle(r)<0 , angle(r)< angle((r(fnz))) );
        
        % 26/10/15 Sinusoidal stop-band root-flipping!
        xt = -3*pi:6*pi/(fnz-2):3*pi;
        tmp(1:fnz-1) = round(0.5+0.87*sin(xt));
        % Repeat the same on top half of the uc
        if strcmp(symtype,'CS')==1 
            % If conjugate symmetric, keep top and bottom the same pattern
            tmp(end:-1:N-fnz+1) = round(0.5+0.87*sin(xt));
        elseif strcmp(symtype,'AF')==1
            % If all-flip, negate the xt term in the sin function.
            tmp(end:-1:N-fnz+1) = round(0.5+0.87*sin(-xt));
        end
        flip = or(flip,tmp');
    end
    
    doflip = flip;

    rt = r;
    rt(doflip == 1) = conj(1./rt(doflip==1));


    % get root-flipped RF
    R = poly(leja_fast(rt)); % get polynomial coefficients back from flipped roots
    R = R/max(abs(freqz(R))); % normalized by max of filter response

    bt = R*sin(pi/2 + atan(d1*2)/2); % scale to target flip
    [~,at] = b2amp(bt);
    rft=islr(at,bt);
    
    if max(abs(rft)) < minpeakrf
        minpeakrf = max(abs(rft));
        bestflips = doflip;
        rfout = rft;
        bout=bt;
        
        % sort the roots by phase angle
        [~,idx] = sort(angle(rt)); 
        rt = rt(idx);
        rout = rt;
        fprintf('Iteration %d of %d. Peak RF: %0.2d rad. Peak init RF: %0.2d rad.\n', ...
                ii,ntrials,minpeakrf,rfinit_max);
    end
end
%%
% sas - call min peak refocusing pulse rf180. Assign outside the MC
% iterations.
if strcmp(symtype,'CS')
    rf180=imag(rfout);
else
    rf180=rfout;
end


N = length(rf180); % # time points

[~,b90di] = rf2ab_sas(rf180,(-N/2:1/2:N/2-1/2)',0); % get beta of 180 pulse
b90d = (b90di.^2)/sqrt(2); % target 90-deg beta profile
b90d= -conj(b90d);

bx=fftshift( fft(ifftshift(b90d))/length(b90d));

[~,ax] = b2amp(bx);
rf90 = -1i*conj(islr(ax,bx));
if CS==1
    rf90=real(rf90);
end

end

