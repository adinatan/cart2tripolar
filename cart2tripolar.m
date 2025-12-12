function [ira, varargout] = cart2tripolar(im, qParams, opts)
% Cartesian to polar pixel transform - with compact representation
% The code starts with a Cartesian image of square pixel size dx*dy=1, and
% outputs a polar image of polar pixels that has a similar amount of
% information in polar coordinates dr*dtheta~=1.  The signal in each polar
% pixel is determined via its fractional overlap with the four surrounding
% Cartesian pixels. The result is a triangular polar representation,
% because the number of polar pixels per radius per angle increase with
% radius, i.e. the larger the radius the higher the angular resolution.
% The code was originally used in the context of analyzing images with
% symmetry around a quadrant so that functionality was kept.
% The code support NaN valued pixels for masking.
%
% Inputs:
%   im      : 2D image (double/float ok), may contain NaN mask
%   qParams : subset of 1:4 (default 1:4)
%   opts    : struct with fields:
%       .include_corners (true/false)  default true
%       .method          'integrate' (default) or 'sample'
%       .Nr              radial quad order  (default 3)
%       .Nt              angular quad order (default 3)
%       .interp          'linear' (default) or 'nearest'
%
% Outputs:
%   ira     : [L x (max(PPR)+1) x 4] triangular polar per quadrant
%   ira_tot : stitched triangular for qParams (optional)
%   spol    : squared polar (optional)
%   area    : area per polar pixel (optional 4th output)

%
% Example:
%------------
% load('testimg.mat');
% [ira ira_tot spol]=cart2tripolar(im);
% figure(2);
% for n=1:4
%     subplot(2,2,n); imagesc(ira(:,:,n));  title(['ira - quadrant # ' num2str(n) ]);
% end
%
% figure('Position',[50 50 900 250]);
% subplot(1,3,1); imagesc(im); axis square; title('Img');
% subplot(1,3,2); imagesc(ira_tot); axis square; title('ira__tot');
% subplot(1,3,3); imagesc([0 2*pi],1:size(spol,1),spol); axis square; title('spol');
%
% __________________________________
%   Adi Natan (natan@stanford.edu)
%   Ver 1.2 , Date: Dec 12th 2025


%% defaults
if nargin < 1 || isempty(im)
    try, load('tim.mat'); catch, im = peaks(127); end
end
if nargin < 2 || isempty(qParams)
    qParams = 1:4;
end
if nargin < 3, opts = struct; end

opts = setDefault(opts, 'include_corners', true);
opts = setDefault(opts, 'method', 'integrate');   % 'integrate' or 'sample'
opts = setDefault(opts, 'Nr', 3);                 % quadrature order in r
opts = setDefault(opts, 'Nt', 3);                 % quadrature order in theta
opts = setDefault(opts, 'interp', 'linear');      % 'linear' or 'nearest'

qParams = qParams(:).';
assert(all(ismember(qParams,1:4)), 'qParams must be subset of 1:4');
im = double(im);

%% geometry / padding
if opts.include_corners
    ims = size(im);
    if ims(1) > ims(2)
        im(:, end+1:ims(1)) = NaN;
    elseif ims(1) < ims(2)
        im(end+1:ims(2), :) = NaN;
    end
    im = padarray(im, ceil(size(im)*(sqrt(2)-1)/2), NaN);
end

nx = size(im,1); ny = size(im,2);
x0 = floor((nx+1)/2);
y0 = floor((ny+1)/2);
L  = min([x0, y0, nx-x0+1, ny-y0+1]);   % quadrant size

% radius centers in "pixel-center" units: rc = r-1
rc = 0:(L-1);

% improved angular bin rule: want arc length ~ 1 at mid-radius (rc+0.5)
% arc ≈ (rc+0.5)*dtheta => choose ntheta ≈ (pi/2)*(rc+0.5)
PPR = max(1, round( (pi/2) * (rc + 0.5) ));     % number of angular bins per quadrant
AngleInc = (0.5*pi) ./ PPR(:);
AngleInc(1) = 0;

%% build quadrants (each LxL), with (1,1)=origin pixel
Q = NaN(L,L,4);
Q(:,:,1) = im(x0:-1:x0-L+1, y0:y0+L-1);
Q(:,:,2) = im(x0:-1:x0-L+1, y0:-1:y0-L+1);
Q(:,:,3) = im(x0:x0+L-1   , y0:-1:y0-L+1);
Q(:,:,4) = im(x0:x0+L-1   , y0:y0+L-1);

%% allocate
maxCols = max(PPR) + 1;               % +1 because we keep qp=0:npr edges
ira  = NaN(L, maxCols, 4);
area = NaN(L, maxCols, 4);

% quadrature nodes/weights on [-1,1]
[xr, wr] = gaussLegendre(opts.Nr);
[xt, wt] = gaussLegendre(opts.Nt);

for a4n = qParams
    a4 = Q(:,:,a4n);

    % interpolant for integration / sampling
    F = griddedInterpolant({1:L, 1:L}, a4, opts.interp, 'none');

    ira(1,1,a4n)  = a4(1,1);
    area(1,1,a4n) = 1;  % arbitrary for origin; not used

    for r = 2:L
        ntheta = PPR(r);            % number of bins
        dtheta = AngleInc(r);

        % radial pixel footprint in r-units (centered at rc=r-1)
        r1 = (r-1) - 0.5;
        r2 = (r-1) + 0.5;
        if r1 < 0, r1 = 0; end

        for k = 1:(ntheta+1)
            % keep qp=0..ntheta (edges included, matches your original)
            % define theta "pixel" around that sample index:
            % for interior points, cell is centered; for endpoints, clamp.
            if k == 1
                t1 = 0;
                t2 = 0.5*dtheta;
            elseif k == ntheta+1
                t1 = 0.5*pi - 0.5*dtheta;
                t2 = 0.5*pi;
            else
                tc = (k-1)*dtheta;
                t1 = tc - 0.5*dtheta;
                t2 = tc + 0.5*dtheta;
            end

            if strcmpi(opts.method,'integrate')
                % --- area-consistent quadrature over [r1,r2]x[t1,t2] with Jacobian r ---
                % map quad nodes from [-1,1] -> [a,b]
                rr = 0.5*(r2-r1)*xr + 0.5*(r2+r1);    % Nr x 1
                tt = 0.5*(t2-t1)*xt + 0.5*(t2+t1);    % Nt x 1

                % tensor grid
                [RR, TT] = ndgrid(rr, tt);

                % quadrant coords (1-based)
                X = RR .* sin(TT) + 1;    % row coordinate
                Y = RR .* cos(TT) + 1;    % col coordinate

                V = F(X, Y);              % may contain NaN

                % weights include Jacobian r and interval scaling
                W = (0.5*(r2-r1)*wr) * (0.5*(t2-t1)*wt).';  % Nr x Nt (outer product)
                W = W .* RR;                                   % Jacobian

                mask = ~isnan(V);
                sw = sum(W(mask));
                if sw > 0
                    ira(r,k,a4n) = sum(W(mask).*V(mask)) / sw; % area-weighted average
                else
                    ira(r,k,a4n) = NaN;
                end

                % analytic area of polar cell (for reference / conservation checks)
                area(r,k,a4n) = 0.5*(r2^2 - r1^2) * (t2 - t1);

            else
                % --- fast sampling fallback (not as accurate as integrate) ---
                tc = (k-1)*dtheta;
                xp = (r-1)*sin(tc) + 1;
                yp = (r-1)*cos(tc) + 1;
                ira(r,k,a4n) = F(xp, yp);
                area(r,k,a4n) = 0.5*(r2^2 - r1^2) * min(dtheta, (t2-t1));
            end
        end
    end
end

%% optional outputs
if nargout > 1
    ira_tot = addquad_fixed(ira, qParams, PPR);
    varargout{1} = ira_tot;

    if nargout > 2
        spol = t2s_fixed(ira_tot, qParams, PPR);
        varargout{2} = spol;
    end
    if nargout > 3
        varargout{3} = area;
    end
end

end

%% helpers

function S = setDefault(S, name, val)
if ~isfield(S,name) || isempty(S.(name))
    S.(name) = val;
end
end

function [x,w] = gaussLegendre(n)
% Gauss-Legendre nodes/weights on [-1,1], n<=8 typical. Hardcode for stability.
switch n
    case 1
        x = 0; w = 2;
    case 2
        x = [-1 1]'/sqrt(3); w = [1 1]';
    case 3
        x = [-sqrt(3/5) 0 sqrt(3/5)]';
        w = [5/9 8/9 5/9]';
    case 4
        x = [-0.8611363116 -0.3399810436 0.3399810436 0.8611363116]';
        w = [0.3478548451 0.6521451549 0.6521451549 0.3478548451]';
    case 5
        x = [-0.9061798459 -0.5384693101 0 0.5384693101 0.9061798459]';
        w = [0.2369268851 0.4786286705 0.5688888889 0.4786286705 0.2369268851]';
    otherwise
        % fallback: use 3 as a reasonable default
        x = [-sqrt(3/5) 0 sqrt(3/5)]';
        w = [5/9 8/9 5/9]';
end
end

function ira_tot = addquad_fixed(ira, qParams, PPR)
Ms = size(ira,1);
nQ = numel(qParams);
ira_tot = NaN(Ms, nQ*(max(PPR)+1));

for r = 2:Ms
    npr = PPR(r) + 1;
    irav = NaN(npr, nQ);

    for k = 1:nQ
        a4n = qParams(k);
        row = ira(r,1:npr,a4n);
        if mod(a4n,2)==0
            row = row(end:-1:1);
        end
        irav(:,k) = row(:);
    end
    ira_tot(r, 1:npr*nQ) = reshape(irav, 1, []);
end

ira_tot(1,1) = mean(ira(1,1,qParams), 3, 'omitnan');
end

function squmat = t2s_fixed(ira_tot, qParams, PPR)
Ms = size(ira_tot,1);
nQ = numel(qParams);
W  = size(ira_tot,2);

squmat = NaN(Ms, W);
squmat(1,:) = ira_tot(1,1);

for r = 2:Ms
    npr = (PPR(r)+1) * nQ;
    row = ira_tot(r, 1:npr);

    if all(isnan(row)), continue; end

    xi  = linspace(1, npr, W);
    idx = max(1, min(npr, round(xi)));

    squmat(r,:) = row(idx);
    squmat(r, isnan(row(idx))) = NaN;
end
end