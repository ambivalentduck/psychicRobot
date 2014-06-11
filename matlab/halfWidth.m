function out = halfWidth(stats,varargin)

%% Code from multcompare.m
gmeans = stats.means(:);
n = stats.n(:);
s = stats.s;
ng = sum(n>0)-1;
df = stats.df;

if nargin>1
    alpha=varargin{1};
else
    alpha=.1;
end

ctype='tukey-kramer';
crit = getcrit(ctype, alpha, df, ng);
gcov = diag((s^2)./n);

%% Switch to one of its subroutines

ng = length(gmeans);
MM = zeros(ng,2);
MM(:,1) = gmeans(:);
MM(:,2) = sqrt(diag(gcov));
MM(isnan(MM(:,1)),2) = NaN;

M = nchoosek(1:ng, 2);      % all pairs of group numbers
M(1,5) = 0;                 % expand M to proper size
g1 = M(:,1);
g2 = M(:,2);
M(:,4) = gmeans(g1) - gmeans(g2);
i12 = sub2ind(size(gcov), g1, g2);
gvar = diag(gcov);
d12 = sqrt(gvar(g1) + gvar(g2) - 2 * gcov(i12));
delta = crit * d12;
M(:,3) = M(:,4) - delta;
M(:,5) = M(:,4) + delta;

d = zeros(ng, ng);
d(i12) = d12;
sum1 = sum(sum(d));
d = d + d';
sum2 = sum(d);
if (ng > 2)
    w = ((ng-1) * sum2 - sum1) ./ ((ng-1)*(ng-2));
else
    w = repmat(sum1, 2, 1) / 2;
end
out = crit * w(:);


function crit = getcrit(ctype, alpha, df, ng)
% Get the minimum of the specified critical values
crit = Inf;
[onetype,ctype] = strtok(ctype);

while(~isempty(onetype))
    if (length(onetype) == 1)
        switch onetype
            case 't', onetype = 'tukey-kramer';
            case 'd', onetype = 'dunn-sidak';
            case 'b', onetype = 'bonferroni';
            case 's', onetype = 'scheffe';
            case 'h', onetype = 'tukey-kramer';
            case 'l', onetype = 'lsd';
        end
    end
    if (isequal(onetype, 'hsd')), onetype = 'tukey-kramer'; end

    switch onetype
        case 'tukey-kramer' % or hsd
            crit1 = stdrinv(1-alpha, df, ng) / sqrt(2);

            % The T-K algorithm is inaccurate for small alpha, so compute
            % an upper bound for it and make sure it's in range.
            ub = getcrit('dunn-sidak', alpha, df, ng);
            if (crit1 > ub), crit1 = ub; end

        case 'dunn-sidak'
            kstar = nchoosek(ng, 2);
            alf = 1-(1-alpha).^(1/kstar);
            if (isinf(df))
                crit1 = norminv(1-alf/2);
            else
                crit1 = tinv(1-alf/2, df);
            end

        case 'bonferroni'
            kstar = nchoosek(ng, 2);
            if (isinf(df))
                crit1 = norminv(1 - alpha / (2*kstar));
            else
                crit1 = tinv(1 - alpha / (2*kstar), df);
            end

        case 'lsd'
            if (isinf(df))
                crit1 = norminv(1 - alpha / 2);
            else
                crit1 = tinv(1 - alpha / 2, df);
            end

        case 'scheffe'
            if (isinf(df))
                tmp = chi2inv(1-alpha, ng-1) / (ng-1);
            else
                tmp = finv(1-alpha, ng-1, df);
            end
            crit1 = sqrt((ng-1) * tmp);

        otherwise
            error('stats:multcompare:BadCType',...
                'Unknown CTYPE (critical value type) %s.', ctype);
    end

    if (~isnan(crit1)), crit = min(crit, crit1); end
    [onetype,ctype] = strtok(ctype);
end

function x = stdrinv(p, v, r)
%STDRINV Compute inverse c.d.f. for Studentized Range statistic
%   STDRINV(P,V,R) is the inverse cumulative distribution function for
%   the Studentized range statistic for R samples and V degrees of
%   freedom, evaluated at P.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.3.2.2 $  $Date: 2004/01/24 09:36:39 $

% Based on Fortran program from statlib, http://lib.stat.cmu.edu
% Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
% Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)

if (length(p)>1 | length(v)>1 | length(r)>1),
   error('stats:stdrinv:NotScalar','Scalar arguments required.'); % for now
end

[err,p,v,r] = distchck(3,p,v,r);
if (err > 0)
   error('stats:stdrinv:InputSizeMismatch',...
         'Non-scalar arguments must match in size.');
end

% Handle illegal or trivial values first.
x = zeros(size(p));
if (length(x) == 0), return; end
ok = (v>0) & (v==round(v)) & (r>1) & (r==round(r) & (p<1));
x(~ok) = NaN;
ok = ok & (p>0);
v = v(ok);
p = p(ok);
r = r(ok);
if (length(v) == 0), return; end
xx = zeros(size(v));

% Define constants
jmax = 20;
pcut = 0.00001;
tiny = 0.000001;
upper = (p > .99);
if (upper)
   uppertail = 'u';
   p0 = 1-p;
else
   uppertail = 'l';
   p0 = p;
end

% Obtain initial values
q1 = qtrng0(p, v, r);
p1 = stdrcdf(q1, v, r, uppertail);
xx = q1;
if (abs(p1-p0) >= pcut*p0)
   if (p1 > p0), p2 = max(.75*p0, p0-.75*(p1-p0)); end
   if (p1 < p0), p2 = p0 + (p0 - p1) .* (1 - p0) ./ (1 - p1) * 0.75; end
   if (upper)
      q2 = qtrng0(1-p2, v, r);
   else
      q2 = qtrng0(p2, v, r);
   end

   % Refine approximation
   for j=2:jmax
      p2 = stdrcdf(q2, v, r, uppertail);
      e1 = p1 - p0;
      e2 = p2 - p0;
      d = e2 - e1;
      xx = (q1 + q2) / 2;
      if (abs(d) > tiny*p0)
         xx = (e2 .* q1 - e1 .* q2) ./ d;
      end
      if (abs(e1) >= abs(e2))
         q1 = q2;
         p1 = p2;
      end
      if (abs(p1 - p0) < pcut*p0), break; end
	   q2 = xx;
   end
end
   
x(ok) = xx;

% ---------------------------------
function x = qtrng0(p, v, r)
% Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
% Calculates an initial quantile p for a studentized range
% distribution having v degrees of freedom and r samples
% for probability p, p.gt.0.80 .and. p.lt.0.995.

t=norminv(0.5 + 0.5 .* p);
if (v < 120), t = t + 0.25 * (t.^3 + t) ./ v; end
q = 0.8843 - 0.2368 .* t;
if (v < 120), q = q - (1.214./v) + (1.208.*t./v); end
x = t .* (q .* log(r-1) + 1.4142);

function xout = stdrcdf(q, v, r, upper)
%STDRCDF Compute c.d.f. for Studentized Range statistic
%   F = STDRCDF(Q,V,R) is the cumulative distribution function for the
%   Studentized range statistic for R samples and V degrees of
%   freedom, evaluated at Q.
%
%   G = STDRCDF(Q,V,R,'upper') is the upper tail probability,
%   G=1-F.  This version computes the upper tail probability
%   directly (not by subtracting it from 1), and is likely to be
%   more accurate if Q is large and therefore F is close to 1.

%   Copyright 1993-2007 The MathWorks, Inc. 
%   $Revision: 1.3.2.3 $  $Date: 2007/12/10 23:08:20 $

% Based on Fortran program from statlib, http://lib.stat.cmu.edu
% Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
% Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)
% Vectorized and simplified for MATLAB.  Added 'upper' option.

if (length(q)>1 | length(v)>1 | length(r)>1),
   error('stats:stdrcdf:NotScalar','Scalar arguments required.'); % for now
end
[err,q,v,r] = distchck(3,q,v,r);
if (err > 0)
   error('stats:stdrcdf:InputSizeMismatch',...
         'Non-scalar arguments must match in size.');
end
uppertail = 0;
if (nargin>3)
   if ~(  isequal(upper,'u') | isequal(upper,'upper') ...
        | isequal(upper,'l') | isequal(upper,'lower'))
      error('stats:stdrcdf:BadUpper',...
            'Fourth argument must be ''upper'' or ''lower''.');
   end
   uppertail = isequal(upper,'u') | isequal(upper,'upper');
end

% Accuracy can be increased by use of a finer grid.  Increase
% jmax, kmax and 1/step proportionally.
jmax = 15;          % controls maximum number of steps
kmax = 15;          % controls maximum number of steps
step = 0.45;        % node spacing
vmax = 120;         % max d.f. for integration over chi-square

% Handle illegal or trivial values first.
xout = zeros(size(q));
if (length(xout) == 0), return; end   
ok = (v>0) & (v==round(v)) & (r>1) & (r==round(r));
xout(~ok) = NaN;
ok = ok & (q > 0);
v = v(ok);
q = q(ok);
r = r(ok);
if (length(v) == 0), return; end
xx = zeros(size(v));

% Compute constants, locate midpoint, adjust steps.
g = step ./ (r .^ 0.2);
if (v > vmax)
   c = log(r .* g ./ sqrt(2*pi));
else
   h = step ./ sqrt(v);
   v2 = v * 0.5;
   c = log(sqrt(2/pi) .* r .* g .* h) - v2 + v2.*log(v2) - gammaln(v2);

   j=(-jmax:jmax)';
   hj = h * j;
   ehj = exp(hj);
   qw = q .* ehj;
   vw = v .* (hj + 0.5 * (1 - ehj .^2));
   C = ones(1,2*kmax+1);         % index to duplicate columns
   R = ones(1,2*jmax+1);         % index to duplicate rows
end

% Compute integral by summing the integrand over a
% two-dimensional grid centered approximately near its maximum.
gk = (0.5 * log(r)) + g * (-kmax:kmax);
w0 = c - 0.5 * gk .^ 2;
pz = normcdf(-gk);
if (~uppertail)
   % For regular cdf, use integrand as in AS 190.
   if (v > vmax)
      % don't integrate over chi-square
      x = normcdf(q - gk) - pz;
      xx = sum(exp(w0) .* (x .^ (r-1)));
   else
      % integrate over chi-square
      x = normcdf(qw(:,C) - gk(R,:)) - pz(R,:);
      xx = sum(sum(exp(w0(R,:) + vw(:,C)) .* (x .^ (r-1))));
   end
else
   % To compute the upper tail probability, we need an integrand that
   % contains the normal probability of a region consisting of a
   % hyper-quadrant minus a rectangular region at the origin of the
   % hyperquadrant.
   if (v > vmax)           % for large d.f., don't integrate over chi-square
      xhq   = (1 - pz) .^ (r-1);
      xrect = (normcdf(q - gk) - pz) .^ (r-1);
      xx = sum(exp(w0) .* (xhq - xrect));
   else                    % for typical cases, integrate over chi-square
      xhq   = (1 - pz) .^ (r-1);
      xrect = (normcdf(qw(:,C) - gk(R,:)) - pz(R,:)) .^ (r-1);
      xx = sum(sum(exp(w0(R,:) + vw(:,C)) .* (xhq(R,:) - xrect)));
   end
end

xout(ok) = xx;
