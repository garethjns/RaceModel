% Tests the race model inequality. Code originally from:
% http://link.springer.com/article/10.3758/BF03193160#
% Behavior Research Methods
% May 2007, Volume 39, Issue 2, pp 291-302
% Testing the race model inequality: An algorithm and computer programs
% Rolf Ulrich, Jeff Miller, Hannes Schröter
%
% Modified to remove artefacts from copying out of pdf file, layout
% clarified, verified to work with Matlab 2015a and a couple of other minor
% edits.
% http://matlaboratory.blogspot.co.uk/2015/03/race-model-inequality-script-r-ulrich-j.html
%
% RT inputs should be:
% - in ms
% - rows (columns work, but will error if n differs between conditions)
% - Rounded (or will error trying to access non-integer indices)
% P can be either colunm or rows, and from >0 P <=1
%
% X,Y,Z are arrays with RTs for conditions Cx, Cy, Cz, respectively.
% P is an array which contains the probabilities for computing
% percentiles.
% If Plot==true, a plot of the result is generated.

function [Xp, Yp, Zp, Bp] = RaceModel(X, Y, Z, P, Plot)

% Step 1: Determine Gx, Gy, and Gz

% Check for ties
[Ux, Rx, Cx] = ties(X);
[Uy, Ry, Cy] = ties(Y);
[Uz, Rz, Cz] = ties(Z);

% Get maximum t value
tmax = max(ceil(max([X, Y, Z])));
% T=1:1:tmax; T unused, removed

% Get function values of G
Gx = CDF(Ux, Rx, Cx, tmax);
Gy = CDF(Uy, Ry, Cy, tmax);
Gz = CDF(Uz, Rz, Cz, tmax);

% Step 2: Compute B = Gx plus Gy
t = 1:tmax;
B = NaN(1, numel(t)); % Preallocation added
for t = t
    B(t) = Gx(t) + Gy(t);
end

stop = 0;
% Check whether requested percentiles can be computed
if ~check(Ux(1), P(1), Gx)
    disp('Not enough X values to compute requested percentiles')
    stop = 1;
end
if ~check(Uy(1), P(1), Gy)
    disp('Not enough Y values to compute requested percentiles')
    stop = 1;
end
if ~check(Uz(1), P(1), Gz)
    disp('Not enough Z values to compute requested percentiles')
end

if stop
    Xp = NaN;
    Yp = NaN;
    Zp = NaN;
    Bp = NaN;
    return
end

% Step 3: Determine percentiles

Xp = GetPercentile(P, Gx, tmax);
Yp = GetPercentile(P, Gy, tmax);
Zp = GetPercentile(P, Gz, tmax);
Bp = GetPercentile(P, B, tmax);

% Generate a plot if requested
if Plot
    plot(Xp, P, 'o-', Yp, P, 'o-', Zp, P, 'o-', Bp, P, 'o-')
    axis([min([Ux Uy Uz])-10, tmax+10, -0.03, 1.03])
    grid on
    title('Test of the Race Model Inequality', 'FontSize', 16)
    xlabel('Time t (ms)','FontSize', 14)
    ylabel('Probability','FontSize', 14)
    legend('G_x(t)', 'G_y(t)', 'G_z(t)', 'G_x(t)+G_y(t)')
end

end


function OK = check(U1, P1, G)

OK = true;
for t = (U1-2):(U1+2);
    if (G(t)>P1) && (G(t-1)==0)
        OK = false;
        return
    end
end

end


function Tp = GetPercentile(P, G, tmax)

% Determine minimum of |G(Tp(i))-P(i)|
np = length(P);
for i = 1:np;
    cc = 100;
    for t = 1:tmax
        if abs(G(t)-P(i)) < cc
            c = t;
            cc = abs(G(t)-P(i));
        end
    end
    if P(i) > G(c)
        Tp(i) = c+(P(i)-G(c))/(G(c+1)-G(c)); %#ok<AGROW>
    else
        Tp(i) = c+(P(i)-G(c))/(G(c)-G(c-1)); %#ok<AGROW>
    end
end

end


function [U, R, C]=ties(W)

% Count number k of unique values and store these values in U.
W = sort(W);
n = length(W);
k = 1;
U(1) = W(1);
for i = 2:n
    if W(i) ~= W(i-1)
        k = k + 1;
        U(k) = W(i); %#ok<AGROW>
    end
end

% Determine number of replications R
R = zeros(1,k);
for i = 1:k
    for j = 1:n
        if U(i) == W(j)
            R(i) = R(i)+1;
        end
    end
end

% Determine the cumulative frequency
C = zeros(1,k);
C(1) = R(1);
for i = 2:k
    C(i) = C(i-1)+R(i);
end

end


function [G] = CDF(U, R, C, maximum)

G = zeros(1, maximum);
k = length(U);
n = C(k);

for i = 1:k
    U(i) = round(U(i));
end


for t = 1:U(1)
    G(t) = 0;
end

for t = U(1):U(2);
    G(t) = (R(1)/2 + (R(1) + R(2))/2 * (t-U(1))/(U(2)-U(1)))/n;
end

for i = 2:(k-1);
    for t = U(i):U(i+1);
        G(t) = (C(i-1) + R(i)/2 + (R(i) + R(i+1))/2 ...
            * (t-U(i))/(U(i+1)-U(i)))/n;
    end
end

for t = U(k):maximum;
    G(t) = 1;
end

end