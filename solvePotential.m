function A0 = solvePotential(X,I,L,h,dim)
PtsPEdge = 30;
mu0 = 4 * pi * 1e-7;
if dim != 1 && dim != 2 && dim != 3
    fprintf("not compatible dimension, dim can be either 1,2 or 3\nunexpet behaviour ")
endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set spire points                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dl = L/PtsPEdge;
X_spire = zeros(3,PtsPEdge*4);
X_spire(1,1:PtsPEdge) = linspace(-L/2+dl/2,L/2-dl/2,PtsPEdge);
X_spire(1,PtsPEdge+1:2*PtsPEdge) = L/2;
X_spire(1,2*PtsPEdge+1:3*PtsPEdge) = linspace(L/2-dl/2,-L/2+dl/2,PtsPEdge);
X_spire(1,3*PtsPEdge+1:4*PtsPEdge) = -L/2;

X_spire(2,1:PtsPEdge) = -L/2;
X_spire(2,PtsPEdge+1:2*PtsPEdge) = linspace(-L/2+dl/2,L/2-dl/2,PtsPEdge);
X_spire(2,2*PtsPEdge+1:3*PtsPEdge) = L/2;
X_spire(2,3*PtsPEdge+1:4*PtsPEdge) = linspace(L/2-dl/2,-L/2+dl/2,PtsPEdge);

X_spire(3,:) = h;
A0 = zeros(1,1:size(X,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            mu0   / jx * dV(Y)     mu0   /  Ix *dl(Y)
% A0x(X)  = ---- * | ---------- =  ---- * | ----------
%           4*pi   /V  |X - Y|     4*pi   /L  |X - Y|
% reference: https://www.feynmanlectures.caltech.edu/II_14.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isy = PtsPEdge*(dim == 2);
isnz = !(dim == 3);
for k = 1:PtsPEdge
    A0 +=  dl*I./norm(X-X_spire(:,isy+k),'cols')            *isnz;
    A0 += -dl*I./norm(X-X_spire(:,isy+2*PtsPEdge+k),'cols') *isnz;
end
A0 *= (mu0/4/pi);