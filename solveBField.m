function B0 = solveBField(X,I,L,h,dim)
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
B0 = zeros(1,1:size(X,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Biot Savart law
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dim == 1
for k = 1:PtsPEdge
    B0 +=  dl*I.*(X(3,:)-X_spire(3,PtsPEdge+k))./(norm(X-X_spire(:,PtsPEdge+k),'cols').^3);
    B0 += -dl*I.*(X(3,:)-X_spire(3,3*PtsPEdge+k))./(norm(X-X_spire(:,3*PtsPEdge+k),'cols').^3);
end
endif

if dim == 2
for k = 1:PtsPEdge
    B0 += -dl*I.*(X(3,:)-X_spire(3,k))./(norm(X-X_spire(:,k),'cols').^3);
    B0 +=  dl*I.*(X(3,:)-X_spire(3,2*PtsPEdge+k))./(norm(X-X_spire(:,2*PtsPEdge+k),'cols').^3);
end
endif

if dim == 3
for k = 1:PtsPEdge
    B0 +=  dl*I.*(X(2,:)-X_spire(2,k))./(norm(X-X_spire(:,k),'cols').^3);
    B0 += -dl*I.*(X(2,:)-X_spire(2,2*PtsPEdge+k))./(norm(X-X_spire(:,2*PtsPEdge+k),'cols').^3);
    B0 += -dl*I.*(X(1,:)-X_spire(1,PtsPEdge+k))./(norm(X-X_spire(:,PtsPEdge+k),'cols').^3);
    B0 +=  dl*I.*(X(1,:)-X_spire(1,3*PtsPEdge+k))./(norm(X-X_spire(:,3*PtsPEdge+k),'cols').^3);
end
endif
B0 *= (mu0/4/pi);