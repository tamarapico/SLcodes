% Spatial to spherical coord
%
%

function readeq = spa2sph_tp(func,maxdegree,long,colat)
N = maxdegree;
n = length(colat);
[x,w] = GaussQuad(n);

alm = []; %alm is matrix where l goes down rows and m columns
dellon = (long(2)-long(1))*pi/180;
% gmx = zeros(length(mvec),length(x));
 gmp = zeros(N+1,length(x));
 for m = 0:N
     exps = exp(-1i.*m.*long'*pi/180);
     gmp(m+1,:) = (func*exps);
     %gmx(m+1,:) = dellon*spatial4*exps; %matrix multiplication sums over all phi
 end
gmf = gmp*dellon;
for l=0:N
     P_lm = legendre_me(l,cos(colat'*pi/180),'me');
     alm =[alm (w*(P_lm.*gmf(1:l+1,:))')/(4*pi)]; %sum over all x's

end

readeq = alm;
% %readeq = flipud(alm);
% readeq = NaN(N+1,N+1); %get into same format as read, need to move up every column
% for l=0:N
%    readeq(1:N+1-l,l+1) = alm(l+1:N+1,l+1);
% end

end