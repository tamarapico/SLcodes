function spa_sum = sph2spa_tp1(a_lm,maxdeg,long,colat)
% Function to read in spherical coordinates and convert to spatial. 
%T.Pico 2015
%
%Read in vector with
% l = 0, m=0, l=1, m=0, m=1, l = 2, m = 0...l, until l = maxdegree. 
%colat should be numberx1, long 1xnumber
% to convert I want to solve 
%U(theta,phi) = sum from l=0 to l= maxdeg of Ul0*Yl0 (all at m=0)
% + sum over m=1 to m = l, and over all l's of 2*Re(UlmYlm)
%ultimately these become Ulm*Plm*e^i*m*phi
%%

%initialize vector of m's from 0 to maxdegree
mvec = 0:maxdeg;
%initialize vector of exponents
evec = exp(-1i*mvec'*long*pi/180);

%matrix to store sum of this (will be output-spatial coord)
%size of spatial grid
spa_sum = zeros(length(colat),length(long));

for l =0:maxdeg
    
Plm = legendre_me(l,cos(colat.*pi/180),'me'); %get legendre polynomials
%Plm is m x colat

% if l = 0, assign coefficient, legendre polynomials 
if l == 0;
 Ulm = a_lm(1);%get coefficient of spherical grid
 Ul0 = Ulm;
 %legendre_me gives colat x 1, for l = 0, in higher degree it gives m x
 %colat, so flip to have consistent sizes
 Pl0 = Plm';
 Plm = Plm';
else
% assign coefficients Ulm
%position in a_lm based on sum(k) = n*(n+1)/2 for l!=0
pos_alm = l*(l+1)/2;
Ulm = a_lm(pos_alm+1:pos_alm+l+1); %get coefficients
Ul0 = Ulm(1);%if m = 0, only take first coeff
Pl0 = Plm(1,:);%get legendre polynomials for m=0
end

%assign emphi for each degree from initialized evec
emphi = evec(1:l+1,:); % size m x long

%take spherical coefficients from m=1 to m=l
Ulmr = Ulm(2:end); % size 1x m

%compute Ul0*Pl0, *emphi just gives correct size (emphi is all 1 for m=0)
UY0 = Ul0*(Pl0'*emphi(1,:)); %size lat x long
%repeat spherical coeff from m =1 to m= l to get same size of emphi (so
%fill by long)
Ulmrep = repmat(Ulmr', 1, length(long),1); %size m x long
%sum Ul0*Pl0 with 2*Re(Plm*Ulm*emphi), this sums all m = 0's, and all m=1
%to m = l, over every single l
%this output is the spatial conversion of the spherical harmonic grid
%first multiply scalar emphi and Ulmrep, then dot product to sum over m's
%with Plm --> gives lat x long matrix
spa_sum = spa_sum+ UY0 + 2*real(Plm(2:end,:)'*(emphi(2:end,:).*Ulmrep));

end

end