function [hs,vrh,ko_go,ustart]=HSBounds(cij_local)
% This function finds the Hashin-Shtrikman bounds for material with elastic
% moduli cij in a 6x6 matrix.  It is the implementation of the variational 
% equations first derived by Hashin and Shtrikman in 1963 and later
% articulated in papers by Peselnick and Meister and Watt and Peselnick. 
% Usage:
%       [hs,vrh,ko_go,ustart]=HSBounds(cij)
% where:
%     cij is a 6x6 matrix of elastic constants of any symmetry  (GPa units)
%     hs is a 2x2 matrix with upper and lower "optimal" bounds for K and G
%     vrh gives, for reference, the Voigt-Reuss-Hill bounds
%     ko_go gives the properties of the reference material at the optimal point
%     ustart gives the starting point used for "go" in the calculations
%
% The following "nested" custom functions are included:
%    1. [xmax,hs]=edgeu('pos' or 'neg'); finds the positive definite
%        boundary at small ko and the negative definite boundary at large ko
%    2. [k,hsl]=edgek(y,'pos' or 'neg'); finds the positive/negative
%        boundary at fixed uo with increasing ko
%    3.  y=lowerbound(x) returns the hs value at the boundaries - used in
%         search for optimal values
%    4.  y=upperbound(x) returns the hs value at the boundaries - used in
%         search for optimal values
%    5.  [hs,value]=hscalc(ko,go,cij) does the Hashin-Shrtikman
%         calculations
%    6.  [K,G]=VRH(C) calculates Voigt-Reuss-Hill bounds
% The built-in to MATLAB search function fminbnd is used to find the optimal bounds
%  JMB 2013

% stop if matrix of moduli is not 6x6
[n,m]=size(cij_local);
if (n~=6 && m~=6)
    error('The input matrix is not 6x6')
end
% stop immediately if the input elastic matrix is not symmetric
 if not(isequal(cij_local,cij_local'))
     error('The input matrix is not symmetric')
 end
% Stop immediately if the input elastic moduli are not positive definite
   [~,tst]=eig(cij_local);
   if sum(sign(diag(tst)))~=6, 
       error('The input matrix of moduli is not positive definite')
   end


% the searches are confined to domains in ko and go as given here
kmin=1;
kmax=1000;
umin=1;
umax=1000;
convg=.00001;
sfty=.01;
sftyu=.01;
sftyl=.01;

[K,G]=VRH(cij_local);
vrh=[K G];

% Find lower bounds on HS moduli
% start by finding uo border for kmin
    [xmax,~]=edgeu('pos');

% Now find the the boundary point (go, ko) where Khs is maximized
    y=fminbnd(@lowerbound,xmax/2, xmax);
%got it, now calculate lower HS limits there and return
    [k,hsl]=edgek(y,'pos');
    ko_go=[k y];
    
% Find upper bounds on HS moduli
% start by finding uo border at kmax
    [xmin,~]=edgeu('neg');

% Now find the the boundary point (go) where Khs is minimized  
    y=fminbnd(@upperbound,xmin, 1.5*xmin);

%got it, now calculate lower HS limits there and return
    [k,hsu]=edgek(y,'neg');
    ko_go=[ko_go;k y];
    hs=[hsu;hsl];
    ustart=[xmax xmin];
    hs=fix(hs*100)/100;
    vrh=fix(vrh*100)/100;
    
    
function y=lowerbound(x)
 % this function is called by fminbnd to find the point on boundary that is
 % the maximum HS value
    [~,hsn]=edgek(x,'pos');
    y=-hsn(1);
end

function y=upperbound(x)
% this function is called by fminbnd to find the point on boundary that is
% the minimum HS value
    [~,hsn]=edgek(x,'neg');
    y=hsn(1);
end


function [un,hs]=edgeu(flg)
% this function finds the boundary between positive(negative) regimes along
% the uo axis
if strcmp(flg,'pos');
    u1=umin;
    u2=umax;
    [~,value]=hscalc(kmin,u1,cij_local);
    if (value==0 || value==-1),
        error('wrong sign at start'),
    end
    [~,value]=hscalc(kmin,u2,cij_local);
        if (value==1),
            error('wrong sign for first kmax'),
        end 
    du=abs(umax-umin)/2;
    uo=umin+du;  
    du=du/2;
    [~,vo]=hscalc(kmin,uo,cij_local);
while du>convg
        if vo==0
            un=uo-du;
        else
            un=uo+du;    
        end
        [hs,vn]=hscalc(kmin,un,cij_local);
        du=abs(un-uo)/2;
        uo=un;
        vo=vn;
end
un=un-sfty;  % push into positive definite area - try to avoid singularities at boundary.

elseif strcmp(flg,'neg')
    u1=umin;
    u2=umax;
    [~,value]=hscalc(kmax,u1,cij_local);
        if (abs(value)==1) ,error('wrong sign at start'),end
    [~,value]=hscalc(kmax,u2,cij_local);
        if value==0,error('wrong sign for kmax at large Go'),end 
    du=abs(umax-umin)/2;
    uo=umin+du;  
    du=du/2;
    [~,vo]=hscalc(kmax,uo,cij_local);

while du>convg
       if vo==0
            un=uo+du;
        else
            un=uo-du;    
        end
        [hs,vn]=hscalc(kmax,un,cij_local);
        du=abs(un-uo)/2;
        uo=un;
        vo=vn;
end
un=un+sfty;  % push into negative definite area - try to avoid singularities at boundary.
end
end

function [kn,hs]=edgek(u,flg)
% this function finds the boundary between positive(negative) regimes along
% the ko axis

if strcmp(flg,'pos');
    k1=kmin;
    k2=kmax;
    [~,value]=hscalc(k1,u,cij_local);
        if value==0,error('wrong sign at start'),end
    [~,value]=hscalc(k2,u,cij_local);
        if abs(value)==1,error('wrong sign for first kmax'),end 
    dk=abs(kmax-kmin)/2;
    ko=kmin+dk;  
    dk=dk/2;
    [~,vo]=hscalc(ko,u,cij_local);
    
while dk>convg
        if vo==0
            kn=ko-dk;
        else
            kn=ko+dk;    
        end
        [hs,vn]=hscalc(kn,u,cij_local);
        dk=abs(kn-ko)/2;
        ko=kn;
        vo=vn;
end
kn=kn-sftyl;  % insure we are fully in positive definite area

elseif strcmp(flg,'neg')
    k1=kmin;
    k2=kmax;
    [~,value]=hscalc(k1,u,cij_local);
        if abs(value)==1,error('wrong sign at start'),end
    [~,value]=hscalc(k2,u,cij_local);
        if value==0,error('wrong sign for first kmax'),end 
    dk=abs(kmax-kmin)/2;
    ko=kmin+dk;  
    dk=dk/2;
    [~,vo]=hscalc(ko,u,cij_local);

while dk>convg
       if vo==0
            kn=ko+dk;
        else
            kn=ko-dk;    
        end
        [hs,vn]=hscalc(kn,u,cij_local);
        dk=abs(kn-ko)/2;
        ko=kn;
        vo=vn;
end
kn=kn+sftyu;  % insure we are fully in positive definite area
end
end

function [hs,value]=hscalc(ko,go,cij)
% Implementation of the calculation of Hashin Shtrikman moduli 
%  for a set of elastic constants of any symmetry.  
% Usage:
%  [hs,value]=hscalc(ko,go,cij)
% where:
%      ko and go are the bulk and shear modulus of the reference isotropic material, 
%      cij is 6x6 matrix of elastic moduli 
%      hs are [khs ghs] and 
%      for the tensor R, the variable "value" is +1 for positive definite R and -1 for negative definite 
% The equations are from Peselnick and Meister 1965 through Watt and Peselnick
% 1980 which are based on Hashin and Shtrikman 1963
%       JMB 8/2013

% the isotropic compliances elements  (eq 10-11 Watt Peselnick 1980)
    alpha=-3/(3*ko+4*go);
    beta=-3*(ko+2*go)/(5*go*(3*ko+4*go));
    gamma=(alpha-3*beta)/9;
% need to track the identity operator for 4th rank elasticity tensors.  Has .5 or 2 in
% the last three diagonal locations.
    I=eye(6,6);
    Iinv=I;
    Iinv(4:6,4:6)=2*I(4:6,4:6);
    I(4:6,4:6)=.5*I(4:6,4:6);
   
% set up isotropic elastic matrix 
    co=2*go*I;
    co(1:3,1:3)=co(1:3,1:3)+(ko-2/3*go)*ones(3,3);
% take difference between anisotropic and isotropic matrices and take
% inverse to get the residual compliance matrix (eq 4 & 8 Watt Peselnick 1980) 
    R=cij-co;
    H=inv(R);
% Subtract the isotropic compliances from the residual compliance matrix
% and invert back to moduli space.  The matrix B is the moduli matrix that
% when multiplied by the strain difference between the reference isotropic 
% response and the isotropic response of the actual material gives the
% average stres <pij> (eq 15 - 16 Watt Peselnick 1980) note identity
% matrix is Iinv
    A=H-beta*Iinv;
    A(1:3,1:3)=A(1:3,1:3)-gamma*ones(3,3);
    B=inv(A);
% now perform the averaging of the B matrix to get two numbers (B1 and B2)
% that are related to the two isotropic moduli - perturbations of the
% reference values (eq 21 and 22 Watt Peselnick 1980)
    sB1=sum(sum(B(1:3,1:3)));
    dB=diag(B);
    sB2=sum(dB(1:3))+2*sum(dB(4:6));
    B1=(2*sB1-sB2)/15; 
    B2=(3*sB2-sB1)/30;  % should be 30
% The Hashin-Shtrikman moduli are then calculated as changes from the
% reference isotropic body. (eq 25 & 27 Watt and Peselnick 1980)
    khs=ko+(3*B1+2*B2)/(3+alpha*(3*B1+2*B2));
    ghs=go+B2/(1+2*beta*B2);
    hs=[khs ghs];
% The moduli are valid in two limits - minimizing (or maximizing) the
% anisotropic difference elastic energy. Think of it as the maximum positive deviations 
% from the reference state or the maximum negative deviations from the reference state.
% These are either in the positive definite or negative definite regime of the matrix R. 
% Here is a test of the properties of R.  A +1 is returned if positive definite, a -1 is
% retunred if negative definite.  0 returned otherwise.  
    [~,D]=eig(R);
    s=sum(sign(diag(D)));
    value=0;
    if s==6, 
         value=1;
    elseif s==-6, 
        value=-1;
    end
end

function [K,G]=VRH(C)
%determines Voigt-Ruess-Hill elastic moduli
%given the matrix of elastic constants
S=inv(C);
Kv=((C(1,1)+C(2,2)+C(3,3))+(2*(C(1,2)+C(2,3)+C(3,1))))/9;
muv=((C(1,1)+C(2,2)+C(3,3))-(C(1,2)+C(2,3)+C(3,1))+3*(C(4,4)+C(5,5)+C(6,6)))/15;
Kr=1/((S(1,1)+S(2,2)+S(3,3))+2*(S(1,2)+S(2,3)+S(3,1)));
mur=15/(4*(S(1,1)+S(2,2)+S(3,3))-4*(S(1,2)+S(2,3)+S(3,1))+3*(S(4,4)+S(5,5)+S(6,6)));
Kh=(Kv+Kr)/2;
muh=(muv+mur)/2;
K=[Kv;Kr;Kh];
G=[muv;mur;muh];
end

end
