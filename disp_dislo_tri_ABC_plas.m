function [u]=disp_dislo_tri_ABC_plas(A,B,C,P,b)

%written by F.Hofmann 9/7/09
%routine to compute the displacement field for a triangular dislocation
%loop ABC at point P

% modified by F.Hofmann 5/11/18

%inputs:
%A,B,C - three column vectors defining the nodes of the dislocation loop
%P - column vector with coordinates of the point at which the displacement
%is evaluated in dimension 1. Then in dimension 2 this is a list of poinsts at which to evaluate the field. 
%b - column vector with 3 burgers vector components

%define local R vectors...
A_big = A*ones(1,size(P,2));
B_big = B*ones(1,size(P,2));
C_big = C*ones(1,size(P,2));

RA = A_big - P;
RB = B_big - P;
RC = C_big - P;

RAn = sqrt(sum(RA.*RA,1));
RBn = sqrt(sum(RB.*RB,1));
RCn = sqrt(sum(RC.*RC,1));

%defining the lambda  unit vectors along R vectors...
lambA = RA./(ones(3,1)*RAn);
lambB = RB./(ones(3,1)*RBn);
lambC = RC./(ones(3,1)*RCn);
lambA(isnan(lambA))=0;
lambB(isnan(lambB))=0;
lambC(isnan(lambC))=0;

%% defining dislocation segment tangent vectors t...
AB=B-A;
% BC=C-B;
CA=A-C;
tAB = safenorm(AB);
% tBC = safenorm(BC);
tCA = safenorm(CA);

%% compute fs
% fAB = cross(b,tAB)*log((RBn/RAn)*((1+dot(lambB,tAB))./(1+dot(lambA,tAB))));
% fBC = cross(b,tBC)*log((RCn/RBn)*((1+dot(lambC,tBC))./(1+dot(lambB,tBC))));
% fCA = cross(b,tCA)*log((RAn/RCn)*((1+dot(lambA,tCA))./(1+dot(lambC,tCA))));

% dot_lamA_tAB = sum(lambA.* (tAB*ones(1,size(P,2))),1);
% dot_lamB_tAB = sum(lambB.* (tAB*ones(1,size(P,2))),1);
% 
% dot_lamB_tBC = sum(lambB.* (tBC*ones(1,size(P,2))),1);
% dot_lamC_tBC = sum(lambC.* (tBC*ones(1,size(P,2))),1);
% 
% dot_lamC_tCA = sum(lambC.* (tCA*ones(1,size(P,2))),1);
% dot_lamA_tCA = sum(lambA.* (tCA*ones(1,size(P,2))),1);

% 
% fAB = (cross(b,tAB)*ones(1,size(P,2)))...
%             .*(ones(3,1)* log((RBn./RAn).*((1+dot_lamB_tAB)./ (1+dot_lamA_tAB)))); % 3 x N size
% fBC = (cross(b,tBC)*ones(1,size(P,2)))...
%             .*(ones(3,1)* log((RCn./RBn).*((1+dot_lamC_tBC)./ (1+dot_lamB_tBC))));
% fCA = (cross(b,tCA)*ones(1,size(P,2)))...
%             .*(ones(3,1)* log((RAn./RCn).*((1+dot_lamA_tCA)./ (1+dot_lamC_tCA))));


%% compute gs
%gAB = (dot(b,cross(lambA,lambB)).*(lambA+lambB))./(1+dot(lambA,lambB));
%gBC = (dot(b,cross(lambB,lambC)).*(lambB+lambC))./(1+dot(lambB,lambC));
%gCA = (dot(b,cross(lambC,lambA)).*(lambC+lambA))./(1+dot(lambC,lambA));

% dot_b_lA_lB = sum(((b*ones(1,size(P,2))).* cross(lambA,lambB,1)),1); % 1 x N size
% dot_b_lB_lC = sum(((b*ones(1,size(P,2))).* cross(lambB,lambC,1)),1);
% dot_b_lC_lA = sum(((b*ones(1,size(P,2))).* cross(lambC,lambA,1)),1);

dot_lA_lB = sum((lambA.*lambB),1); % 1xN size
dot_lB_lC = sum((lambB.*lambC),1);
dot_lC_lA = sum((lambC.*lambA),1);

% gAB = (ones(3,1)*(dot_b_lA_lB./(1+dot_lA_lB))).*(lambA+lambB); % 3 x N size
% gBC = (ones(3,1)*(dot_b_lB_lC./(1+dot_lB_lC))).*(lambB+lambC);
% gCA = (ones(3,1)*(dot_b_lC_lA./(1+dot_lC_lA))).*(lambC+lambA);

%% compute solid angle
%compute a,b,c and s
% aa = acos(dot(lambB, lambC));
% bb = acos(dot(lambC, lambA));
% cc = acos(dot(lambA, lambB));
% ss = (aa+bb+cc)/2;

aa = acos(dot_lB_lC); % 1 x N size
bb = acos(dot_lC_lA);
cc = acos(dot_lA_lB);
ss = (aa+bb+cc)/2;

%compute omega value without sign
%ome_v = real(4*atan(sqrt(tan(ss/2)*tan((ss-aa)/2)*tan((ss-bb)/2)*tan((ss-cc)/2))));
%ome = -sign(dot(lambA,cross(tCA,tAB)))*ome_v;

ome_v = real(4.*atan(...
                sqrt(tan(ss/2).*tan((ss-aa)./2).*tan((ss-bb)./2).*tan((ss-cc)./2))...
                    )); % 1 x N size
ome = - sign(sum((lambA.*(cross(tCA,tAB)*ones(1,size(P,2)))),1)).* ome_v; % 1 x N size


%% assembly of the displacement vector u
%u = -(b*ome)./(4*pi)
%                    -(1-2*nu)/(8*pi*(1-nu))*(fAB+fBC+fCA)
%                                                         +1/(8*pi*(1-nu))*(gAB+gBC+gCA);

u = -((b*ones(1,size(P,2))) .* (ones(3,1)*ome))/(4*pi);...
%     -((1-2*nu)/(8*pi*(1-nu))) * (fAB+fBC+fCA)...
%     +(1/(8*pi*(1-nu)))*(gAB+gBC+gCA);   % 3 x N size overall...


% %experimental bit, to get rid of discontinuity, when computing strain, add
% %burgers vector to displacement of the measurement point is below the
% %plane...
% if (sign(dot(lambA,cross(tAB,(-tCA))))<=0);
%     u = u + b;
% end

end
function unitR = safenorm(R)

normR = sqrt(R(1)*R(1) + R(2)*R(2)+R(3)*R(3));

if normR > eps
    unitR = R/normR;
else
    unitR = zeros(3,1);
end

end
