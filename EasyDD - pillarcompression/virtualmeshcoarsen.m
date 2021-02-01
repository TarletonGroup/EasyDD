%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether virtual nodes can be eliminated
% Francesco Ferroni
% Defects Group / Materials for Fusion and Fission Power Group
% Department of Materials, University of Oxford
% francesco.ferroni@materials.ox.ac.uk
% February 2014
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen(rn,links,connectivity,linksinconnect,fseg,lmin,lmax,MU,NU,a,Ec)
rnnew=rn;
%[lrn, lrn2]=size(rn);
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
linksnew=links;
fsegnew=fseg;
connectivitynew=connectivity;
linksinconnectnew=linksinconnect;
areamin2=areamin*areamin;
delta=1e-16;
i=1;
while i<=length(rnnew(:,1)) %the loop looks each nodes M
    %fprintf('%i \n',i);
    if (connectivitynew(i,1)==2) && rnnew(i,end)==67 %virtual nodes
        % the discretization node is normal so set up the conditions to check for whether is should be coarsened away
        % looking for the environnment of 'node i': 2 links and 2 other nodes M. 
        link1=connectivitynew(i,2);
        link2=connectivitynew(i,4);
        posi_link1=connectivitynew(i,3);
        posi_link2=connectivitynew(i,5);
        posnoti_link1=3-posi_link1;
        posnoti_link2=3-posi_link2;
        link1_nodenoti=linksnew(link1,posnoti_link1);
        link2_nodenoti=linksnew(link2,posnoti_link2);
        if rnnew(link1_nodenoti,end)~=67 || rnnew(link2_nodenoti,end)~=67
        % if neighbour nodes are not virtual nodes then continue M
        % if at least one neighbour node is a virtual node : see after the loop
            i=i+1;
            continue;
        end
        vec1=rnnew(link1_nodenoti,1:3)-rnnew(i,1:3);
        vec2=rnnew(link2_nodenoti,1:3)-rnnew(i,1:3);
        vec3=vec2-vec1;
        r1=sqrt(vec1*vec1');
        r2=sqrt(vec2*vec2');
        r3=sqrt(vec3*vec3');
        if r3<lmax %if the coarsening would result in a link length larger than lmax don't coarsen
            s=0.5*(r1+r2+r3);
            area2=(s*(s-r1)*(s-r2)*(s-r3));
            dvec1dt=rnnew(link1_nodenoti,4:6)-rnnew(i,4:6);
            dvec2dt=rnnew(link2_nodenoti,4:6)-rnnew(i,4:6);
            dvec3dt=dvec2dt-dvec1dt;
            dr1dt=(vec1*dvec1dt')/(r1+delta);
            dr2dt=(vec2*dvec2dt')/(r2+delta);
            dr3dt=(vec3*dvec3dt')/(r3+delta);
            dsdt=0.5*(dr1dt+dr2dt+dr3dt);
            darea2dt=dsdt*(s-r1)*(s-r2)*(s-r3);
            darea2dt=darea2dt+s*(dsdt-dr1dt)*(s-r2)*(s-r3);
            darea2dt=darea2dt+s*(s-r1)*(dsdt-dr2dt)*(s-r3);
            darea2dt=darea2dt+s*(s-r1)*(s-r2)*(dsdt-dr3dt);
            if ((area2<areamin2)&&(darea2dt<0.0d0))||((r1<lmin)||(r2<lmin))%remove single node critierion
                %the area is less than minimum and shrinking or one of the arms is less than lmin and shrinking
                [rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,~]=mergenodes(rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,link2_nodenoti,i,MU,NU,a,Ec);   
            else %((area2>=areamin2)|(dareadt>=0.0d0)) & ((r1>=lmin)|(dr1dt>=lmin)) & ((r2>=lmin)|(dr2dt>=lmin))
                i=i+1;
            end 
        else % r3>=lmax
            i=i+1;
        end 
    else % connectivitynew(i,1)>2
        i=i+1;
    end 
end %while loop