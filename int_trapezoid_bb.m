function [rn,vn,dt,fn,fseg]=int_trapezoid_bb(rn,dt,dt0,MU,NU,a,Ec,links,connectivity,...
    rmax,rntol,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d)

global unique_errmag_point unique_distmag_point

unique_errmag_point = 0;
unique_distmag_point = 0;

%Implicit numerical integrator using the Euler-trapezoid method adapted
%from [Cai & Bulatov, Algorithm 10.2, pg. 216]. The tiemstep size is
%controlled so that it can increase suitably quickly and not decrease
%excessively whilst remaining within acceptable tolerence limits
%Written by B.Bromage and D.Celis-Garza 05/11/2018

%Convert rn into a single column of coordinates and store the node flags
rnvec0=[rn(:,1);rn(:,2);rn(:,3)]; flag=rn(:,4);
    
%Calculate the current nodal velocities
[vnvec0,fn,fseg]=drndt(rnvec0,flag,MU,NU,a,Ec,links,connectivity,...
        mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
    
maxiter=10;       %Maximum number of times the timestep can be increased
counter=1;        %Counter variable for the while loop
dt_old=0;         %Variable for the maximumum acceptable timestep
maxchange=1.2;    %Maximum factor the timestep can be increased by
exponent=20;      %Variable that controls the degree that the calculated error affects the cange in timestep
dt_old_good = 0;  %Logical which flags whether an acceptable timestep has been calculated
convergent=0;     %Logical for the operation of the while loop
large_scaling = 0;%HY20190509
stat_nodes = logical(rn(:,end)==7 | rn(:,end)==67);
stat_nodes=[stat_nodes;stat_nodes;stat_nodes];
    
    while(~convergent)
        
        unique_errmag_point = 0;
        unique_distmag_point = 0;
    
        rnvec1=rnvec0+vnvec0*dt; %Euler forward method [Cai & Bulatov, eq. 10.43]
        
        if isempty(rnvec1)  %If there are no dislocations then use the maximum possible timestep
            convergent=1;
        end
        
        %Calculate the nodal velocities for the next timestep accordin to
        %Euler forward method
        [vnvec1,fn,fseg]=drndt(rnvec1,flag,MU,NU,a,Ec,links,connectivity,...
            mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
        
        distvec=rnvec1-rnvec0;
%         distmag=max(abs(distvec)); %Largest distance any node moves
        
        err=distvec-(vnvec1+vnvec0)/2*dt; %Euler-trapezoid method
%         errmag=max(abs(err));  %Largest difference (error) calculated by the Euler-trapezoid method
        
        %HY20190515:
        %Search for unusually fast accelerating/moving nodes
        errtol=err(stat_nodes==0);
        errtol=mean(abs(errtol))+std(abs(errtol))*6;
        disttol=distvec(stat_nodes==0);
        disttol=mean(abs(disttol))+std(abs(disttol))*6;
        errsing=find(abs(err)>errtol);
        distsing=find(abs(distvec)>disttol);
        err(errsing)=err(errsing)/20;
        distvec(errsing)=distvec(errsing)/20;
        err(distsing)=err(distsing)/20;
        distvec(distsing)=distvec(distsing)/20;
        
        if ~isempty(errsing)
            unique_errmag_point = 1;
        end
        
        if ~isempty(distsing)
            unique_distmag_point = 1;
        end
%         [val1, ind1] = sort(abs(err),'descend');
%         if val1(1)>mean(val1)*20
% %             disp('unique errmag point')
%             unique_errmag_point = 1;
% %             val(1)/val(4);;
%             err(ind1(1))=0;
%             distvec(ind1(1))=0;
%         end
%         [val2, ind2] = sort(abs(distvec),'descend');
%         if val2(1)>mean(val2)*20
% %             disp('unique distmag point')
%             unique_distmag_point = 1;
% %             val2(1)/val2(4);
%             err(ind2(1))=0;
%             distvec(ind2(1))=0;
%         end
        
        errmag=max(abs(err));
        distmag=max(abs(distvec));
        
%         %HY20190508: added by HY to highlight the nodes with maximum error
%         [val ind] = sort(abs(err),'descend');
%         val(1:3)
%         ind(1:3)
%         max_nodes(ind(1:3))
%         plotnodesHIGHLIGHT(rn,links,plim,vertices,max_nodes)
        
        if isempty(errmag) %If the Euler-trapzoid method yields no results use maximum time step and end loop
            dt=dt0;
            break
        end
        
        
        if(errmag<rntol) && (distmag<rmax)  %If error and max distance move are in acceptable limits
            dt_old = dt;                    %Store current timestep as maximum acceptable timestep
            factor=maxchange*(1/(1+(maxchange^exponent-1)*(errmag/rntol)))^(1/exponent);
            
            %HY20190508:
            large_scaling = 0;
%             if (errmag<rntol*1E-2)
%                 factor = 1E0;
%                 counter=counter-1;
%                 large_scaling = 1;
%             elseif (errmag<rntol*0.1)
%                 factor = 5;
%                 counter=counter-1;
%                 large_scaling = 1;
% %             elseif (errmag<rntol*0.5)
% %                 factor = 1.9;
% %                 counter=counter-1;
% %                 large_scaling = 1;
%             end
            
            dt=min(dt*factor,dt0);          %Increase timestep depending on the magnitude of the error
            dt_old_good = 1;                %Flag acceptable timestep calculated
            counter=counter+1;              %Proceed to next iteration
        else
            if dt_old_good == 1 && large_scaling ==0            %If current timestep is too large, use largest acceptable timestep
                dt = dt_old;
                counter=maxiter;
            else
                dt=dt/2;                    %If no acceptable timestep has been calculated, halve timestep and try again
            end
        end
        
%         dt
%         counter
%         pause(0.2)
        
        if counter>maxiter || dt==dt0       %End loop if maximum number of iterations is reached or curren timestep is maximum timestep
           convergent=1; 
        end
        
    end
    
%     counter
    
    %HY20190508: added by HY to highlight the nodes with maximum error
%         [val ind] = sort(abs(err),'descend');
%         num = 10;
%         val(1:num)
%         ind(1:num)
%         max_index=(ind(1:num));
%         for i = 1:num
%             if mod(max_index(i),3)==0
%                 max_nodes(i) = max_index(i)/3;
%             else
%                 max_nodes(i) = floor(max_index(i)/3)+1;
%             end
%         end
%         figure(100)
%         plotnodesHIGHLIGHT(rn,links,4.2017e+04,vertices,max_nodes)
    
    %HY20190516:
    if unique_errmag_point==1 && unique_distmag_point==1
       srn=size(rn,1);
       index1=rem(errsing,srn);
       index1(index1(:)==0)=srn;
       index2=rem(distsing,srn);
       index2(index2(:)==0)=srn;
       index0=[index1;index2];
       vnvec0(index0) = vnvec0(index0)/20;
       vnvec0(srn+index0) = vnvec0(srn+index0)/20;
       vnvec0(2*srn+index0) = vnvec0(2*srn+index0)/20;
       vnvec1(index0) = vnvec1(index0)/20;
       vnvec1(srn+index0) = vnvec1(srn+index0)/20;
       vnvec1(2*srn+index0) = vnvec1(2*srn+index0)/20;
       rnvec1=rnvec0+vnvec0*dt;
       
    elseif unique_errmag_point==1 
        srn=size(rn,1);
        index1=rem(errsing,srn);
%         index1=rem(ind1(1),srn);
        index1(index1(:)==0)=srn;
        vnvec0(index1) = vnvec0(index1)/20;
        vnvec0(srn+index1) = vnvec0(srn+index1)/20;
        vnvec0(2*srn+index1) = vnvec0(2*srn+index1)/20;
        vnvec1(index1) = vnvec1(index1)/20;
        vnvec1(srn+index1) = vnvec1(srn+index1)/20;
        vnvec1(2*srn+index1) = vnvec1(2*srn+index1)/20;
        rnvec1=rnvec0+vnvec0*dt;
        
    elseif unique_distmag_point==1 
        srn=size(rn,1);
        index2=rem(distsing,srn);
%         index2=rem(ind2(1),srn);
        index2(index2(:)==0)=srn;
        vnvec0(index2) = vnvec0(index2)/20;
        vnvec0(srn+index2) = vnvec0(srn+index2)/20;
        vnvec0(2*srn+index2) = vnvec0(2*srn+index2)/20;
        vnvec1(index2) = vnvec1(index2)/20;
        vnvec1(srn+index2) = vnvec1(srn+index2)/20;
        vnvec1(2*srn+index2) = vnvec1(2*srn+index2)/20;
        rnvec1=rnvec0+vnvec0*dt; 
    end

    %Rearrange rnvec and vnvec into rn and vn matrices
    rn=[reshape(rnvec1,length(rnvec1)/3,3),flag];
    vn=reshape(vnvec1,length(vnvec1)/3,3);

