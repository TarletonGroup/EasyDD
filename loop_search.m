function loop_list = loop_search(rn,links,maxconnections)

[~,~,node_connectivity]=genconnectivity(rn,links,maxconnections); 

num_nodes = size(rn,1);

%disp('Searching for all loops in network');
maxconn = size(node_connectivity,2);

loop_counter=1;
loop_list = [];

for i=1:num_nodes
    vector = zeros(1,3);
    
    for j=1:(maxconn-1)
        vector(1) = node_connectivity(i,j);
        if vector(1) == 0 %trivial
            continue;
        end
        for k=(j+1):maxconn;
            vector(end) = node_connectivity(i,k);
            if vector(end) == 0 %trivial
                continue;
            elseif vector(end)<vector(1) %swap
                temp=vector(end);   
                vector(end)=vector(1);
                vector(1)=temp;
            end
            %Is B connected to C?
            if isempty(find(node_connectivity(vector(1),:)==vector(end), 1))
                %not connected.
                %disp('not connected')
                vector(2)=i;
                candidates=vector;
                while not(isempty(candidates))
                    for p=1:size(candidates,1)
                        %append new nodes onto candidates
                        counter=1;
                        for l=1:size(vector,1)
                                new_nodes = node_connectivity(candidates(p,end),:)';
                                new_nodes = new_nodes(all(new_nodes,2),:); %remove zeros
                                candidates = repmat(candidates,size(new_nodes,1),1);
                                new_nodes = repmat(new_nodes,size(candidates,1)/size(new_nodes,1),1);
                                candidates = horzcat(candidates,new_nodes);
                                counter=counter+1;
                        end
                    end
                    index=true(1,size(candidates,1));
                    for p=1:size(candidates,1)
                        if not(isempty(find(node_connectivity(candidates(p,2:(end-2),:),:)==candidates(p,end), 1)))
                            index(p) = false;
                            continue;
                        elseif not(isempty(find(candidates(p,1:(end-1))==candidates(p,end), 1)))
                            index(p) = false;
                            continue;
                        elseif isempty(find(node_connectivity(candidates(p,1),:)==candidates(p,end), 1))
                            %if it's not connected, don't discard!
                            continue;
                        else
                            %better way IO?
                            loop = candidates(p,:)';
                            ID = repmat(loop_counter,length(loop),1);
                            loop_list = vertcat(loop_list , horzcat(loop,ID));
                            %once outputed, discard
                            index(p) = false;
                            loop_counter=loop_counter+1;
                        end
                    end
                    candidates=candidates(index,:);
                end 
            else
                %better way IO?
                loop = vertcat(vector(1),i,vector(end));
                ID = repmat(loop_counter,3,1);
                loop_list = vertcat( loop_list , horzcat(loop,ID) );
                loop_counter=loop_counter+1;
            end
        end
    end
    
    node_connectivity(i,:)=zeros(1,maxconn);
    [X,Y]=find(node_connectivity==i);
    for j=1:size(X,1);
        node_connectivity(X(j),Y(j))=0;
    end
        
end