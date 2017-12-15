function [coords, burgers, n_dln] = extract_dislocation_nodes(dln_nodes,...
                                        dln_node_cnct)
    %% Generate dislocation line segments.
    dln   = constructsegmentlist(dln_nodes, dln_node_cnct)';
    n_dln = size(dln, 2);

    %% Extract dislocation line nodes and burgers vectors.
    coords = zeros(3 * n_dln, 2);
    coords(:, 1) = reshape(dln(6:8 , :), 3 * n_dln, 1);
    coords(:, 2) = reshape(dln(9:11, :), 3 * n_dln, 1);
    burgers      = reshape(dln(3:5,  :), 3 * n_dln, 1);
    
    %% Cleanup.
    clear dln;
end %function