function f_tilda = calculateNumericTractions(rn, links, gamma_dln, xnodes, a, MU, NU, ...
        f_tilda, x3x6, n_nodes, n_nodes_t, n_se, idxi, f_tilda_node, ...
        f_tilda_se, use_gpu, n_threads, para_scheme, para_tol)

    [segments, ~] = constructsegmentlist(rn, links, true);
    f_tilda = numeric_traction(gamma_dln, segments, xnodes, a, MU, NU, f_tilda);

end
