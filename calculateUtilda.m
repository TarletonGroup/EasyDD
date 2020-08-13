function utilda = calculateUtilda(rn, links, gn, NU, xnodes, dx, dy,...
    dz, mx, my, mz, utilda)
%=========================================================================%
% Wrapper for Utilda calculation by Bruce Bromage.
%=========================================================================%

[Ux, Uy, Uz] = Utilda_bb_vec(rn, links, gn, NU, xnodes, dx, dy, dz, mx,...
    my, mz);

utilda(3*gn - 2) = Ux;
utilda(3*gn - 1) = Uy;
utilda(3*gn    ) = Uz;
end