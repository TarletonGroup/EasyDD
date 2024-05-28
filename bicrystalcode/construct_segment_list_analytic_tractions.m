function segments=construct_segment_list_analytic_tractions(rn, links)
    %%===================================================================%%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 15/03/2018.
    %---------------------------------------------------------------------%
    % Function constructs a segment list for analytic tractions. Works the
    % same way as constructsegmentslist(rn,links) but does not take virtual
    % segments into account.

[LINKMAX,~]=size(links);

segments=zeros(LINKMAX,14);
nseg=0;
for i=1:LINKMAX
    n0=links(i,1);
    n1=links(i,2);
    if(     (        n0     ~= 0  &&    n1     ~= 0  ))
          if (rn(n0, end) == 67 || rn(n1, end) == 67 )
              continue
          end
        nseg=nseg+1;
        segments(nseg,:)=[links(i,1:5),rn(n0,1:3),rn(n1,1:3),links(i,6:8)];
    end
end
segments=segments(1:nseg,:);
