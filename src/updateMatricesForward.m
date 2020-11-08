function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = updateMatricesForward(rnnew, ...
        vn, links, connectivity, linksinconnect, fseg)
    rnnew = [rnnew(:, 1:3) vn rnnew(:, 4)];
    linksnew = links;
    connectivitynew = connectivity;
    linksinconnectnew = linksinconnect;
    fsegnew = fseg;
end
