function [rn, vn, links, connectivity, linksinconnect, fseg] = updateMatricesBackward(rnnew, ...
        linksnew, connectivitynew, linksinconnectnew, fsegnew)
    rn = [rnnew(:, 1:3) rnnew(:, 7)];
    vn = rnnew(:, 4:6);
    links = linksnew;
    connectivity = connectivitynew;
    linksinconnect = linksinconnectnew;
    fseg = fsegnew;
end
