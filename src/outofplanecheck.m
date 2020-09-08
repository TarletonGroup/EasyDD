function [index] = outofplanecheck(rn, links)
    index = false(size(links, 1), 1);

    for i = 1:size(links)

        if rn(links(i, 1), end) == 67 || rn(links(i, 2), end) == 67
            continue
        end

        linvec = rn(links(i, 2), 1:3) - rn(links(i, 1), 1:3);
        linvec = linvec / norm(linvec);
        normvec = links(i, 6:8);
        dotprod = linvec * normvec';

        if abs(dotprod) > 0.01
            index(i) = true;
        end

    end

end
