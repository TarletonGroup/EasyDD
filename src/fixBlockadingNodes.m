function rnnew = fixBlockadingNodes(rnnew,connectivitynew)

    %find blockading fixed nodes and allow code to deal with them
    for i = 1:size(rnnew, 1)

        if rnnew(i, end) ~= 7
            continue
        end

        if connectivitynew(i, 1) > 7
            rnnew(i, end) = 0;
        end

    end

end
