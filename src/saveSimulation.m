function saveSimulation(simName, curstep, saveFreq)

    if mod(curstep, saveFreq) == 0
        close all;
        save(sprintf('../output/%s_%d', simName, curstep));
    end

end
