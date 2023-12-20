% Найдём максимальный фазовый сдвиг для каждой h, d
for h_ind = 1:length(h)
    for d_ind = 1:length(d)
        numRays = length(t.(d_ind){h_ind,1}{1,1});
        for cnt = 1:numRays
            rayPhaseShifts(cnt,:) = t.(d_ind){h_ind,1}{1, 1}(1,cnt).PhaseShift;
        end
        maxPhaseShift(h_ind, d_ind) = max(rayPhaseShifts);
    end
end