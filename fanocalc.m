
function [fftally] = fanocalc(spkTim,b)

%numpoints=1000;
%b = logspace(-1,2.5,numpoints); %spacing of points in units of spike time

fftally = zeros(1,length(b));

%parpool(7)
for gg=1:length(b)
    values = spkTim(1):b(gg):spkTim(end);
    if numel(values)>1
        temprate = histcounts(spkTim,values);
    else
        temprate = nan;
    end
    fftally(gg) = var(temprate)./mean(temprate);
end

%delete(gcp('nocreate'))

end


