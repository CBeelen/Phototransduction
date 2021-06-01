% Resample vector of PDE
% input args: endtime and time increment of new vector, old time and PDE
% output args: new time and PDE vectors

function [new_time, new_PDE] = resamplePDE(endtime, increment, old_time, old_PDE)
    n_time = endtime/increment;
    new_time = linspace(0,endtime,n_time+1);
    new_PDE = zeros(1,n_time+1);
    k=1;
    for i=1:n_time+1
        if new_time(i) < old_time(1)
            new_PDE(i) = 0;
        elseif new_time(i) < old_time(k)
            new_PDE(i) = old_PDE(k-1);
        elseif new_time(i) >= old_time(length(old_time))
            new_PDE(i) = old_PDE(length(old_time));
        else
            while new_time(i) >= old_time(k)
                k = k+1;
            end
            new_PDE(i) = old_PDE(k-1);
        end
    end
end
