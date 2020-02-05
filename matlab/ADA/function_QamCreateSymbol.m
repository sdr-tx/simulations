function [symbolVector] = function_QamCreateSymbol( pwmSize, duty1, duty2, duty3 )
    symbolVector = zeros(1, 3*pwmSize);
    
    for i = 1 : 3*pwmSize
        if i <= pwmSize
            if i <= duty1
                symbolVector(i) = 1;
            end
        elseif i <= 2*pwmSize
            if i <= pwmSize + duty2
                symbolVector(i) = 1;
            end
        else
            if i <= 2*pwmSize + duty3
                symbolVector(i) = 1;
            end    
        end
    end   
end



