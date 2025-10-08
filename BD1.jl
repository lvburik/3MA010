"""Modeling exercise 1"""



function random_walk(steps)
    """
    performs a random walk of given number of steps
    returns the final position as a tuple (x,y)
    """
    x, y = 0, 0

    for i in 1:steps
        direction = rand(1:4) #1=up, 2=down, 3=left, 4=right

        if direction == 1
            y += 1
        elseif direction == 2
            y -= 1
        elseif direction == 3
            x -= 1
        else
            x += 1
        end
    end

    return [x, y]
end

function mean_squared_displacement(walks)
    """
    calculates the mean squared displacement from a list of walks
    each walk is a tuple (x,y)
    """
    n = length(walks)
    
    total_squared_displacement = sum(reduce(hcat, walks).^2)
    
    return total_squared_displacement / n
end


displacement = [mean_squared_displacement([random_walk(10^ii) for _ in 1:1000]) for ii in 1:5]