"""
Diffusion Limited Aggregation (DLA) model simulation in Julia.
This code simulates the DLA process on a 2D lattice with periodic boundary conditions.
A seed particle is placed at the center of the lattice, and particles are spawned at the edges,
which then perform a random walk until they stick to the cluster.
The simulation continues until a particle reaches the boundary of the lattice.
The final lattice configuration is displayed in the console.
"""

#import necessary libraries
using Statistics #for mean function


#parameters
lattice_dimensions = [40, 40]

function DLA_Simulation(lattice_dimensions)
    """
    runs a DLA simulation on a lattice of given dimensions
    """

    # initialize lattice
    lattice = zeros(lattice_dimensions[1], lattice_dimensions[2])

    #put seed at center
    center = round.(Int, lattice_dimensions/2)
    lattice[center[1],center[2]] = 1;

    ## sim
    while check_lattice(lattice, lattice_dimensions) == false
        
        particle = spawn_particle(lattice_dimensions)
        while particle.freeze == false
            move_particle(particle, lattice_dimensions)
            check_particle(particle, lattice, lattice_dimensions)
        end

        #freeze particle in lattice
        lattice[particle.x, particle.y] = 1
    end

    return lattice
end

mutable struct Particle
    x :: Int
    y :: Int
    freeze :: Bool
end

function spawn_particle(lattice_dimensions)
    """spawn particle at boundary of lattice"""
    edge = rand(1:4) #pick which boundary

    if edge%2 == 0
        index = rand(1:lattice_dimensions[1])
    else
        index = rand(1:lattice_dimensions[2])
    end

    if edge == 1
        particle = Particle(1, index, false)
    elseif edge == 2
        particle = Particle(index, 1, false)
    elseif edge == 3
        particle = Particle(lattice_dimensions[1], index, false)
    else
        particle = Particle(index, lattice_dimensions[2], false)
    end
    
    return particle
end

function move_particle(particle, lattice_dimensions)
    """Move particle with periodic boundary conditions"""
    direction = rand(1:4)
    if direction == 1
        particle.x += 1
    elseif direction == 2
        particle.y += 1
    elseif direction == 3
        particle.x -= 1
    else
        particle.y -= 1
    end
    # Apply periodic boundary conditions
    particle.x = ((particle.x - 1) % lattice_dimensions[1] + lattice_dimensions[1]) % lattice_dimensions[1] + 1
    particle.y = ((particle.y - 1) % lattice_dimensions[2] + lattice_dimensions[2]) % lattice_dimensions[2] + 1
end


function check_particle(particle, lattice, lattice_dimensions)
    """check if particle should be frozen"""

    #check if particle at edge
    neighbors = [
        (particle.x + 1, particle.y),
        (particle.x - 1, particle.y),
        (particle.x, particle.y + 1),
        (particle.x, particle.y - 1)
    ]

    for (nx, ny) in neighbors
        # Apply periodic boundary conditions
        nx = ((nx - 1) % lattice_dimensions[1] + lattice_dimensions[1]) % lattice_dimensions[1] + 1
        ny = ((ny - 1) % lattice_dimensions[2] + lattice_dimensions[2]) % lattice_dimensions[2] + 1
        if lattice[nx, ny] == 1
            particle.freeze = true
            break
        end
    end

end



function check_lattice(lattice, lattice_dimensions)
    """
    check if lattice has frozen paricle at boundary
        """
    is_finished = false

    if 1 in lattice[1, :]
        is_finished = true
    elseif 1 in lattice[:, 1]
        is_finished = true
    elseif 1 in lattice[lattice_dimensions[1], :]
        is_finished = true
    elseif 1 in lattice[:, lattice_dimensions[2]]
        is_finished = true
    end

    return is_finished
end

function display_lattice(lattice)
    # Display the lattice in the console
    for i in 1:size(lattice, 1)
        for j in 1:size(lattice, 2)
            print(lattice[i, j] == 1 ? "██" : "  ")
        end
        println()
    end
end

function find_cluster_mass(lattice)
    """find mass of cluster"""
    return sum(lattice)
end

function find_cluster_length_scale(lattice)
    """"find length scale of cluster using radius of gyration"""
    coords = findall(x -> x == 1, lattice)
    n = length(coords)

    x_coords = [coord[1] for coord in coords]
    y_coords = [coord[2] for coord in coords]
    x_center = mean(x_coords)
    y_center = mean(y_coords)
    rg_squared = sum((x - x_center)^2 + (y - y_center)^2 for (x, y) in zip(x_coords, y_coords)) / n
    return sqrt(rg_squared)
end

function find__fractal_dimension(masses, length_scales)
    """find fractal dimension using linear regression"""
    log_masses = log.(masses)
    log_length_scales = log.(length_scales)
    
    # Linear regression: logM = D * logR + c
    X = hcat(ones(length(log_length_scales)), log_length_scales)
    coeffs = X \ log_masses
    D = coeffs[2]  # Slope
    return D
end

function run_DLA_and_analyze(lattice_dimensions, number_of_simulations=1, display_each_lattice=false)
    """run DLA simulation and find mean cluster properties"""

    
    fractal_masses = Vector{Int}(undef, number_of_simulations)
    fractal_length_scales = Vector{Float64}(undef, number_of_simulations)

    for ii in 1:number_of_simulations
        lattice = DLA_Simulation(lattice_dimensions)
        
        if display_each_lattice
            display_lattice(lattice)
        end

        fractal_masses[ii] = find_cluster_mass(lattice)
        fractal_length_scales[ii] = find_cluster_length_scale(lattice)
        
        
    end

    #return mean results
    return find__fractal_dimension(fractal_masses, fractal_length_scales), mean(fractal_masses), mean(fractal_length_scales)
end

rundla = run_DLA_and_analyze(lattice_dimensions, 1, true)
println("Fractal Dimension: ", rundla[1])
println("Mean Cluster Mass: ", rundla[2])
println("Mean Cluster Length Scale: ", rundla[3])

