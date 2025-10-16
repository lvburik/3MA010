##simulates brownian motion of a polymer chain using the rouse Modeling

### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------
### Functions and Structs
### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------


struct PhysicalParameters
    N::Int                  #number of particles
    D::Int                  #number of dimensions
    K::Float64              # spring constant (dimensionless)
    ζ::Float64              # friction coefficient (dimensionless)
    kT::Float64             # thermal energy (dimensionless)
end

struct SimulationParameters
    dt::Float64                 # time step
    equilibration_steps::Int    # number of equilibration_steps
    nsteps::Int                 # number of simulation steps
    save_nsteps::Int            # save r every N steps
    use_randomforce::Bool
end

function InitializeBeads(N::Int , D::Int)
    #initialize all beads at the origin
    bead_locations = zeros(D, N)

    #spread beads along positice x-axis
    for bead in 1:N
        bead_locations[1, bead] = bead-1
    end

    return bead_locations
end

function CalculatePotentialForces(bead_locations::Matrix{Float64}, D::Int, N::Int, K::Float64)
    #initialize potentials as 0
    F = zeros(D, N)

    l = bead_locations[:, 1:N-1] .- bead_locations[:, 2:N]  # 3 x (N-1) array

    # terminal particle 1
    F[:, 1] = -K * l[:, 1]

    # interior particles 2..N-1
    for i in 2:N-1
        F[:, i] = K * (l[:, i-1] - l[:, i])
    end

    # terminal particle N
    F[:, N] = K * l[:, N-1]

    return F
end

function GenerateRandomForces(kT::Float64, dt::Float64, ζ::Float64,D::Int, N::Int)
    return sqrt(2 * kT * dt / ζ) * randn(D, N)
end

function SimulationStep(bead_locations, physicalparameters)
    F = CalculatePotentialForces(bead_locations, physicalparameters.D, physicalparameters.N, physicalparameters.K)

        
    if simulationparameters.use_randomforce
        R = GenerateRandomForces(physicalparameters.kT, simulationparameters.dt, physicalparameters.ζ, physicalparameters.D, physicalparameters.N)
    else
        R = zeros(physicalparameters.D, physicalparameters.N)
    end

    #update and return new location
    return bead_locations .+ (simulationparameters.dt / physicalparameters.ζ) .* F .+ R
end

function Simulate(physicalparameters, simulationparameters)
    #initialize bead locations
    bead_locations = InitializeBeads(physicalparameters.N, physicalparameters.D)
    saved_bead_locations = []
    
    #equilibration loop
    
    for step in 1:simulationparameters.equilibration_steps
        bead_locations = SimulationStep(bead_locations, physicalparameters)
    end
    

    #save t = 0
    push!(saved_bead_locations, copy(bead_locations))

    #simulation loop
    for step in 1:simulationparameters.nsteps
        bead_locations = SimulationStep(bead_locations, physicalparameters)

        if step%simulationparameters.save_nsteps == 0
            push!(saved_bead_locations, copy(bead_locations))
        end
    end
    return saved_bead_locations
end

### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------
### Simulate
### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------

#define parameters
physicalparameters = PhysicalParameters(
    16,                     # number of particles
    3,                      # number of dimensions
    1.0,                    # spring constant (dimensionless)
    1.0,                    # friction coefficient (dimensionless)
    1.0,                    # thermal energy (dimensionless)
)

simulationparameters = SimulationParameters(
    .001,                   # timestep
    10000,                    # number of equilibration steps
    1000000,                   # number of steps
    500,                      # save location every n steps
    true,                   # use random force
)

#run and time simulation
@time bead_locations = Simulate(physicalparameters, simulationparameters)

### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------
### Analyze results
### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------
using Printf
using Plots

function CompareAvarageBondLenthSquared(kT, K, bead_locations)
    theoretical_average_bond_length_squared = 3 * kT / K

    # Calculate average bond length from simulation
    total_bond_length = 0.0
    num_bonds = 0
    
    # Loop through each timestep
    for positions in bead_locations
        # Calculate bond lengths for this timestep
        for ii in 1:size(positions,2)-1
            bond_vector = positions[:,ii] - positions[:,ii+1]
            total_bond_length += sum(bond_vector.^2)
            num_bonds += 1
        end
    end
    
    simulated_average_bond_length_squared = total_bond_length / num_bonds

    # Print results
    @printf("Theoretical average bond length^2: %.4f\n", theoretical_average_bond_length_squared)
    @printf("Simulated average bond length^2:   %.4f\n", simulated_average_bond_length_squared)
    @printf("Relative difference: %.2f%%\n", 100 * abs(theoretical_average_bond_length_squared - simulated_average_bond_length_squared) / theoretical_average_bond_length_squared)
    
    return theoretical_average_bond_length_squared, simulated_average_bond_length_squared
end

CompareAvarageBondLenthSquared(physicalparameters.kT, physicalparameters.K, bead_locations)

function PlotPolymerChain(bead_locations, timestep::Int=1)
    # Get the positions at the specified timestep
    positions = bead_locations[timestep]
    
    # Get the number of dimensions
    D = size(positions, 1)
    
    if D == 2
        # 2D plot
        p = plot(positions[1,:], positions[2,:],
                 label="Polymer chain",
                 marker=:circle,
                 markersize=4,
                 line=:solid,
                 xlabel="x",
                 ylabel="y",
                 title="Polymer Chain Configuration at timestep $timestep",
                 aspect_ratio=:equal)
                 
    elseif D == 3
        # 3D plot
        p = plot3d(positions[1,:], positions[2,:], positions[3,:],
                   label="Polymer chain",
                   marker=:circle,
                   markersize=4,
                   line=:solid,
                   xlabel="x",
                   ylabel="y",
                   zlabel="z",
                   title="Polymer Chain Configuration at timestep $timestep")
    else
        error("Can only plot 2D or 3D configurations")
    end
    
    return p
end

# Example usage:
# Plot initial configuration
p1 = PlotPolymerChain(bead_locations, 1)
display(p1)
savefig(p1, "polymer_initial.png")

# Plot final configuration
p2 = PlotPolymerChain(bead_locations, length(bead_locations))
display(p2)
savefig(p2, "polymer_final.png")

# Example usage:
# Plot initial configuration
fig1 = PlotPolymerChain(bead_locations, 1)
display(fig1)


# Plot final configuration
fig2 = PlotPolymerChain(bead_locations, length(bead_locations))
display(fig2)


