##simulates brownian motion of a polymer chain using the rouse Modeling

### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------
### Functions and Structs
### ---------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------
using Plots

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






