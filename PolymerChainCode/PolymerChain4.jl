include("PolymerchainModel.jl")

using Plots
using Printf
using Statistics

#define parameters
physicalparameters = PhysicalParameters(
    16,                     # number of particles
    3,                      # number of dimensions
    3,                    # spring constant (dimensionless)
    1.0,                    # friction coefficient (dimensionless)
    1.0,                    # thermal energy (dimensionless)
)

simulationparameters = SimulationParameters(
    .01,                   # timestep
    2500,                    # number of equilibration steps
    10000,                   # number of steps
    1,                      # save location every n steps
    true,                   # use random force
)

#run and time simulation
@time bead_locations = Simulate(physicalparameters, simulationparameters)

#4a

function CalculateBondVectors(bead_locations)
    num_timesteps = length(bead_locations)
    D = size(bead_locations[1], 1)  # number of dimensions
    N = size(bead_locations[1], 2)  # number of beads

    #calculate bond vectors
    bond_vectors = Array{Float64}(undef, (num_timesteps, D, N-1))
    for (t_idx, positions) in enumerate(bead_locations)
        for bond_idx in 1:N-1
            bond_vectors[t_idx, :, bond_idx] = positions[:, bond_idx+1] - positions[:, bond_idx]
        end
    end

    return bond_vectors
end

bond_vectors = CalculateBondVectors(bead_locations)

function CalculateNormalModes(bond_vectors, p, dim = 1)
    #project onto x axis if not done so already
    if length(size(bond_vectors)) == 3
        bond_vectors = bond_vectors[:, dim, :]
    end

    (num_timesteps, N_bonds) = size(bond_vectors)
    
    if p > N_bonds
        throw("not a possible mode, p should be smaller than the number of beads")
    end

    #calculate normal modes
    normal_modes = zeros(num_timesteps)
    
    for (t_idx, bonds) in enumerate(eachrow(bond_vectors))
        for (bond_idx, bond) in enumerate(bonds)
            normal_modes[t_idx] += sin(pi * bond_idx * p / (N_bonds+1))bond
        end
    end
    normal_modes .*= sqrt(2/(N_bonds+1))

    return normal_modes
end

q = CalculateNormalModes(bond_vectors, 15)

function InnerProductNormalModes(bond_vectors)
    (num_timesteps, D, N) = size(bond_vectors)
    
    Q = Array{Float64}(undef, num_timesteps, N)

    for p = 1:N
        Q[:, p] = CalculateNormalModes(bond_vectors, p)
    end

    Q_abs = Array{Float64}(undef, num_timesteps) 
    for (timestep, q) in enumerate(eachrow(Q))
        Q_abs[timestep] = q'q

    end
    return Q_abs
end

Q_abs = InnerProductNormalModes(bond_vectors)

U_abs = [bond_vectors[ii, 1, :]'bond_vectors[ii, 1, :] for ii in 1:size(bond_vectors)[1]]

diffs = Q_abs - U_abs

#4b
function CorellateMode(bond_vectors, p)
    q = CalculateNormalModes(bond_vectors, p)
    
    (num_timesteps, D, N_bonds) = size(bond_vectors)

    #τ = physicalparameters.ζ/(2*physicalparameters.K(1-cos(pi*p/(N_bonds+1))))
    
    correlation = Array{Float64}(undef, num_timesteps)
    for t = 1:num_timesteps
        correlation[t] = 1/(num_timesteps-t)*q[1:num_timesteps-t]'q[t+1:num_timesteps]/(q'q/num_timesteps)
    end

    return correlation
end

@profview [display(plot(CorellateMode(bond_vectors, ii))) for ii in 1:14]
