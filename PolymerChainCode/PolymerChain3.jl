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

#plot final configuration to check for stability
finalconfigplot = PlotPolymerChain(bead_locations, length(bead_locations))


#3d

#part 1: calculate bond lentghs over time and plot:
function CalculateBondLengths(bead_locations, show_plot = false)
    # Calculate instantaneous bond lengths for all timesteps
    num_timesteps = length(bead_locations)
    N = size(bead_locations[1], 2)  # number of beads
    num_bonds = N - 1
    
    # Initialize array to store bond lengths: (num_bonds x num_timesteps)
    bond_lengths = zeros(num_bonds, num_timesteps)
    
    # Loop through each timestep
    for (t_idx, positions) in enumerate(bead_locations)
        # Calculate bond lengths for this timestep
        for bond_idx in 1:num_bonds
            bond_vector = positions[:, bond_idx] - positions[:, bond_idx+1]
            bond_lengths[bond_idx, t_idx] = sqrt(sum(bond_vector.^2))
        end
    end
    
    if show_plot == true
        # Plot 5 random bonds over time
        num_bonds_to_plot = min(2, num_bonds)
        random_bonds = rand(1:num_bonds, num_bonds_to_plot)
        
        timesteps = 0:length(bead_locations)-1
        
        p = plot(xlabel="Timestep", ylabel="Bond Length", 
                 title="Bond Lengths Over Time", legend=:best,
                 framestyle=:box, tick_direction=:in)
        
        for bond_idx in random_bonds
            plot!(p, timesteps, bond_lengths[bond_idx, :], 
                  label="Bond $bond_idx - $(bond_idx+1)", alpha=0.7)
        end
        
        display(p)
        return bond_lengths, p
    end
    
    return bond_lengths
end

bond_lengths, bondlenthsplot = CalculateBondLengths(bead_locations, true)
savefig(bondlenthsplot, "bond_length_fluctuations.png")

#part 2: calculate bondlength squared and compare
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

temp = CompareAvarageBondLenthSquared(physicalparameters.kT, physicalparameters.K, bead_locations)

#3e

function CalculateSquaredEndToEndLength(bead_locations, makeplot = false)
    num_timesteps = length(bead_locations)

    end_to_end_length_squared = Array{Float64}(undef, num_timesteps)
    for (t_idx, positions) in enumerate(bead_locations)
        length_vector = positions[:, end] - positions[:,1]
        end_to_end_length_squared[t_idx] = sum(length_vector.^2)
    end
    
    if makeplot==false
        return end_to_end_length_squared
    end

    p = plot(xlabel="Timestep", ylabel="Square end to end Length", 
                 title="End To End Length Over Time", legend=:best,
                 framestyle=:box, tick_direction=:in)

    plot!(p, end_to_end_length_squared, label = "N=$(physicalparameters.N)")
    display(p)
    return end_to_end_length_squared, p
end

sq_end2end, end2endplot = CalculateSquaredEndToEndLength(bead_locations, true)
mean_sq_end2end = mean(sq_end2end)
std_sq_end2end = std(sq_end2end)
savefig(end2endplot, "end2endplot.png")

