include("PolymerchainModel.jl")

#define parameters
physicalparameters = PhysicalParameters(
    16,                     # number of particles
    3,                      # number of dimensions
    3,                    # spring constant (dimensionless)
    1.0,                    # friction coefficient (dimensionless)
    1.0,                    # thermal energy (dimensionless)
)

simulationparameters = SimulationParameters(
    .1,                   # timestep
    0,                    # number of equilibration steps
    100000,                   # number of steps
    1,                      # save location every n steps
    false,                   # use random force
)

#run and time simulation
@time bead_locations = Simulate(physicalparameters, simulationparameters)

PlotPolymerChain(bead_locations, length(bead_locations))
