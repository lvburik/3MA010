
#Modeling exercise 2

using LinearAlgebra
using Statistics

#------------------
#--- Functions ---
#------------------

function Initialize(N)
    R = zeros(3, N) #particle positions
    # initialize chain in a straight line along x-axis
    for i in 1:N
        R[:, i] = [(i-1), 0.0, 0.0]
    end
    return R
end

function calculate_potential_forces(R, params)
    N = params.N
    K = params.K

    F = zeros(3, N)

    l = R[:, 1:N-1] .- R[:, 2:N]  # 3 x (N-1) array
    
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

function Simulate(params)
    r = Initialize(params.N)
    savedr = [r]

    for step in 1:params.nsteps
        #potential force term
        F = calculate_potential_forces(r, params)

        # random force term
        R = sqrt(2 * params.kT * params.dt / params.ζ) * randn(3, params.N)

        # update positions
        r .= r .+ (params.dt / params.ζ) .* F .+ R

        if step%params.save_nsteps == 0
            push!(savedr, r)
        end
    end

    return savedr 
end

#------------------
#--- Parameters ---
#------------------
struct parameters
    N::Int #number of particles
    K::Float64            # spring constant (dimensionless)
    ζ::Float64            # friction coefficient (dimensionless)
    kT::Float64            # thermal energy (dimensionless)
    dt::Float64            # time step
    nsteps::Int      # number of simulation steps
    save_nsteps::Int     # save r every N steps
end
#------------------
#---- Simulate ----
#------------------
params = parameters(16, 1.0, 1.0, 1.0, 0.01, 100000, 1)
@time savedr = Simulate(params)

#------------------
#- Accuracy Tests -
#------------------

function mean_bond_length_squared(savedr)
    total = 0.0
    count = 0
    for r in savedr
        l = r[:, 1:end-1] .- r[:, 2:end]
        total += sum(norm.(eachcol(l)).^2)
        count += size(l, 2)
    end
    return total / count
end

println("⟨l²⟩ =", mean_bond_length_squared(savedr))
println("Expected =", 3 * params.kT / params.K)

function mean_square_displacement(savedr, dt)
    nframes = length(savedr)
    r0 = savedr[1]
    msd = Float64[]
    for i in 2:nframes
        Δr = savedr[i] - r0
        push!(msd, mean(sum(Δr.^2, dims=1)))  # average over particles
    end
    t = dt * collect(1:nframes-1)
    return t, msd
end

t, msd = mean_square_displacement(savedr, params.dt * params.save_nsteps)
println("Diffusion slope (approx) =", msd[end]/t[end])
println("Expected D =", params.kT / params.ζ)

function average_potential_energy(savedr, K)
    totalU = 0.0
    count = 0
    for r in savedr
        l = r[:, 1:end-1] .- r[:, 2:end]
        totalU += 0.5 * K * sum(norm.(eachcol(l)).^2)
        count += 1
    end
    return totalU / count
end

println("Average U =", average_potential_energy(savedr, params.K))
println("Expected =", 1.5 * (params.N - 1) * params.kT)





