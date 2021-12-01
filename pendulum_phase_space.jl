include("./integrators.jl")
using DelimitedFiles

function pendulum_f(p::Float64, q::Float64; ω::Float64)
    """
        Pendulum's grad V. See Box 2.3 of Tremaine.
    """
    return ω^2 * sin(q)
end

# set up initial conditions
p0_array = range(-π, π, length = 50)
q0_array = range(-π, π, length = 50)

# integration time
n = 1000

# first, the modified Euler
println("Modified Euler")
io = open("output/pendulum_mod_euler.out", "w");

for p0 ∈ p0_array
    for q0 ∈ q0_array
        p, q = modified_euler(p0, q0, n, 1.0, pendulum_f, ω=1.0)
        writedlm(io, [p q])
    end
end

close(io)

# now, leapfrog
println("Leapfrog")
io = open("output/pendulum_leapfrog.out", "w");

for p0 ∈ p0_array
    for q0 ∈ q0_array
        p, q = leapfrog(p0, q0, n, 1.0, pendulum_f, ω=1.0)
        writedlm(io, [p q])
    end
end

close(io)

# finally, SS3
println("SS3")
io = open("output/pendulum_SS3.out", "w");
α_1 = 1.0/(2.0 - (2.0^(1/3) ))
αs = [α_1, 1-(2 * α_1), α_1]

for p0 ∈ p0_array
    for q0 ∈ q0_array
        p, q = SS3(p0, q0, αs, n, 0.75, pendulum_f, ω=1.0)
        writedlm(io, [p q])
    end
end

close(io)