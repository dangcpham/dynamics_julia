include("./particle_struct.jl")
include("./utils.jl")

######################################################################

# gravitational constant
G = 1.0

# set up the primary particle (i.e. star in the middle)
primary = particle()
primary.m = 1.0

# orbital elements for the secondary particle
m = 3.0
a, e, inc, Ω, ω, f = 1.2, 0.5, 0.3, 0.2, 0.1, 1.5

######################################################################

println("Initial orbital elements:")
println("a: ", a)
println("e: ", e)
println("i: ", inc)
println("Ω: ", Ω)
println("ω: ", ω)
println("f: ", f)

secondary, err = orbit_to_particle(G=G, primary=primary,
                                    m=m, a=a, e=e, inc=inc, 
                                    Ω=Ω, ω=ω, f=f)
println("\nParticle's position and velocity from orbital elements:")
println("Error code ", err)
println("x: ", secondary.x)
println("y: ", secondary.y)
println("z: ", secondary.z)
println("vx: ", secondary.vx)
println("vy: ", secondary.vy)
println("vz: ", secondary.vz)

orbit_out, err = particle_to_orbit(G=G, primary=primary, 
                                          p=secondary)
println("\nOrbital elements from particle's position and velocity:")
println("Error code ", err)
println("a: ", orbit_out.a)
println("e: ", orbit_out.e)
println("i: ", orbit_out.inc)
println("Ω: ",orbit_out.Ω)
println("ω: ",orbit_out.ω)
println("f: ",orbit_out.f)