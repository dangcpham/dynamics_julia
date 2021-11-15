mutable struct particle
    # mass
    m::Float64
    
    # cartesian position, velocity, acceleration
    x::Float64
    y::Float64
    z::Float64
    
    vx::Float64
    vy::Float64
    vz::Float64
    
    ax::Float64
    ay::Float64
    az::Float64
    
    particle() = new()
end

particle_NaN = particle()
particle_NaN.m = NaN

mutable struct orbit
    # distance from primary
    d::Float64
    # velocity rel. to primary
    v::Float64
    # mutual hill radius
    rhill::Float64
    # semi-major axis
    a::Float64
    # angular momentum
    h::Float64

    # eccentricity
    e::Float64
    # mean motion
    n::Float64
    # period
    P::Float64
    # inclination
    inc::Float64

    # longitude of ascending node
    Ω::Float64
    # mean anomaly
    M::Float64
    # true longitude
    θ::Float64
    # longitude of periapsis
    ϖ::Float64
    # argument of perihelion
    ω::Float64
    # true anomaly
    f::Float64
    # mean longitude
    l::Float64


    # time of pericenter passage
    T::Float64

    orbit() = new()

end

orbit_NaN = orbit()
orbit_NaN.d = NaN