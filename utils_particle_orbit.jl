# Define some constants needed for later calculations

# tiny number for orbit calculation
ϵ = 1.E-308

# Below this inclination, the broken angles ϖ and θ 
# equal the corresponding unbroken angles to within 
# machine precision, so a practical boundary 
# for planar orbits
MIN_INC = 1.e-8

# Below this eccentricity, corrections at order e^2 
# are below machine precision, so we use stable 
# expressions accurate to O(e) for the mean longitude 
# below for near-circular orbits.
MIN_ECC = 1.e-8

function acos2(num::Float64, denom::Float64, disambiguator::Float64)::Float64
    """
        Finds the arccos within [0, π] using disambiguator
        to tell which quadrant the angle is in.

        Returns 0 or π appropriately if num > denom
        Returns 0 if denom == 0.0

        From https://github.com/hannorein/rebound/blob/main/src/tools.c#L935
    """

    cosine = num/denom

    if cosine > -1.0 && cosine < 1.0
        # get the normal acos
        val = acos(cosine)
        
        # if disambiguator
        # use the value on the other side
        if disambiguator < 0
            val = -val
        end
    else
        # returns π if cosine <= -1.0
        # 0 other wise
        val = ( cosine <= -1.0 ) ? π : 0.0
    end

    return val
end

function mod2pi(f::Float64)::Float64
    """
        Puts angles into [0, 2π] range.
        From https://github.com/hannorein/rebound/blob/main/src/tools.c#L457
    """
    return mod(2*π + mod(f, 2*π), 2*π)
end

function orbit_to_particle(;G::Float64, primary::particle, m::Float64, a::Float64, 
    e::Float64, inc::Float64, Ω::Float64, ω::Float64, 
    f::Float64)
    """
        Gets Cartesian position and velocity from orbital elements.

        From https://github.com/hannorein/rebound/blob/main/src/tools.c#L838

        I try to keep the flow here same as the implementation in REBOUND
        for easier comparison. Some comments are copied over for clarifications.
    """
    
    # error code
    # 0 is no error
    err = UInt(0)

    if e == 1.0
        # radial orbit does not work
        err = 1
        return particle_NaN, err
    end
    
    if e < 0.0
        # no negative eccentricity either
        err = 2
        return particle_NaN, err
    end
    
    if e > 1.0
        # if bound orbit (a > 0)
        # eccentricity must be < 1
        if a > 0.0            
            err = 3
            return particle_NaN, err
        end
    else
        # if unbound orbit (a < 0)
        # eccentricity must be > 1
        if a < 0.0
            err = 4
            return particle_NaN, err
        end
    end
    
    if e*cos(f) < -1.0
        err = 5
        return particle_NaN, err
    end
    
    if primary.m < ϵ
        # if primary particle has no mass
        err = 6
        return particle_NaN, err
    end
    
    # create new particle
    p = particle()
    p.m = m

    r = a * ( 1 - e^2 ) / ( 1 + e * cos(f) )
    v0 = sqrt( G * ( m + primary.m ) / a / ( 1. - e^2 ) )
    
    cO = cos(Ω)
    sO = sin(Ω)
    
    co = cos(ω)
    so = sin(ω)
    
    cf = cos(f)
    sf = sin(f)
    
    ci = cos(inc)
    si = sin(inc)
    
    # Murray & Dermott 2.122
    p.x = primary.x + r * ( cO * (co*cf-so*sf) - sO * (so*cf+co*sf) * ci)
    p.y = primary.y + r * ( sO * (co*cf-so*sf) + cO * (so*cf+co*sf) * ci)
    p.z = primary.z + r * ( so*cf + co*sf) * si

    # Murray & Dermott 2.36 after applying the 3 rotation matrices 
    # from Sec. 2.8 to the velocities in the orbital plane
    p.vx = primary.vx + v0 * ((e+cf) * (-ci*co*sO - cO*so) - sf * (co*cO - ci*so*sO))
    p.vy = primary.vy + v0 * ((e+cf) * (ci*co*cO - sO*so)  - sf * (co*sO + ci*so*cO))
    p.vz = primary.vz + v0 * ((e+cf) * co * si - sf * si * so)

    p.ax = 0
    p.ay = 0
    p.az = 0
    
    return p, err
    
end

function particle_to_orbit(;G::Float64, primary::particle, p::particle)
    """
        Gets orbital elements from particle position and velocity.
        Note: T (time of pericenter passage) is not set here
            since it needs a time t (c.f. lines 335 - 337 here 
            versus 1073 - 1078 in REBOUND).
        Note: needs a primary particle.

        From https://github.com/hannorein/rebound/blob/main/src/tools.c#L950

        I try to keep the flow here same as the implementation in REBOUND
        for easier comparison. Some comments are copied over for clarifications.
    """
    err = UInt8(0)
    o = orbit()

    if primary.m < ϵ
        # if the primary has no mass
        err = 1
        return orbit_NaN, err
    end

    # reduced mass
    mu = G * ( p.m + primary.m )

    # distance to primary
    dx = p.x - primary.x
    dy = p.y - primary.y
    dz = p.z - primary.z
    
    # velocity rel. to primary
    dvx = p.vx - primary.vx
    dvy = p.vy - primary.vy
    dvz = p.vz - primary.vz

    # angular momentum of particle p
    hx = (dy*dvz - dz*dvy)
    hy = (dz*dvx - dx*dvz)
    hz = (dx*dvy - dy*dvx)

    # set orbital distance
    o.d = sqrt( dx^2 + dy^2 + dz^2 )

    vsquared = dvx^2 + dvy^2 + dvz^2
    vcircsquared = mu/o.d

    # set orbital info 
    # (velocity, acc, hill radius, angular momentum)
    o.v = sqrt(vsquared)
    o.a = - mu / ( vsquared - 2.0 * vcircsquared )	
    o.rhill = o.a * cbrt( p.m / (3.0 * primary.m) )
    o.h = sqrt( hx^2 + hy^2 + hz^2 )

    if o.d <= ϵ
        # if the particle as at the same location as primary	
        err = 2
        return reb_orbit_nan, err
    end

    vdiffsquared = vsquared - vcircsquared
    vr = ( dx*dvx + dy*dvy + dz*dvz )/o.d
    rvr = o.d * vr
    muinv = 1.0 / mu

    # get eccentricity
    ex = muinv * ( vdiffsquared*dx - rvr*dvx )
    ey = muinv * ( vdiffsquared*dy - rvr*dvy )
    ez = muinv * ( vdiffsquared*dz - rvr*dvz )

    # set eccentricity, mean motion, period
    o.e = sqrt( ex^2 + ey^2 + ez^2 )
    # n and P are negative if hyperbolic
    o.n = o.a/abs(o.a)*sqrt( abs(mu / (o.a*o.a*o.a)) )
    o.P = 2*π/o.n

    # ̂h ⋅ ̂z = cos(i) ∈ [0,π] so use acos2
    o.inc = acos2(hz, o.h, 1.0)

    # vector pointing along ascending node
    # ̂z × h
    nx = -hy
    ny =  hx
    n = sqrt( nx^2 + ny^2 )

    # ̂x ⋅ ̂n = cos(Ω) (0 if i == 0)
    # use y components as disambiguator
    # since Ω, ϖ, θ are all measured from
    # the x axis
    o.Ω = acos2(nx, n, ny)

    if o.e < 1.0
        # eccentric anomaly
        # if vr < 0, apocenter → pericenter
        # ⟹ ea ∈ [π, 2π] ⟹ ea = -acos(cos(ea))
        ea = acos2(1.0 - o.d/o.a, o.e, vr)
        # mean anomaly (Kepler's equation)
        o.M = ea - o.e*sin(ea)
    else
        ea = acosh( (1.0 - o.d / o.a) / o.e)
        if vr < 0.0
            # approaching pericenter ⟹ negative ea
            ea = -ea
        end
        # hyperbolic Kepler's
        o.M = o.e*sinh(ea) - ea
    end

    # if near-planar, true longitude is well-defined for position
    # and ϖ for pericenter if e ≠ 0
    # so get those, and get the other angles from them
    if o.inc < MIN_INC || o.inc > π - MIN_INC
        # ̂x ⋅ ̂r = cos(θ) (true longitude)
        o.θ = acos2(dx, o.d, dy)
        # ̂x ⋅ ̂e = cos(ϖ) (longitude of periapsis)
        o.ϖ = acos2(ex, o.e, ey)

        if o.inc < π/2.0
            o.ω = o.ϖ - o.Ω
            o.f = o.θ - o.ϖ

            # calculate the mean longitude
            if o.e > MIN_ECC
                # ϖ is well-defined when e ≠ 0
                # can get the mean longitude
                o.l = o.ϖ + o.M
            else
                # otherwise, e << 1 and ϖ not well-defined
                # use l = θ+(M-f) since M-f is O(e) so well behaved
                # M-f from Murray & Dermott Eq 2.93. 
                # This way l → θ smoothly as e → 0
                o.l = o.θ - 2.0*o.e*sin(o.f)
            end
        else
            o.ω = o.Ω - o.ϖ
            o.f = o.ϖ - o.θ

            # get mean longitude
            # same as above, but retrograde changes sign
            if o.e > MIN_ECC
                o.l = o.ϖ - o.M
            else
                o.l = o.θ + 2.0*o.e*sin(o.f)
            end

        end

    # if non-planar, can't get broken angles like above
    # ω+f is always well-defined.
    # ω well-defined if e ≠ 0
    else
        # omega plus f. 
        # Both angles measured in orbital plane. Well defined if i ≠ 0.
        ωpf = acos2(nx*dx + ny*dy, n*o.d, dz)
        o.ω = acos2(nx*ex + ny*ey, n*o.e, ez)
        
        # calculate the mean longitude
        if o.inc < π/2.0
            o.ϖ = o.Ω + o.ω
            o.f = ωpf - o.ω
            o.θ = o.Ω + ωpf

            if o.e > MIN_ECC
                # ϖ is well-defined e ≠ 0
                o.l = o.ϖ + o.M
            else
                # use l = θ+(M-f). M-f is O(e) so well behaved
                # Murray & Dermott eqn 2.93
                o.l = o.θ - 2.0*o.e*sin(o.f)
            end

        # same as above, but retrograde changes sign
        else
            o.ϖ = o.Ω - o.ω
            o.f = ωpf - o.ω
            o.θ = o.Ω - ωpf

            if o.e > MIN_ECC
                # ϖ is well-defined e ≠ 0
                o.l = o.ϖ - o.M
            else
                # same as above
                o.l = o.θ + 2.0*o.e*sin(o.f)
            end
        end
    end

    # code to set T (time of pericenter passage) here
    # need time t. M = n(t-T).
    # o.T = t - o.M/abs(o.n)

    # move angles to [0, 2π]
    o.f = mod2pi(o.f)
    o.l = mod2pi(o.l)
    o.M = mod2pi(o.M)
    o.θ = mod2pi(o.θ)
    o.ω = mod2pi(o.ω)

    return o, err
end