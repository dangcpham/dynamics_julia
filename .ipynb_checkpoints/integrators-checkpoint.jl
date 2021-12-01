function modified_euler(q0::Float64, p0::Float64, 
    n::Int, h::Float64, f::Function; args...)

    """ Modified Euler integrator.

    Parameters:
        p0, q0: initial momentum and position
        n: how long to integrate
        h: time step size
        f: function to change momentum (grad V)
        args: arguments for f

    Returns:
        (momentum p array, position q array)

    """

    # initialize arrays for p,q
    p = Array{Float64}(undef, n)
    q = Array{Float64}(undef, n)
    
    # set the initial condition
    p[1], q[1] = p0, q0
    
    # perform modified euler steps
    for i ∈ 1:n-1
        p[i+1] = p[i] - h*f(p[i], q[i]; args...)
        q[i+1] = q[i] + h*p[i+1]
    end

    return p, q
end

function kdk_base(p0::Float64, q0::Float64, hprime::Float64,
    f::Function; args...)

    """ Base kick-drift-kick function. Can be used with leapfrog or other
    composition method (e.g. SS3). 

    Parameters:
        p0, q0: initial momentum and position
        hprime: time step size
        f: function to change momentum (grad V)
        args: arguments for f

    Returns:
        (new momentum p, new position q)
    """
    # kick 1/2
    phalf = p0 - (hprime/2.)*f(p0, q0; args...)
    # drift 1
    qnew  = q0 + hprime*phalf
    # kick 1/2
    pnew  = phalf - (hprime/2.)*f(phalf, qnew; args...)

    return pnew, qnew
end

function leapfrog(p0::Float64, q0::Float64, n::Int, h::Float64,
    f::Function; args...)

    """ Leapfrog integrator.

    Parameters:
        p0, q0: initial momentum and position
        n: how long to integrate
        h: time step size
        f: function to change momentum (aka force, grad V)
        args: arguments for f

    Returns:
        (momentum p array, position q array)

    """

    # initialize arrays for p, q
    p = Array{Float64}(undef, n)
    q = Array{Float64}(undef, n)

    # set initial p, q
    p[1], q[1] = p0, q0

    for i ∈ 1:n-1
        p[i+1], q[i+1] = kdk_base(p[i], q[i], h, f; args...)
    end

    return p, q

end

function SS3(p0::Float64, q0::Float64, 
    αs::Array, n::Int, h::Float64,
    f::Function; args...)

    """ Leapfrog composition integrator SS3[4] (c.f. Blanes & Cases, page 89)

    Parameters:
        p0, q0: initial momentum and position
        αs: list of α for each composition step
        n: how long to integrate
        h: time step size
        f: function to change momentum (aka force, grad V)
        args: arguments for f

    Returns:
        (momentum p array, position q array)

    """
    # initialize arrays for p, q
    p = Array{Float64}(undef, n)
    q = Array{Float64}(undef, n)

    # set initial p, q
    p[1], q[1] = p0, q0

    for i ∈ 1:n-1
        # get current p, q
        pnew, qnew = p[i], q[i]

        # do leapfrog kick-drift-kick for each alpha
        for α ∈ αs
            pnew, qnew = kdk_base(pnew, qnew, α*h, f; args...)
        end

        # update next p, q
        p[i+1], q[i+1] = pnew, qnew
    end

    return p, q
end
        