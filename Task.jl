include("CustomLogger.jl")

"""
	f1(t::Float64, y::Vector{Float64})

Return value of ``f_1`` by given t and ``(y_1, y_2)^T``.

# Arguments
- `t::Float64`: current value of the argument;
- `y::Vector{Float64}`: vector of function values (2 components).
"""
function f1(t::Float64, y::Vector{Float64})
    max_val = max(y[2], 10^(-3))
    return 2 * t * y[1] * log(max_val)
end  # f1

"""
	f2(t::Float64, y::Vector{Float64})

Return value of ``f_2`` by given t and ``(y_1, y_2)^T``.

# Arguments
- `t::Float64`: current value of the argument;
- `y::Vector{Float64}`: vector of function values (2 components).
"""
function f2(t::Float64, y::Vector{Float64})
    max_val = max(y[1], 10^(-3))
    return -2 * t * y[2] * log(max_val)
end  # f2

"""
	f(t::Float64, y::Vector{Float64})

Return value of ``f`` as vector of 2 values by given t and ``(y_1, y_2)^T``.

# Arguments
- `t::Float64`: current value of the argument;
- `y::Vector{Float64}`: vector of function values (2 components).
"""
function f(t::Float64, y::Vector{Float64})
    f_vec = Vector{Float64}()
    append!(f_vec, f1(t, y))
    append!(f_vec, f2(t, y))
    return f_vec
end  # f

"""
    coeffs_rk()

Return tuple with vectors of coefficients for 4-th order Runge - Kutta method. Structure of 
returned tuple: (c, b, a).
"""
function coeffs_rk()
    c_coeffs = [0, 1 / 2, 1 / 2, 1]
    b_coeffs = [1 / 6, 2 / 6, 2 / 6, 1 / 6]
    a_coeffs = [[1 / 2], [0, 1 / 2], [0, 0, 1]]

    return (c_coeffs, b_coeffs, a_coeffs)
end  # coeffs_rk

"""
    coeffs_dp()

Return tuple with vectors of coefficients for Dormand - Prince method. Structure of returned 
tuple: (c, b, bhat, a).
"""
function coeffs_dp()
    c_coeffs = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1]
    b_coeffs = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0]
    bhat_coeffs = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 
        1 / 40]
    a_coeffs = [[1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9], 
                [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
                [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
                [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]]

    return (c_coeffs, b_coeffs, bhat_coeffs, a_coeffs)
end  # coeffs_dp

"""
    exact_sol(t::Float64)

Returns exact solution as a vector.

# Arguments
`t::Float64`: point in which ``y(t)`` should be returned.
"""
function exact_sol(t::Float64)
    sol = Vector{Float64}()

    append!(sol, exp(sin(t^2)))
    append!(sol, exp(cos(t^2)))

    return sol
end # exact_sol
