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

Return tuple with vectors of test method coefficients. Structure of returned tuple:
(c, b, a).
"""
function coeffs_rk()
    c_coeffs = [0, 1 / 2, 1 / 2, 1]
    b_coeffs = [1 / 6, 2 / 6, 2 / 6, 1 / 4]
    a_coeffs = [[1 / 2], [0, 1 / 2], [0, 0, 1]]

    return (c_coeffs, b_coeffs, a_coeffs)
end  # coeffs_rk

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
