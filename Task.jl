"""
	f1(t::Real, y::Vector{Real})

Return value of ``f_1`` by given t and ``(y_1, y_2)^T``.

# Arguments
- `t::Real`: current value of the argument;
- `y::Vector{Real}`: vector of function values (2 components).
"""
function f1(t::Real, y::Vector{Real})
    max_val = max(y[2], 10^(-3))
    return 2 * t * y[1] * log(max_val)
end  # f1

"""
	f2(t::Real, y::Vector{Real})

Return value of ``f_2`` by given t and ``(y_1, y_2)^T``.

# Arguments
- `t::Real`: current value of the argument;
- `y::Vector{Real}`: vector of function values (2 components).
"""
function f2(t::Real, y::Vector{Real})
    max_val = max(y[1], 10^(-3))
    return -2 * t * y[2] * log(max_val)
end  # f2

"""
	f(t::Real, y::Vector{Real})

Return value of ``f`` as vector of 2 values by given t and ``(y_1, y_2)^T``.

# Arguments
- `t::Real`: current value of the argument;
- `y::Vector{Real}`: vector of function values (2 components).
"""
function f(t::Real, y::Vector{Real})
    f_vec = Vector{Real}()
    append!(f_vec, f1(t, y))
    append!(f_vec, f2(t, y))
    return f_vec
end  # f

"""
    coeffs_rk()

Return tuple with tuples of test method coefficients. Structure of returned tuple:
(c, b, a).
"""
function coeffs_rk()
    c_coeffs = (0, 1 / 2, 1 / 2, 1)
    b_coeffs = (1 / 6, 2 / 6, 2 / 6, 1 / 4)
    a_coeffs = ((1 / 2), (0, 1 / 2), (0, 0, 1))

    return (c_coeffs, b_coeffs, a_coeffs)
end  # coeffs_rk
