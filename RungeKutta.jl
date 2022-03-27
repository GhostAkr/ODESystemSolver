"""
	get_kval_rk(curr_stage::Integer, t::Real, curr_step::real, curr_val::Real, 
		prev_kvals::Vector{Real}, c_coeffs::Tuple{Real}, a_coeffs::Tuple{Real})

Return value of ``k_i`` for given ``i`` (`curr_stage`).

No additional checks are performed for efficiency (assuming that these checks are made
inside `solve_rk()` function.

# Arguments
- `f::Function`: right part of the system of ODEs;
- `curr_stage::Integer`: current stage (``i`` value);
- `t::Real`: current value of the argument;
- `curr_step::Real`: current step value;
- `curr_val::Real`: current value;
- `prev_kvals::Real`: values of ``k`` on previous stages;
- `c_coeffs::Tuple{Real}`: tuple with ``c`` coefficients of the method;
- `a_coeffs::Tuple{Tuple{Real}}`: ``a`` coefficients of the method.
"""
function get_kval_rk(f::Function, curr_stage::Integer, t::Real, curr_step::Real, 
	curr_val::Real, prev_kvals::Vector{Real}, c_coeffs::Tuple{Real}, 
    a_coeffs::Tuple{Tuple{Real}}
)
    # Get relevant c coefficient
    c = 0
    if curr_stage > 1
        с = c_coeffs[curr_stage - 1]
    end

    # Sub-sum for the second argument in f
    subsum = 0
    for term_ix in 1:(curr_stage - 1)
        subsum += (a_coeffs[curr_stage - 1][term_ix] * prev_kvals[term_ix])
    end
    subsum *= curr_step

    return f(t + c * curr_step, curr_val + subsum)
end  # get_kval_rk
