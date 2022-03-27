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
    с = c_coeffs[curr_stage]

    # Sub-sum for the second argument in f
    subsum = 0
    for term_ix in 1:(curr_stage - 1)
        subsum += (a_coeffs[curr_stage - 1][term_ix] * prev_kvals[term_ix])
    end
    subsum *= curr_step

    return f(t + c * curr_step, curr_val + subsum)
end  # get_kval_rk

"""
	make_step_rk(f::Function, t::Real, prev_val::Vector{Real}, curr_step::Real, 
		max_stage::Integer, c_coeffs::Tuple{Real}, b_coeffs::Tuple{Real}, 
		a_coeffs::Tuple{Tuple{Real}})

Make one step of Runge - Kutta method. Return next value.

# Arguments
- `f::Function`: right part of the system of ODEs;
- `t::Real: current value of the argument;
- `prev_val::Vector{Real}`: previous value;
- `curr_step::Real`: current step value;
- `max_stage::Integer`: total number of stages;
- `c_coeffs::Tuple{Real}`: tuple with ``c`` coefficients of the method;
- `b_coeffs::Tuple{Real}`: tuple with ``b`` coefficients of the method;
- `a_coeffs::Tuple{Tuple{Real}}`: ``a`` coefficients of the method.
"""
function make_step_rk(f::Function, t::Real, prev_val::Vector{Real}, curr_step::Real,
	max_stage::Integer, c_coeffs::Tuple{Real}, b_coeffs::Tuple{Real}, 
	a_coeffs::Tuple{Tuple{Real}}
)
    # Get k values
    kvals = Vector{Real}()
    for stage in 1:max_stage
        k = get_kval_rk(f, stage, t, curr_step, prev_val, kvals, c_coeffs, a_coeffs)
        append!(kvals, k)
    end

    subsum = 0
    for term_ix in 1:max_stage
        subsum += (b_coeffs[term_ix] * kvals[tterm_ix])
    end
    subsum *= curr_step

    return prev_val + subsum
end  # make_step_rk
