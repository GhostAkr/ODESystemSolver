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

"""
    maxnorm(vec::Vector{Real})

Computes 'max' norm of given vector `vec`. Max norm can be calculated as follows:
``\\underset{i = 1, \\dots n}{\\text{max}} |v_i|``. Return norm of `vec` as one real
number.

# Arguments
`vec::Vector{Real}`: vector which norm should be calculated.
"""
function maxnorm(vec::Vector{Real})
    if isempty(vec)
        @error("Trying to calculate norm of empty vector")
    end

    max_component = abs(vec[1])
    for component in vec
        abs_comp = abs(component)

        if abs_comp > max_component
            max_component = abs_comp
        end
    end

    return max_component
end  # maxnorm

"""
	solve_rk(f::Function, t_limits::Tuple{Real}, initial_step::Real, 
		initial_vals::Tuple{Real}, max_stage::Integer, c_coeffs::Tuple{Real}, 
		b_coeffs::Tuple{Real}, a_coeffs::Tuple{Tuple{Real}})

Solve system of ODEs using Runge - Kutta method. System of ODE can be represented as
follows: ``y' = f(t, y), y(t_0) = y_0`` where y is a vector. To achieve more accurancy
this method uses Richardson's extrapolation to control step during the computational 
process.

Return dictionary in following format: `<t>: <y(t)`.

Function `f` takes ``t`` and vector ``y`` and returns vector of real values which length
is equal to the length of ``y``.

``c`` and ``b`` coefficients are stored as usual tuples where each element of tuple
denotes coefficient for appropriate stage. ``a`` coefficients are represented as nested 
tuples. Each internal tuple contains coeffcients of one method's stage. For example for 
4-staged method we consider following structure of ``a`` coefficients:
((a21), (a31, a32), (a41, a42, a43)). Length of tuples with coefficients should 
correspond to `max_stage` value.

# Arguments
- `f::Function`: right part of the system of ODEs;
- `t_limits::Tuple{Real}`: system argument limits;
- `initial_step::Real`: initial step;
- `initial_val::Tuple{Real}`: initial values (``y_0``);
- `tol::Real`: target tolerance;
- `max_stage::Integer`: total number of stages;
- `c_coeffs::Tuple{Real}`: tuple with``c`` coefficients of the method;
- `b_coeffs::Tuple{Real}`: tuple with ``b`` coefficients of the method;
- `a_coeffs::Tuple{Tuple{Real}}`: ``a`` coefficients of the method.
"""
function solve_rk(f::Function, t_limits::Tuple{Real}, initial_step::Real, 
	initial_val::Tuple{Real}, tol::Real, max_stage::Integer, c_coeffs::Tuple{Real},
	b_coeffs::Tuple{Real}, a_coeffs::Tuple{Tuple{Real}}
)
    curr_val = initial_val
    curr_step = initial_step
    t = t_limits[1]  # Current argument value

    solution = Dict{Real, Real}()

    while t < t_limits[2]
        near_end = false
        if t + curr_step * 2 > t_limits[2]
            near_end = true
        end

        if near_end
            curr_step = (t_limits[2] - t) / 2
        end


        while (true)
            # Make 2 usual steps
            y1 = make_step_rk(f, t, curr_val, curr_step, max_stage, c_coeffs, b_coeffs, 
                a_coeffs)
            y2 = make_step_rk(f, t + curr_step, y1, curr_step, max_stage, c_coeffs, 
                b_coeffs, a_coeffs)
                
            # If we are near the end of interval then we don't need to perform Richardson's 
            # extrapolation
            if (near_end)
                solution[t + curr_step] = y1
                solution[t + curr_step * 2] = y2
                return solution
            end

            # Make 1 double step
            w = make_step_rk(f, t, curr_val, curr_step * 2, max_stage, c_coeffs, b_coeffs, 
                a_coeffs)

            err = 1 / (2^p - 1) * maxnorm(y2 - w)

            fac = 0.9
            facmin = 0.5
            facmax = 3
            new_step = curr_step * min(facmax, max(facmin, fac * 
                (tol / err)^(1 / (max_stage))))

            if err <= tol
                t += (curr_step * 2)
                curr_val = y2
                curr_step = new_step
                break
            else
                curr_step = new_step
            end
        end

        solution[t - curr_step] = y1
        solution[t] = y2
    end

    return solution
end  # solve_rk
