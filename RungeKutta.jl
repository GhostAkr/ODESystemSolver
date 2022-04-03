using DelimitedFiles

include("Task.jl")

"""
	get_kval_rk(curr_stage::Integer, t::Float64, curr_step::real, curr_val::Vector{Float64}, 
		prev_kvals::Vector{Vector{Float64}}, c_coeffs::Vector{Float64}, 
        a_coeffs::Vector{Float64})

Return value of ``k_i`` for given ``i`` (`curr_stage`).

No additional checks are performed for efficiency (assuming that these checks are made
inside `solve_rk()` function.

# Arguments
- `f::Function`: right part of the system of ODEs;
- `curr_stage::Integer`: current stage (``i`` value);
- `t::Float64`: current value of the argument;
- `curr_step::Float64`: current step value;
- `curr_val::Vector{Float64}`: current value;
- `prev_kvals::Vector{Vector{Float64}}`: values of ``k`` on previous stages;
- `c_coeffs::Vector{Float64}`: tuple with ``c`` coefficients of the method;
- `a_coeffs::Vector{Vector{Float64}}`: ``a`` coefficients of the method.
"""
function get_kval_rk(f::Function, curr_stage::Integer, t::Float64, curr_step::Float64, 
	curr_val::Vector{Float64}, prev_kvals::Vector{Vector{Float64}}, 
    c_coeffs::Vector{Float64}, a_coeffs::Vector{Vector{Float64}}
)
    # Get relevant c coefficient
    с = c_coeffs[curr_stage]

    # Sub-sum for the second argument in f
    subsum = zeros(Float64, length(curr_val))
    for term_ix in 1:(curr_stage - 1)
        subsum += (a_coeffs[curr_stage - 1][term_ix] * prev_kvals[term_ix])
    end
    subsum *= curr_step

    return f(t + с * curr_step, curr_val + subsum)
end  # get_kval_rk

"""
	make_step_rk(f::Function, t::Float64, prev_val::Vector{Float64}, curr_step::Float64, 
		max_stage::Integer, c_coeffs::Vector{Float64}, b_coeffs::Vector{Float64}, 
		a_coeffs::Vector{Tuple{Float64}})

Make one step of Runge - Kutta method. Return next value.

# Arguments
- `f::Function`: right part of the system of ODEs;
- `t::Float64: current value of the argument;
- `prev_val::Vector{Float64}`: previous value;
- `curr_step::Float64`: current step value;
- `max_stage::Integer`: total number of stages;
- `c_coeffs::Vector{Float64}`: tuple with ``c`` coefficients of the method;
- `b_coeffs::Vector{Float64}`: tuple with ``b`` coefficients of the method;
- `a_coeffs::Vector{Vector{Float64}}`: ``a`` coefficients of the method.
"""
function make_step_rk(f::Function, t::Float64, prev_val::Vector{Float64}, curr_step::Float64,
	max_stage::Integer, c_coeffs::Vector{Float64}, b_coeffs::Vector{Float64}, 
	a_coeffs::Vector{Vector{Float64}}
)
    # Get k values
    kvals = Vector{Vector{Float64}}()
    for stage in 1:max_stage
        k = get_kval_rk(f, stage, t, curr_step, prev_val, kvals, c_coeffs, a_coeffs)
        push!(kvals, k)
    end

    subsum = zeros(Float64, length(prev_val))
    for term_ix in 1:max_stage
        subsum += (b_coeffs[term_ix] * kvals[term_ix])
    end
    subsum *= curr_step

    return prev_val + subsum
end  # make_step_rk

"""
    maxnorm(vec::Vector{Float64})

Computes 'max' norm of given vector `vec`. Max norm can be calculated as follows:
``\\underset{i = 1, \\dots n}{\\text{max}} |v_i|``. Return norm of `vec` as one real
number.

# Arguments
`vec::Vector{Float64}`: vector which norm should be calculated.
"""
function maxnorm(vec::Vector{Float64})
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
    export_monitor_info(local_err::Vector, global_err::Vector, steps::Vector, 
        vals::Vector)

Export vectors with monitoring information to external files.

# Arguments
`local_err::Vector`: local errors;
`global_err::Vector`: global errors;
`steps::Vector`: method steps;
`vals::Vector`: function values.
"""
function export_monitor_info(local_err::Vector, global_err::Vector, steps::Vector, 
    vals::Vector
)
    monitorpath = "monitoring"
    mkpath(monitorpath)

    open(monitorpath * "/local_err.dat", "w") do outfile
        writedlm(outfile, local_err)
    end

    open(monitorpath * "/global_err.dat", "w") do outfile
        writedlm(outfile, global_err)
    end

    open(monitorpath * "/steps.dat", "w") do outfile
        writedlm(outfile, steps)
    end

    open(monitorpath * "/vals.dat", "w") do outfile
        writedlm(outfile, vals)
    end
end

"""
    solve_rk(f::Function, t_limits::Tuple{Float64, Float64}, initial_step::Float64, 
	    initial_val::Vector{Float64}, tol::Float64, max_stage::Integer, 
        c_coeffs::Vector{Float64}, b_coeffs::Vector{Float64}, 
        a_coeffs::Vector{Vector{Float64}}, fac::Float64, facmin::Float64, facmax::Float64,
        enable_monitoring::Bool = false)

Solve system of ODEs using Runge - Kutta method. System of ODE can be represented as
follows: ``y' = f(t, y), y(t_0) = y_0`` where y is a vector. To achieve more accurancy
this method uses Richardson's extrapolation to control step during the computational 
process.

Return several values packed in the tuple: 
`(<solution>, <total steps count>, <rejected steps count>, <computation time>)`
Solution is represented as a dictionary in following format: `<t>: <y(t)>` where `t` is 
Float64 and `y(t)` is Vector{Float64}.

Function `f` takes ``t`` and vector ``y`` and returns vector of real values which length
is equal to the length of ``y``.

``c`` and ``b`` coefficients are stored as vectors where each element of vector denotes 
coefficient for appropriate stage. ``a`` coefficients are represented as nested vectors. 
Each internal vector contains coeffcients of one method's stage. For example for 4-staged 
method we consider following structure of ``a`` coefficients:
[[a21], [a31, a32], [a41, a42, a43]]. Length of vectors with coefficients should correspond 
to `max_stage` value.

If `enable_monitoring::Bool` is set to `true` then method outputs some internal paraments to
a file. Parametrs are exported as sets of pairs. First value in each pair is the current 
value of argument. Second value in pair can be:

- local error;
- global error;
- current step;
- current function value.

# Arguments
- `f::Function`: right part of the system of ODEs;
- `t_limits::Tuple{Float64, Float64}`: system argument limits;
- `initial_step::Float64`: initial step;
- `initial_val::Vector{Float64}`: initial values (``y_0``);
- `tol::Float64`: target tolerance;
- `max_stage::Integer`: total number of stages;
- `c_coeffs::Vector{Float64}`: tuple with``c`` coefficients of the method;
- `b_coeffs::Vector{Float64}`: tuple with ``b`` coefficients of the method;
- `a_coeffs::Vector{Vector{Float64}}`: ``a`` coefficients of the method;
- `fac::Float64`: garant factor;
- `facmin::Float64`: max coefficient for step reduce;
- `facmax::Float64`: max coefficient for step induce;
- `enable_monitoring::Bool`: true if algorithm monitoring should be enabled.
"""
function solve_rk(f::Function, t_limits::Tuple{Float64, Float64}, initial_step::Float64, 
	initial_val::Vector{Float64}, tol::Float64, max_stage::Integer, 
    c_coeffs::Vector{Float64}, b_coeffs::Vector{Float64}, a_coeffs::Vector{Vector{Float64}},
    fac::Float64, facmin::Float64, facmax::Float64, enable_monitoring::Bool = false
)
    curr_val = initial_val
    curr_step = initial_step
    t = t_limits[1]  # Current argument value

    solution = Dict{Float64, Vector{Float64}}()

    solution[t] = curr_val

    total_steps = 0
    rejected_steps = 0

    # Monitoring arrays
    local_err = []
    global_err = []
    steps = []
    vals = []

    if enable_monitoring
        # Initial values
        push!(local_err, [t, 0])
        push!(global_err, [t, 0])
        push!(steps, [t, curr_step])
        push!(vals, [t, curr_val])
    end

    total_time = @elapsed begin
    while t < t_limits[2]
        near_end = false
        if t + curr_step * 2 > t_limits[2]
            near_end = true
        end

        if near_end
            curr_step = (t_limits[2] - t) / 2
        end

        y1 = 0.
        y2 = 0.
        err = 0.
        glob_err_y1 = 0.
        glob_err_y2 = 0.

        init_facmax = facmax

        while true
            total_steps += 1

            # Make 2 usual steps
            y1 = make_step_rk(f, t, curr_val, curr_step, max_stage, c_coeffs, b_coeffs, 
                a_coeffs)
            y2 = make_step_rk(f, t + curr_step, y1, curr_step, max_stage, c_coeffs, 
                b_coeffs, a_coeffs)

            # Make 1 double step
            w = make_step_rk(f, t, curr_val, curr_step * 2, max_stage, c_coeffs, b_coeffs, 
                a_coeffs)

            err = 1 / (2^max_stage - 1) * maxnorm(y2 - w)

            if enable_monitoring
                glob_err_y1 = maxnorm(y1 - exact_sol(t + curr_step))
                glob_err_y2 = maxnorm(y2 - exact_sol(t + curr_step * 2))
            end

            # If we are near the end of interval then we don't need to perform Richardson's 
            # extrapolation
            if near_end
                solution[t + curr_step] = y1
                solution[t + curr_step * 2] = y2
                break
            end

            new_step = curr_step * min(facmax, max(facmin, fac * 
                (tol / err)^(1 / (max_stage))))

            if err <= tol
                t += (curr_step * 2)
                curr_val = y2
                curr_step = new_step
                break
            else
                rejected_steps += 1
                curr_step = new_step
                facmax = 1
            end
        end

        curr_t1 = 0
        curr_t2 = 0
        if near_end
            curr_t1 = t + curr_step
            curr_t2 = t + curr_step * 2
        else
            curr_t1 = t - curr_step
            curr_t2 = t
        end

        if enable_monitoring
            push!(local_err, [curr_t1, err])
            push!(local_err, [curr_t2, err])
            push!(global_err, [curr_t1, glob_err_y1])
            push!(global_err, [curr_t2, glob_err_y2])
            push!(steps, [curr_t1, curr_step])
            push!(steps, [curr_t2, curr_step])
            push!(vals, [curr_t1, y1])
            push!(vals, [curr_t2, y2])
        end

        if near_end
            break
        end

        facmax = init_facmax

        solution[t - curr_step] = y1
        solution[t] = y2
    end
    end

    if enable_monitoring
        export_monitor_info(local_err, global_err, steps, vals)
    end

    return solution, total_steps, rejected_steps, total_time
end  # solve_rk

"""
	solve_nested_rk(f::Function, t_limits::Tuple{Float64, Float64}, initial_step::Float64,
		initial_val::Vector{Float64}, tol::Float64, max_stage::Integer, 
        c_coeffs::Vector{Float64}, b_coeffs::Vector{Float64}, bhat_coeffs::Vector{Float64},
        a_coeffs::Vector{Vector{Float64}}, enable_monitoring::Bool = false)

Solve system of ODEs using nested Runge - Kutta method. System of ODE can be represented
as follows: ``y' = f(t, y), y(t_0) = y_0`` where y is a vector. Nested Runge - Kutta 
method makes 2 calculations at one step. Then it adapts step according to difference
between these 2 values.

Function `f` takes ``t`` and vector ``y`` and returns vector of real values which length
is equal to the length of ``y``.

``c``, ``b`` and ``\\hat{b}`` coefficients are stored as usual vectors where each element
of vector denotes coefficient for appropriate stage. ``a`` coefficients are represented 
as nested vectors. Each internal vector contains coeffcients of one method's stage. 
For example for 4-staged method we consider following structure of ``a`` coefficients:
[[a21], [a31, a32], [a41, a42, a43]]. Length of tuples with coefficients should 
correspond to `max_stage` value.

If `enable_monitoring::Bool` is set to `true` then method outputs some internal paraments to
a file. Parametrs are exported as sets of pairs. First value in each pair is the current 
value of argument. Second value in pair can be:

- local error;
- global error;
- current step;
- current function value.

# Arguments
- `f::Function`: right part of the system of ODEs;
- `t_limits::Tuple{Float64, Float64}`: system argument limits;
- `initial_step::Float64`: initial step;
- `initial_val::Vector{Float64}`: initial values (``y_0``);
- `tol::Float64`: target tolerance;
- `max_stage::Integer`: total number of stages;
- `c_coeffs::Vector{Float64}`: tuple with``c`` coefficients of the method;
- `b_coeffs::Vector{Float64}`: tuple with ``b`` coefficients of the method;
- `bhat_coeffs::Vector{Float64}`: tuple with ``\\hat{b}`` coefficients of the method;
- `a_coeffs::Vector{Vector{Float64}}`: ``a`` coefficients of the method;
- `enable_monitoring::Bool`: true if algorithm monitoring should be enabled.
"""
function solve_nested_rk(f::Function, t_limits::Tuple{Float64, Float64}, 
    initial_step::Float64, initial_val::Vector{Float64}, tol::Float64, max_stage::Integer, 
    c_coeffs::Vector{Float64}, b_coeffs::Vector{Float64}, bhat_coeffs::Vector{Float64},  
    a_coeffs::Vector{Vector{Float64}}, enable_monitoring::Bool = false
)
    curr_val = initial_val
    curr_step = initial_step
    t = t_limits[1]

    solution = Dict{Float64, Vector{Float64}}()

    solution[t] = curr_val

    # Monitoring arrays
    local_err = []
    global_err = []
    steps = []
    vals = []

    if enable_monitoring
        # Initial values
        push!(local_err, [t, 0])
        push!(global_err, [t, 0])
        push!(steps, [t, curr_step])
        push!(vals, [t, curr_val])
    end

    while t < t_limits[2]
        near_end = false
        if t + curr_step > t_limits[2]
            near_end = true
        end

        if near_end
            curr_step = t_limits[2] - t
        end

        y1 = 0.
        err = 0.
        glob_err_y1 = 0.
        fac = 0.9
        facmin = 0.5
        facmax = 3

        while true
            # Main step
            y1 = make_step_rk(f, t, curr_val, curr_step, max_stage, c_coeffs, b_coeffs, 
                a_coeffs)


            # Step for accuracy control
            y1hat = make_step_rk(f, t, curr_val, curr_step, max_stage, c_coeffs, 
                bhat_coeffs, a_coeffs)

            err = maxnorm(y1 - y1hat)

            if enable_monitoring
                glob_err_y1 = maxnorm(y1 - exact_sol(t + curr_step))
            end

            if near_end
                solution[t + curr_step] = y1
                break
            end

            new_step = curr_step * min(facmax, max(facmin, fac * 
                (tol / err)^(1 / (max_stage))))

            if err <= tol
                t += curr_step
                curr_val = y1
                curr_step = new_step
                break
            else
                curr_step = new_step
                facmax = 1
            end
        end

        curr_t = 0
        if near_end
            curr_t = t + curr_step
        else
            curr_t = t
        end

        if enable_monitoring
            push!(local_err, [curr_t, err])
            push!(global_err, [curr_t, glob_err_y1])
            push!(steps, [curr_t, curr_step])
            push!(vals, [curr_t, y1])
        end

        if near_end
            break
        end

        solution[t] = y1
    end

    if enable_monitoring
        export_monitor_info(local_err, global_err, steps, vals)
    end

    return solution
end  # solve_nested_rk
