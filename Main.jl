using DelimitedFiles

include("RungeKutta.jl")
include("Task.jl")

"""
Possible arguments which are used as iteratables when printing results to files.
"""
@enum ResultArgument begin
    InitStep
    Fac
    FacMax
    FacMin
end

"""
    test_rk()

Take problem and initial condiditions from `Task.jl` and solve it with Runge - Kutta method.
This method also validates the answer.
"""
function test_rk()
    t_limits = (0.1, 4.1)
    initial_step = 0.05
    initial_val = exact_sol(t_limits[1])
    tol = 1e-4
    max_stage = 4
    coeffs = coeffs_rk()
    fac = 0.9
    facmin = 0.5
    facmax = 3.

    @info("Solving ODE with Runge - Kutta method")
    (num_sol, total_steps, rejected_steps, total_time) = solve_rk(f, t_limits, initial_step, 
        initial_val, tol, max_stage, coeffs[1], coeffs[2], coeffs[3], fac, facmin, facmax,
        true)

    @info ("Validating result")
    err_vec = Vector{Real}()
    for (t, y) in num_sol
        exact_val = exact_sol(t)
        append!(err_vec, maxnorm(exact_val - y))
    end

    max_err = 0
    for err in err_vec
        if err > max_err
            max_err = err
        end
    end

    println("Max error = ", max_err)
    println("Total steps = ", total_steps)
    println("Rejected steps = ", rejected_steps)
    println("Total time = ", total_time)
end  # test_rk

"""
    test_dp()

Take problem and initial condiditions from `Task.jl` and solve it with Dormand - Prince 
method. This method also validates the answer.
"""
function test_dp()
    t_limits = (0.1, 4.1)
    initial_step = 0.05
    initial_val = exact_sol(t_limits[1])
    tol = 1e-4
    max_stage = 7
    coeffs = coeffs_dp()
    fac = 0.9
    facmin = 0.5
    facmax = 3.

    @info("Solving ODE with Dormand - Prince method")
    (num_sol, total_steps, rejected_steps, total_time) = solve_nested_rk(f, t_limits, 
        initial_step, initial_val, tol, max_stage, coeffs[1], coeffs[2], coeffs[3], 
        coeffs[4], fac, facmin, facmax, true)

    @info ("Validating result")
    err_vec = Vector{Real}()
    for (t, y) in num_sol
        exact_val = exact_sol(t)
        append!(err_vec, maxnorm(exact_val - y))
    end

    max_err = 0
    for err in err_vec
        if err > max_err
            max_err = err
        end
    end

    println("Max error = ", max_err)
    println("Total steps = ", total_steps)
    println("Rejected steps = ", rejected_steps)
    println("Total time = ", total_time)
end  # test_dp

"""
    output_by_argument_rk()

Output Runge - Kutta calculation results as a tables. Argument of each table is `arg`. Table
values: total count of steps, count of rejected steps, total calculation time, global error.

# Arguments
`arg::ResultArgument`: parameter which should be the argument of tables.
"""
function output_by_argument_rk(arg::ResultArgument)
    t_limits = (0.1, 4.1)
    initial_steps = [step for step in 0.01:0.001:0.2]
    initial_val = exact_sol(t_limits[1])
    tol = 1e-4
    max_stage = 4
    coeffs = coeffs_rk()
    facs = [fac for fac in 0.8:0.001:0.9]
    facmins = [facmin for facmin in 0.2:0.001:0.66]
    facmaxes = [facmax for facmax in 1.5:0.01:5]

    respath = "results"
    mkpath(respath)

    total_steps_data = []
    rejected_steps_data = []
    total_time_data = []
    global_err_data = []

    param_vec = []
    file_prefix = ""
    if arg == InitStep
        file_prefix = "initstep"
        param_vec = initial_steps
    elseif arg == Fac
        file_prefix = "fac"
        param_vec = facs
    elseif arg == FacMax
        file_prefix = "facmax"
        param_vec = facmaxes
    elseif arg == FacMin
        file_prefix = "facmin"
        param_vec = facmins
    end

    for param in param_vec
        res_tuple = Tuple

        if arg == InitStep
            res_tuple = solve_rk(f, t_limits, param, initial_val, tol, max_stage, 
                coeffs[1], coeffs[2], coeffs[3], facs[1], facmins[1], facmaxes[1])
        elseif arg == Fac
            res_tuple = solve_rk(f, t_limits, initial_steps[1], initial_val, tol, max_stage, 
                coeffs[1], coeffs[2], coeffs[3], param, facmins[1], facmaxes[1])
        elseif arg == FacMax
            res_tuple = solve_rk(f, t_limits, initial_steps[1], initial_val, tol, max_stage, 
                coeffs[1], coeffs[2], coeffs[3], facs[1], facmins[1], param)
        elseif arg == FacMin
            res_tuple = solve_rk(f, t_limits, initial_steps[1], initial_val, tol, max_stage, 
                coeffs[1], coeffs[2], coeffs[3], facs[1], param, facmaxes[1])
        end

        num_sol = res_tuple[1]
        total_steps = res_tuple[2]
        rejected_steps = res_tuple[3]
        total_time = res_tuple[4]

        err_vec = Vector{Real}()
        for (t, y) in num_sol
            exact_val = exact_sol(t)
            append!(err_vec, maxnorm(exact_val - y))
        end

        push!(total_steps_data, [param, total_steps])
        push!(rejected_steps_data, [param, rejected_steps])
        push!(total_time_data, [param, total_time])
        push!(global_err_data, [param, max(err_vec...)])
    end

    file_name = respath * "/" * file_prefix * "_n.dat"
    open(file_name, "w") do outfile
        writedlm(outfile, total_steps_data)
    end

    file_name = respath * "/" * file_prefix * "_nrej.dat"
    open(file_name, "w") do outfile
        writedlm(outfile, rejected_steps_data)
    end

    file_name = respath * "/" * file_prefix * "_time.dat"
    open(file_name, "w") do outfile
        writedlm(outfile, total_time_data)
    end

    file_name = respath * "/" * file_prefix * "_err.dat"
    open(file_name, "w") do outfile
        writedlm(outfile, global_err_data)
    end
end  # output_by_argument_rk

"""
    output_results_rk()

Perform calculations with Runge - Kutta method and output results to external files.
"""
function output_results_rk()

    @info "Solving for different initial steps..."
    output_by_argument_rk(InitStep)
    
    @info "Solving for different facs..."
    output_by_argument_rk(Fac)

    @info "Solving for different facmaxes..."
    output_by_argument_rk(FacMax)

    @info "Solving for different facmins..."
    output_by_argument_rk(FacMin)
end  # output_results_rk

# test_rk()
test_dp()
# output_results_rk()
