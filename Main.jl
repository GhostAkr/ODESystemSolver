include("RungeKutta.jl")
include("Task.jl")

"""
    test_rk()

Take problem and initial condiditions from `Task.jl` and solve it with Runge - Kutta method.
This method also validates the answer.
"""
function test_rk()
    t_limits = (0.1, 4.1)
    initial_step = 0.05
    initial_val = exact_sol(t_limits[1])
    tol = 10^(-4)
    max_stage = 4
    coeffs = coeffs_rk()

    @info("Solving ODE with Runge - Kutta method")
    num_sol = solve_rk(f, t_limits, initial_step, initial_val, tol, max_stage, coeffs[1], 
        coeffs[2], coeffs[3])

    # Validation
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
end  # test_rk
