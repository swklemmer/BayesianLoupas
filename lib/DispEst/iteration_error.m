function stop = iteration_error(x, optimValues, state)
%ITERATION_ERROR Runs on every iteration to track error versus elapsed time
% x = current solution
% optimValues = optimization data for current iteration

persistent max_iter u_true t elapsed_t err fun_eval

% Save measurements to base
switch state
    case 'init'
        % Begin timer
        tic();
        max_iter = evalin('base', 'opt_param.MaxIterations');
        elapsed_t = zeros(max_iter, 1);
    
        % Begin iteration count
        t = 1;
            
        % Load true displacement
        u_true = evalin('base', 'u_tru');
        err = zeros(max_iter, 3);

        % Begin function evaluation count
        fun_eval = zeros(max_iter, 1);

        fprintf('    INIT\n');

    case 'iter'
        % Track elapsed time and assign to base
        elapsed_t(t) = toc();
        
        % Calculate error metrics 
        err(t, :) = error_metrics(x, u_true, 1);

        % Track number of function evaluations
        fun_eval(t) = optimValues.funccount;
        
        % Increment iteration count
        fprintf('    Iter. %d\n', t);
        if t < max_iter; t = t + 1; end

    case 'done'
        % Assign variables to base
        assignin('base', 'elapsed_t_j', elapsed_t);
        assignin('base', 'err_j', err);
        assignin('base', 'fun_eval_j', fun_eval);
end

stop = false;
end
