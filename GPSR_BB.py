def GPSR_BB(y, A, tau):
    # test for the number of required parameters

    # flag for initial x (can take any values except 0, 1, 2)
    Initial_X_supplied = 3333

    # Set the defaults for the optional parameters
    stopCriterion       = 3
    tolA                = 1e-6
    tolD                = 0.0001
    debias              = 0
    maxiter             = 10000
    maxiter_debias      = 500
    miniter             = 5
    miniter_debias      = 5
    init                = 0
    enforceMonotone     = 1
    alphamin            = 1e-30
    alphamax            = 1e30
    compute_mse         = 0
    AT                  = 0
    verbose             = 0
    continuation        = 0
    cont_steps          = -1
    firstTauFactorGiven = 0
    x                   = np.mat(np.zeros((1000, 1)))


    # Set the defaults for outputs that may not be computed
    debias_start = 0
    x_debias     = None
    mses         = None

    # # Read the optional parameters
    # for key, value in kwargs.items():
    #     key_upper = key.upper()
    #     if key_upper == 'STOPCRITERION':
    #         stopCriterion = value
    #     elif key_upper == 'TOLERANCEA':
    #         tolA = value
    #     elif key_upper == 'TOLERANCED':
    #         tolD = value
    #     elif key_upper == 'DEBIAS':
    #         debias = value
    #     elif key_upper == 'MAXITERA':
    #         maxiter = value
    #     elif key_upper == 'MAXITERD':
    #         maxiter_debias = value
    #     elif key_upper == 'MINITERA':
    #         miniter = value
    #     elif key_upper == 'MINITERD':
    #         miniter_debias = value
    #     elif key_upper == 'INITIALIZATION':
    #         if np.prod(np.array(value).shape) > 1:
    #             # initial x supplied as an array
    #             init = Initial_X_supplied
    #             x = np.array(value)
    #         else:
    #             init = value
    #     elif key_upper == 'MONOTONE':
    #         enforceMonotone = value
    #     elif key_upper == 'CONTINUATION':
    #         continuation = value
    #     elif key_upper == 'CONTINUATIONSTEPS':
    #         cont_steps = value
    #     elif key_upper == 'FIRSTTAUFACTOR':
    #         firstTauFactor = value
    #         firstTauFactorGiven = 1
    #     elif key_upper == 'TRUE_X':
    #         compute_mse = 1
    #         true = np.array(value)
    #     elif key_upper == 'ALPHAMIN':
    #         alphamin = value
    #     elif key_upper == 'ALPHAMAX':
    #         alphamax = value
    #     elif key_upper == 'VERBOSE':
    #         verbose = value
    #     else:
    #         # Hmmm, something wrong with the parameter string
    #         raise ValueError(f'Unrecognized option: {key}')
        
    # Check stopping criterion
    if np.sum(np.isin(stopCriterion, [0, 1, 2, 3, 4, 5])) == 0:
        raise ValueError('Unknown stopping criterion')

    


    # Precompute A'*y
    


    # Initialization
    if init == 0:
        x = np.matmul(np.transpose(A),np.zeros_like(y))
        # x = np.mat(np.concatenate([np.transpose(A), np.zeros_like(y).reshape(-1, 1)], axis=1))

    
    # Check if tau is an array
    if np.prod(np.array(tau).shape) > 1:
        try:
            dummy = np.matmul(x,tau)
        except:
            raise ValueError('Parameter tau has wrong dimensions; it should be scalar or size(x)')

    # # Check the size of the true x if given
    # if compute_mse and true.shape != x.shape:
    #     raise ValueError('Initial x has incompatible size')

    # Check the value of tau for scalar tau
    if np.prod(np.array(tau).shape) == 1:
        aux = np.matmul(np.transpose(A),y)
        max_tau = np.max(np.abs(aux))
        if tau >= max_tau:
            x = np.zeros_like(aux)
            if debias:
                x_debias = x
            objective[0] = 0.5 * np.vdot(y,y)
            times[0] = 0
            # if compute_mse:
            #     mses[0] = np.sum(true.ravel() ** 2)
            # return x, x_debias, objective, times, debias_start, mses, taus

    # Initialize u and v
    u = np.maximum(x, 0)
    v = -np.minimum(x, 0)

    # Define the indicator vector or matrix of nonzeros in x
    # nz_x = (x != 0.0)
    nz_x     =  ~np.isin(x,0.0)
    num_nz_x = np.sum(nz_x)

    
    # Start the clock
    t0 = time.time()

    # Store given tau for continuation procedure
    final_tau = tau

    # Store given stopping criterion and threshold for continuation procedure
    final_stopCriterion = stopCriterion
    final_tolA = tolA

    # Set continuation factors
    if continuation and (cont_steps > 1):
        if np.prod(np.array(tau).shape) == 1:
            if (firstTauFactorGiven == 0) or (firstTauFactor * tau >= max_tau):
                firstTauFactor = 0.5 * max_tau / tau
                if verbose:
                    print('\n setting parameter FirstTauFactor\n')
        cont_factors = 10 ** np.linspace(np.log10(firstTauFactor),0,num=cont_steps)
    if not continuation:
        cont_factors = np.array([1])
        cont_steps = 1

    iter = 1
    if compute_mse:
        mses[iter] = np.sum((x - true) ** 2)


       
    keep_continuation = 1
    cont_loop = 1
    iter = 1
    taus = []

    # loop for continuation
    while keep_continuation:
        # Compute and store initial value of the objective function
        resid = y - np.matmul(A,x)

        if cont_steps == -1:
            gradq = np.matmul(np.transpose(A),resid)
            tau = max(final_tau, 0.2 * np.max(np.abs(gradq)))
            if tau == final_tau:
                stopCriterion = final_stopCriterion
                tolA = final_tolA
                keep_continuation = 0
            else:
                stopCriterion = 1
                tolA = 1e-5
        else:
            tau = final_tau * cont_factors[cont_loop - 1]
            if cont_loop == cont_steps:
                stopCriterion = final_stopCriterion
                tolA = final_tolA
                keep_continuation = 0
            else:
                stopCriterion = 1
                tolA = 1e-5

        taus.append(tau)

        if verbose:
            print(f'\nSetting tau = {tau}')

        # if in the first continuation iteration, compute and store
        # initial value of the objective function
        if cont_loop == 1:
            alpha = 1.0
            f = 0.5 * np.vdot(resid, resid) + np.sum(tau * u) + np.sum(tau * v)

            objective = f
            if compute_mse:
                tempp = x - true
                mses = np.vdot(tempp,tempp)
            if verbose:
                print(f'Initial obj={f}, alpha={alpha}, nonzeros={num_nz_x}')

        # Compute the initial gradient and the useful
        # quantity resid_base
        resid_base = y - resid

        # control variable for the outer loop and iteration counter
        keep_going = 1

        if verbose:
            print(f'\nInitial obj={f}, nonzeros={num_nz_x}')

        while keep_going:
            # compute gradient
            temp = np.matmul(np.transpose(A),resid_base)

            term  =  temp - np.matmul(np.transpose(A),y)
            gradu =  term + tau
            gradv = -term + tau

            # projection and computation of the search direction vector
            du = np.maximum(u - alpha * gradu, 0.0) - u
            dv = np.maximum(v - alpha * gradv, 0.0) - v
            dx = du - dv
            old_u = u
            old_v = v

            # calculate useful matrix-vector product involving dx
            auv = np.matmul(A,dx)
            dGd = np.vdot(auv,auv)

            if enforceMonotone == 1:
                # monotone variant: calculate the minimizer along the direction (du, dv)
                lambda0 = - (np.vdot(gradu,du) + np.vdot(gradv, dv)) / (np.finfo(float).eps + dGd)
                if lambda0 < 0:
                    print(f'ERROR: lambda0 = {lambda0} negative. Quit')
                    return
                lambda_val = min(lambda0, 1)
            else:
                # nonmonotone variant: choose lambda=1
                lambda_val = 1

            u = old_u + lambda_val * du
            v = old_v + lambda_val * dv
            uvmin = np.minimum(u, v)
            u = u - uvmin
            v = v - uvmin
            x = u - v


           
            # calculate nonzero pattern and number of nonzeros (do this *always*)
            nz_x_prev = np.copy(nz_x)
            # nz_x = (x != 0.0)
            nz_x     =  ~np.isin(x,0.0)
            num_nz_x = np.sum(nz_x)

            # update residual and function
            resid = y - resid_base - lambda_val * auv
            prev_f = f
            f = 0.5 * np.vdot(resid, resid) + np.sum(tau * u) + np.sum(tau * v)


            # compute new alpha
            dd = np.vdot(du, du) + np.vdot(dv, dv)
            if dGd <= 0:
                # something wrong if we get to here
                print(f'dGd={dGd}, nonpositive curvature detected')
                alpha = alphamax
            else:
                alpha = min(alphamax, max(alphamin, dd / dGd))
            resid_base = resid_base + lambda_val * auv

            # print out stuff
            if verbose:
                print(f'It={iter}, obj={f}, alpha={alpha}, nz={num_nz_x}')

            # update iteration counts, store results, and times
            iter += 1
            objective= f
            t2 = time.time() - t0

            if compute_mse:
                err = true - x
                mses.append(np.vdot(err, err))

            # Stopping criteria
            if stopCriterion == 0:
                # compute the stopping criterion based on the change
                # of the number of non-zero components of the estimate
                num_changes_active = np.sum(nz_x != nz_x_prev)
                if num_nz_x >= 1:
                    criterionActiveSet = num_changes_active
                else:
                    criterionActiveSet = tolA / 2
                keep_going = criterionActiveSet > tolA
                if verbose:
                    print(f'Delta n-zeros = {num_changes_active} (target = {tolA})')
            elif stopCriterion == 1:
                # compute the stopping criterion based on the relative
                # variation of the objective function.
                criterionObjective = abs(f - prev_f) / prev_f
                keep_going = criterionObjective > tolA
                if verbose:
                    print(f'Delta obj. = {criterionObjective} (target = {tolA})')
            elif stopCriterion == 2:
                # stopping criterion based on relative norm of step taken
                delta_x_criterion = np.linalg.norm(dx) / np.linalg.norm(x)
                keep_going = delta_x_criterion > tolA
                if verbose:
                    print(f'Norm(delta x)/norm(x) = {delta_x_criterion} (target = {tolA})')
            elif stopCriterion == 3:
                # compute the "LCP" stopping criterion
                w = np.concatenate([np.minimum(gradu, old_u), np.minimum(gradv, old_v)])
                criterionLCP = np.linalg.norm(w, np.inf)
                criterionLCP = criterionLCP / max([1.0e-6, np.linalg.norm(old_u, np.inf), np.linalg.norm(old_v, np.inf)])
                keep_going = criterionLCP > tolA
                if verbose:
                    print(f'LCP = {criterionLCP} (target = {tolA})')
            elif stopCriterion == 4:
                # continue if not yet reached target value tolA
                keep_going = f > tolA
                if verbose:
                    print(f'Objective = {f} (target = {tolA})')
            elif stopCriterion == 5:
                # stopping criterion based on relative norm of step taken
                delta_x_criterion = np.sqrt(dd) / np.sqrt(np.dot(x.ravel(), x.ravel()))
                keep_going = delta_x_criterion > tolA
                if verbose:
                    print(f'Norm(delta x)/norm(x) = {delta_x_criterion} (target = {tolA})')
            else:
                raise ValueError('Unknown stopping criterion')

            # take no less than miniter...
            if iter <= miniter:
                keep_going = True
            elif iter > maxiter:  # and no more than maxiter iterations
                keep_going = False

        # increment continuation loop counter
        cont_loop += 1

    if verbose:
        print('\nFinished the main algorithm!\nResults:')
        print('||A x - y ||_2^2 = {:.3e}'.format(np.vdot(resid, resid)))
        print('||x||_1 = {:.3e}'.format(np.sum(np.abs(x))))
        print('Objective function = {:.3e}'.format(f))
        # nz_x = (x != 0.0)
        nz_x     =  ~np.isin(x,0.0)
        num_nz_x = np.sum(nz_x)
        print('Number of non-zero components = {}'.format(num_nz_x))
        print('CPU time so far = {:.3e}'.format(t2))
        print('\n')

    # If the 'Debias' option is set to 1, try to remove the bias from the l1 penalty
    # if debias and np.sum(x != 0) != 0:
    if debias and (np.sum(~np.isin(x,0.0)) != 0):
        if num_nz_x > len(y):
            if verbose:
                print('\nDebiasing requested, but not performed')
                print('There are too many nonzeros in x\n')
                print('nonzeros in x: {}, length of y: {}'.format(num_nz_x, len(y)))
        elif num_nz_x == 0:
            if verbose:
                print('\nDebiasing requested, but not performed')
                print('x has no nonzeros\n')
        else:
            if verbose:
                print('\nStarting the debiasing phase...\n')

            x_debias = np.copy(x)
            # zeroind = (x_debias != 0)
            zeroind = ~np.isin(x_debias,0.0)
            cont_debias_cg = 1
            debias_start = iter

            # calculate initial residual
            resid = np.matmul(A,x_debias)
            resid = resid - y
            resid_prev = np.finfo(float).eps * np.ones(resid.shape)

            rvec = np.matmul(np.transpose(A),resid)
            rvec = np.matmul(rvec,zeroind)
            rTr_cg = np.vdot(rvec, rvec)

            tol_debias = tolD * rTr_cg

            pvec = -rvec

            # main loop
            while cont_debias_cg:
                RWpvec = np.matmul(A,pvec)
                Apvec = np.matmul(np.transpose(A),RWpvec)
                Apvec = np.matmul(Apvec,zeroind)

                alpha_cg = rTr_cg / (np.vdot(pvec, Apvec))

                x_debias = x_debias + alpha_cg * pvec
                resid = resid + alpha_cg * RWpvec
                rvec = rvec + alpha_cg * Apvec

                rTr_cg_plus = np.vdot(rvec, rvec)
                beta_cg = rTr_cg_plus / rTr_cg
                pvec = -rvec + beta_cg * pvec
                rTr_cg = rTr_cg_plus

                iter += 1

                objective.append(0.5 * np.vdot(resid, resid) +
                                 np.sum(np.matmul(tau,np.abs(x_debias.ravel()))))
                times.append(time.time() - t0)

                if compute_mse:
                    err = true - x_debias
                    mses.append(np.vdot(err, err))

                if verbose:
                    print(' Iter = {:5d}, debias resid = {:.8e}, convergence = {:.3e}'.format(
                        iter, np.vdot(resid, resid), rTr_cg / tol_debias))

                cont_debias_cg = (iter - debias_start <= miniter_debias) or \
                                 ((rTr_cg > tol_debias) and (iter - debias_start <= maxiter_debias))

            if verbose:
                print('\nFinished the debiasing phase!\nResults:')
                print('||A x - y ||_2^2 = {:.3e}'.format(np.vdot(resid, resid)))
                print('||x||_1 = {:.3e}'.format(np.sum(np.abs(x_debias))))
                print('Objective function = {:.3e}'.format(objective[-1]))
                nz = (x_debias != 0.0)
                print('Number of non-zero components = {}'.format(np.sum(nz)))
                print('CPU time so far = {:.3e}'.format(times[-1]))
                print('\n')

    if compute_mse:
        mses = mses / len(true)
    
    return x
