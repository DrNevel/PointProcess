function HESS = eval_hess(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn, eta)

    % init
    HESS=zeros(length(update)); 

        for i=1:length(update)
            incr=update;
            decr=update;
            step=abs(update(i))*(nthroot(eps,3));

            incr(i) = incr(i) + step;
            decr(i) = decr(i) - step;

            GRAD_incr = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, incr, wn, xn, eta);
            GRAD_decr = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, decr, wn, xn, eta);
            
            HESS(i,:)= (GRAD_incr - GRAD_decr ) ./ (2*step);

    %         for j=1:length(update)          
    %             HESS(i,j)= (GRAD_incr(j) - GRAD_decr(j) ) ./ (2*step) ;
    %         end

        end
    end