
N = [2.0, 3.0];

%Lc= 3.513830719;
for n = N
[ U0,exact, rel_error]=general_cases(n, 0.01,0.1, 100)
n
U0
exact
rel_error
end



function[ U0,exact, rel_error] = general_cases(N,convergence_tolerance_step, convergence_tolerance_abs, max_iterations)
    h= 1/N;
    Xn= [h:h:1-h];
    lambda = 3.513830719;
      denominator=h^2;
    %denominator=2*log(cosh(h));

    %last_tolerance=0;
    
   % N=3;

   %Calculate the exact value
   %theta = 4.7986906881905903435512277458101
   syms x
   theta= vpasolve(0==(sqrt(2*lambda)*cosh(x/4))-x , x);
   exact = -2*log(cosh((Xn-1/2)*theta/2)/cosh(theta/4))';

   %U0=exact

   %return

    derivatives=zeros([0 N+1]);
    


     for j= 0:N-1
        row = [zeros([1 j]),1, -2, 1, zeros([1 N-j-2])];
        derivatives= [derivatives;row(1:N+1)]; 
     end
     derivatives = derivatives(1:N-1, 2:N)./denominator;

         U0= ones([1 N-1])'; %[1.0311, 1.0311]'; %
f = derivatives*U0 + lambda*exp(U0);
  
    for iterations = 1:max_iterations
    %while true
            %second_derivative= derivatives*U0
          
          
    %return
    %tmp = ones([1 N-1]);

    exponentials= lambda*exp(U0);

    exponential_j = diag(exponentials);

    jacobian= derivatives + exponential_j;

    previous_U0 = U0;

    U0 = U0 - jacobian\f;

    previous_f=f;
  f = derivatives*U0 + lambda*exp(U0);
 
    abs_error = norm(exact-U0);
    rel_error = abs_error/norm(exact);
        
%return
        if ~isfinite(U0)
            U0=nan;
            return;
        end
        %if abs(norm(f)<convergence_tolerance_step
         %   changeF=norm(f-previous_f)
        %end
        if   norm(U0-previous_U0)<convergence_tolerance_step && norm(f)<convergence_tolerance_abs%rel_error <convergence_tolerance_step && norm(f)<convergence_tolerance_abs %abs(norm(U0-previous_U0))<convergence_tolerance_step &&%
            f;
            U0;
            return;
        end
    end
    %U0 = nan;
    return;
end
