
N = [2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0, 200.0, 500.0, 1000.0];

%Lc= 3.513830719;
for n = N
[ U0, rel_error]=general_cases(n, 0.01,2)

end


function[ U0, rel_error] = general_cases(N,convergence_tolerance_step, convergence_tolerance_abs , max_iterations)

    h= 1/N; % h is the inverse of N
    Xn= [h:h:1-h]; % Produces vector of Xi values to calulcate exact values used to calculate relative error
    lambda = 3.513830719; % critical lambda value used to calculate exact value and function

      denominator=h^2; % standard difference
    %denominator=2*log(cosh(h)); % Mickens non-standard difference

   %Calculate the exact value

   syms x
   theta= vpasolve(0==(sqrt(2*lambda)*cosh(x/4))-x , x);
   exact = -2*log(cosh((Xn-1/2)*theta/2)/cosh(theta/4))';
% exact is a vector of U(Xi) of i from 1 => N-1


% using centered divided difference method to estimate second derivative in
% system on non linear equations
    derivatives=zeros([0 N+1]);
    


     for j= 0:N-1 % builds matrix for i from 0=>N
        row = [zeros([1 j]),1, -2, 1, zeros([1 N-j-2])];
        derivatives= [derivatives;row(1:N+1)]; 
     end
     derivatives = derivatives(1:N-1, 2:N)./denominator; % crops matrix to imply boundary conditions, U(0)=U(1)=0

%initial guess of system root, Ui, i from 1 to N-1
U0= ones([1 N-1])'; 
% initial evaluation of function at initial guess
f = derivatives*U0 + lambda*exp(U0);
  for iterations = 1:max_iterations
    %while true

          
          

%calculating Jacobian matrix using the linear cofficients to represent constants and the exponentials
%along the diagnol
    exponentials= lambda*exp(U0);
    exponential_j = diag(exponentials);
    jacobian= derivatives + exponential_j;
% recording the previous aproximation
    previous_U0 = U0;
    U0 = U0 - jacobian\f; % solve for next aproximation of the root by using the inverse multiplication backslash functionality of matlab

    previous_f=f; % recording the previous function output
  f = derivatives*U0 + lambda*exp(U0); %calulcating the funtion output at the new aproximation
 

    abs_error = norm(exact-U0); %calulating the total error between the exact value and newest functional value
    rel_error = abs_error/norm(exact); % calculates the relative error
        

        if ~isfinite(U0)
            U0=nan;
            return;
        end

        if  rel_error <convergence_tolerance_step %norm(f)<convergence_tolerance_abs && %
            f;
            U0;
            return;
        end
    end
    %U0 = nan;
    return;
end
