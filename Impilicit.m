
landa = input("Thermal Conductivity: ");
rho = input("Density: ");
Cp = input("Specific Heat: ");
alpha=landa/(rho*Cp);

L = input("Length of Beam: ");
Time = input("Total Time: ");
deltax = input("DeltaX: ");
deltat = input("Deltat: ");
tol = input("Tolerance: ");
max_iter = input("Max iteratio: ");

N = L/deltax;
M = Time/deltat;

T = ones(1,N)*30;
T(1) = 100;
T(N) = 20;


C = zeros(N,N);
B = zeros (N,1); 

C(1,1) = 1;
C(N,N) = 1; 

for i=2:N-1

    C(i,i-1) = -alpha*(deltat/deltax^2); 
    C(i,i) = 1+2*alpha*(deltat/deltax^2);
    C(i,i+1) =  -alpha*(deltat/deltax^2); 

end

T1 = zeros (M,N);
T1 (1,:) = T; 
x0 = ones(N,1)*30;

function x = gauss_seidel(C, B, x0, tol, max_iter)

    x = x0;
    N = length(B);

    for iter = 1:max_iter
        x0 = x;
        for i = 1:N
           
            x(i) = (B(i) - C(i, 1:i-1) * x(1:i-1) - C(i, i+1:end) * x(i+1:end)) / C(i, i);
        end

    
        if norm(x - x0, inf) < tol
            break;
        end
    end
end


for t = 2:M
    
    B(:) = T;
    B(1) = 100;  
    B(N) = 20;  
    
   
    T_new = gauss_seidel(C, B, T', tol, max_iter);
    
    
    T = T_new';
    T1(t, :) = T;
end


figure; 
plot(0:deltat:Time-deltat,T1(:,N/2)); 
xlabel('t(s)');
ylabel('T (x=0.5m)');
title('Temperature distribution over time at x = 0.5 m');
legend('show');
grid on; 


figure; 
plot(0:deltax:L-deltax,T1(800/deltat,:));
xlabel('x (m)');
ylabel('T (t = 800 s)');
title('Temperature distribution along beam at t = 800 s');
legend('show');
grid on; 

disp ("T(x=0.75,t=600): ");
disp(T1(600/deltat,0.75/deltax));


