function plot_results(n,p,N,junc,t,x_est)
% Based on the type of initial node and end node,
% the number of state variables will be different,
% we defined and used a variable "a" to apply this change

[row_p,col_p] = size(p);
junc_size = max(size(junc));
c=0;

for i = 1:row_p
    if n(p(i,1),2) == 1 && n(p(i,2),2) == 2
        a=1;
    elseif n(p(i,1),2) == 1 && n(p(i,2),2) == 3
        a=0;
    elseif n(p(i,1),2) == 2 && n(p(i,2),2) == 2
        a=2;
    elseif n(p(i,1),2) == 2 && n(p(i,2),2) == 3
        a=1;
    end
figure(i)
subplot(2,1,1)
plot(t,x_est(:,(1+c):((N(i)+1)-a)+c))
str1 = sprintf('Pressure for inner nodes and nodes of pipe %g (except for junctions)', i);
title(str1)
xlabel('t (s)')
ylabel('Pressure (Pa)')
subplot(2,1,2)
plot(t,x_est(:,(N(i)+2-a+c):(2*(N(i)+1)-a)+c))
str2 = sprintf('Mass Flow Rate for inner nodes and nodes of pipe %g', i);
title(str2)
c=c+2*(N(i)+1)-a;
xlabel('t (s)')
ylabel('Mass flow rate (kg/s)')
end

figure(i+1)
str3 = sprintf('Pressure for junctions');
title(str3)
hold on
for j = 1:junc_size
plot(t,x_est(:,c+j))
xlabel('t (s)')
ylabel('Pressure (Pa)')
node_index = junc(j);
str{j} = sprintf('Junction Node Index = %g', node_index);
end
legend(str)

end