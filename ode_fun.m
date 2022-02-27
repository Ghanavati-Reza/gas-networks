function dydt = ode_fun(t,y,N,p,n,dx,sg,d,F,T)
%% Defining ODE equations

A = pi*d^2/4; % m^2
R = 8.314; % J/mol.K
v = 0.1; % time constant will be used in ode equation of junctions

[row_n,col_n] = size(n);
[row_p,col_p] = size(p);

junc_indexes = [];
c_1 = 1;
for i = 1:row_n
    if n(i,2) == 2
        junc_indexes(c_1) = i;
        c_1 = c_1+1;
    end
end

if ~isempty(junc_indexes)
length_junc = max(size(junc_indexes));
juncs_data = [junc_indexes' zeros(length_junc,1)];
end

%% Calculating number of states

c_2 = 0;
for i = 1:row_p
    for j = 1:col_p
        for k = 1:length_junc
            if p(i,j) == junc_indexes(k)
                c_2 = c_2+1;
            end
        end
    end
end

state_number = 2*(sum(N)+row_p)-c_2+length_junc;

%% Calculate C^2 = 1000 * z * R * T / Mw; C unit(m/s)
% for gas in pipeline we have C^2 = p/density
c_3 = 0;
c_4 = 0;
for i = 1:row_p
    if n(p(i,1),2) == 1 && n(p(i,2),2) == 2
        a_1=1;
    elseif n(p(i,1),2) == 1 && n(p(i,2),2) == 3
        a_1=0;
    elseif n(p(i,1),2) == 2 && n(p(i,2),2) == 2
        a_1=2;
    elseif n(p(i,1),2) == 2 && n(p(i,2),2) == 3
        a_1=1;
    end
    % definig AGA compresibility factor
    for j = 1:N(i)+1-a_1
        P = y(j+c_3);
        tc2 = 169 + 314 * sg;
        tc = (273.15 + ((tc2) - 459.67 - 32)/1.8);
        pc2 = 708.75 - 57.5 * sg;
        pc = ((pc2)*1e5/14.5037738);
        Z = 1 + (P / pc) * (0.257 - 0.533 / (T / tc));
        Mw = sg*28.969;
        C2(j+c_4,1) = 1000*Z*R*T/Mw;
    end
    c_3 = c_3+j+(N(i)+1);
    c_4 = c_4+j;
end
for i = 1:length_junc
        P = y(i+c_3);
        juncs_data(i,3) = i+c_3;
        tc2 = 169 + 314 * sg;
        tc = (273.15 + ((tc2) - 459.67 - 32)/1.8);
        pc2 = 708.75 - 57.5 * sg;
        pc = ((pc2)*1e5/14.5037738);
        Z = 1 + (P / pc) * (0.257 - 0.533 / (T / tc));
        Mw = sg*28.969;
        C2(i+c_4,1) = 1000*Z*R*T/Mw;
end
      
%% Define dydt for all states except for junctions pressure   
% I divided equations to 4 different parts based on the type of the
% pipeline: 1) start with supply node and end to junction node 2)start 
% with supply node and end to delivery node 3)start with junction node and 
% end to junction node 4) start with junction node and end to delivery node
dydt = zeros(state_number,1);

c_5 = 0;
c_6 = 0;
for i = 1:row_p
      if n(p(i,1),2) == 1 && n(p(i,2),2) == 2
            a_1 = 1;
          for j = (1+c_5):((N(i)+1)-a_1+c_5)
            if j == 1+c_5 
                dydt(j,1) = 0;                
            elseif j == 2+c_5
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j-a_1)-y(N(i)+j-a_1));
            elseif j == N(i)+c_5
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j-a_1)-y(N(i)+j-a_1));
                for k = 1:length_junc
                        if n(p(i,2),1) == juncs_data(k,1)
                            juncs_data(k,2) = juncs_data(k,2)+y(N(i)+1+j-a_1);
                            P_j = y(juncs_data(k,3));
                        end
                end
            else
              dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j-a_1)-y(N(i)+j-a_1));
            end
          end
         for j = (c_5+(N(i)+2)-a_1):(c_5+(2*(N(i)+1))-a_1)
            if j == N(i)+2-a_1+c_5
              y(j-N(i)-1+a_1)=n(p(i,1),3);
              dydt(j,1) = (-A*(y(j-N(i)+a_1)-y(j-N(i)-1+a_1))/(dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/((y(j-N(i)-1+a_1)))));
            elseif j == N(i)+3-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i)+a_1)-y(j-N(i)-2+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-N(i)-1+a_1))));
            elseif j == 2*(N(i)+1)-1-a_1+c_5
              dydt(j,1) = (-A*(P_j-y(j-(N(i)+2)+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+1)+a_1))));
            elseif j == 2*(N(i)+1)-a_1+c_5
              dydt(j,1) = (-A*(P_j-y(j-(N(i)+2)+a_1))/(dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(P_j)));
            else
              dydt(j,1) = (-A*(y(j-N(i)+a_1)-y(j-(N(i)+2)+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+1)+a_1))));
            end
        end
      elseif n(p(i,1),2) == 1 && n(p(i,2),2)==3
        a_1 = 0;
          for j = (1+c_5):((N(i)+1)-a_1+c_5)
            if j == 1+c_5 
                    dydt(j,1) = 0;                
            elseif j == 2+c_5
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j-a_1)-y(N(i)+j-a_1));
            elseif j == N(i)+c_5
                y(N(i)+2+j-a_1) = n(p(i,2),3);
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j-a_1)-y(N(i)+j-a_1));
            elseif j == N(i)+1+c_5
                dydt(j,1) = (-(C2(j-c_6))/(A*dx(i)))*(y(N(i)+1+j-a_1)-y(N(i)+j-a_1));
            else
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j-a_1)-y(N(i)+j-a_1));
            end
          end
         for j = (c_5+(N(i)+2)-a_1):(c_5+(2*(N(i)+1))-a_1)
            if j == N(i)+2-a_1+c_5
              y(j-N(i)-1+a_1) = n(p(i,1),3);
              dydt(j,1) = (-A*(y(j-N(i)+a_1)-y(j-N(i)-1+a_1))/(dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/((y(j-N(i)-1+a_1)))));
            elseif j == N(i)+3-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i)+a_1)-y(j-N(i)-2+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-N(i)-1+a_1))));
            elseif j == 2*(N(i)+1)-1-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i)+a_1)-y(j-(N(i)+2)+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+1)+a_1))));
            elseif j == 2*(N(i)+1)-a_1+c_5
              dydt(j,1) = 0;
            else
              dydt(j,1) = (-A*(y(j-N(i)+a_1)-y(j-(N(i)+2)+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+1)+a_1))));
            end
        end
      elseif n(p(i,1),2) == 2 && n(p(i,2),2) == 2
        a_1=2;
        for j = (1+c_5):((N(i)+1)-a_1+c_5)               
            if j == 1+c_5
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+3+j-a_1)-y(N(i)+1+j-a_1));
                for k=1:length_junc
                        if n(p(i,1),1) == juncs_data(k,1)
                            juncs_data(k,2) = juncs_data(k,2)-y(N(i)+1+j-a_1);
                            P_j1 = y(juncs_data(k,3));
                        end
                end
            elseif j == N(i)-1+c_5
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+3+j-a_1)-y(N(i)+1+j-a_1));
                for k = 1:length_junc
                        if n(p(i,2),1) == juncs_data(k,1)
                            juncs_data(k,2) = juncs_data(k,2)+y(N(i)+3+j-a_1);
                            P_j2 = y(juncs_data(k,3));
                        end
                end
            else
              dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+3+j-a_1)-y(N(i)+1+j-a_1));
            end
         end
         for j = (c_5+(N(i)+2)-a_1):(c_5+(2*(N(i)+1))-a_1)
            if j == N(i)+2-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i)-1+a_1)-P_j1)/(dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/((P_j1))));
            elseif j == N(i)+3-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i)-1+a_1)-P_j1)/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1))+a_1)/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-N(i)-2+a_1))));
            elseif j == 2*(N(i)+1)-1-a_1+c_5
              dydt(j,1) = (-A*(P_j2-y(j-(N(i)+3)+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+2)+a_1))));
            elseif j == 2*(N(i)+1)-a_1+c_5
              dydt(j,1) = (-A*(P_j2-y(j-(N(i)+3)+a_1))/(dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(P_j2)));
            else
              dydt(j,1) = (-A*(y(j-N(i)-1+a_1)-y(j-(N(i)+3)+a_1))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+2)+a_1))));
            end
        end
      elseif n(p(i,1),2)==2 && n(p(i,2),2)==3
        a_1 = 1;
        for j = (1+c_5):((N(i)+1)-a_1+c_5)
            if j == 1+c_5 
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+3+j-a_1)-y(N(i)+1+j-a_1));
                for k = 1:length_junc
                        if n(p(i,1),1) == juncs_data(k,1)
                            juncs_data(k,2) = juncs_data(k,2)-y(N(i)+1+j-a_1);
                            P_jj = y(juncs_data(k,3));
                        end
                end
            elseif j == N(i)-a_1+c_5
                y(N(i)+2+j) = n(p(i,2),3);
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j)-y(N(i)+j));
            elseif j == N(i)+1-a_1+c_5
                dydt(j,1) = (-(C2(j-c_6))/(A*dx(i)))*(y(N(i)+1+j)-y(N(i)+j));
            else
                dydt(j,1) = (-(C2(j-c_6))/(2*A*dx(i)))*(y(N(i)+2+j)-y(N(i)+j));
            end
        end
         for j = (c_5+(N(i)+2)-a_1):(c_5+(2*(N(i)+1))-a_1)
            if j == N(i)+2-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i))-P_jj)/(dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/((P_jj))));
            elseif j == N(i)+3-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i))-P_jj)/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-N(i)-1))));
            elseif j == 2*(N(i)+1)-1-a_1+c_5
              dydt(j,1) = (-A*(y(j-N(i))-y(j-(N(i)+2)))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+1)))));
            elseif j == 2*(N(i)+1)-a_1+c_5
              dydt(j,1) = 0;
            else
              dydt(j,1) = (-A*(y(j-N(i))-y(j-(N(i)+2)))/(2*dx(i)))-((F*(C2(j-c_6-(N(i)+1)+a_1))/(2*(A)*d))*(((y(j))*abs((y(j))))/(y(j-(N(i)+1)))));
            end
        end
      end

    c_5 = j;
    c_6 = c_6+N(i)+1;
end

%% Define dydt for junction pressure states
if ~isempty(junc_indexes)

for i = 1:length_junc
     dydt(juncs_data(i,3)) = (C2(c_4+i)/v)*(juncs_data(i,2));
end

end
end