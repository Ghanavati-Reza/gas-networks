function  [t_total,y_total,junc_data] = ode_solver(tstart,tend,t_e,N,p,n,dx,sg,d,F,T)

[row_n,col_n] = size(n);
[row_p,col_p] = size(p);

%% finding number of junctions and state variables

c1 = 1;
junc_data = [];
for i = 1:row_n
    if n(i,2) == 2
        junc_data(c1) = i;
        c1 = c1+1;
    end
end
if ~isempty(junc_data)
jun=max(size(junc_data));
end
c2 = 0;
for i = 1:row_p
    for j = 1:col_p
        for k = 1:jun
            if p(i,j) == junc_data(k)
                c2 = c2+1;
            end
        end
    end
end

state_number = 2*(sum(N)+row_p)-c2+jun;

%% Initializing states
x0 = zeros(state_number,1);
c3 = 0;
a1 = 2*(sum(N)+row_p)-c2;
c4 = 0;
c5 = 1;
node_state = [];
for i=1:row_p
    if n(p(i,1),2) == 1
        if n(p(i,2),2) == 3
            x0(1+c3:(N(i)+1+c3),1) = n(p(i,1),3);
            node_state(c5,1) = n(p(i,1),1);
            node_state(c5,2) = 1+c3;
            c5 = c5 + 1;
            c4 = N(i)+1;
        elseif n(p(i,2),2) == 2
            x0(1+c3:(N(i)+c3),1) = n(p(i,1),3);
            node_state(c5,1) = n(p(i,1),1);
            node_state(c5,2) = 1+c3;
            c5 = c5 + 1;
            c4 = N(i);
        end
    elseif n(p(i,1),2) == 2
        if n(p(i,2),2) == 3
            x0(1+c3:(N(i)+c3),1) = 0.8*max(n(:,3));
            c4 = N(i);
        elseif n(p(i,2),2) == 2
            x0(1+c3:(N(i)-1+c3),1) = 0.8*max(n(:,3));
            c4 = N(i)-1;
        end
    end
    if n(p(i,2),2) == 3
        x0(c4+1+c3:c4+(N(i)+1)+c3-1,1) = 0;
        x0(c4+(N(i)+1)+c3,1) = n(p(i,2),3);
        node_state(c5,1) = n(p(i,2),1);
        node_state(c5,2) = c4+(N(i)+1)+c3;
        c5 = c5 + 1;
    else
        x0(c4+1+c3:c4+(N(i)+1)+c3,1) = 0;
    end
  c3 = c3+c4+(N(i)+1);
end
x0(c3+1:state_number,1) = 0.8*max(n(:,3));
%% solving the ODEs
c6=1;
[row_event,col_event] = size(t_e);
tic
for i = 1:row_event
t_event_current = t_e(i,1);
tspan = [tstart tend];
options=odeset('Events',@(t,y) event(t,y,t_event_current));
[t,Y,te,ye,ie] = ode23s(@(t,y) ode_fun(t,y,N,p,n,dx,sg,d,F,T), tspan, x0, options);
% Defining the event

n(t_e(i,2),3) = n(t_e(i,2),3)+t_e(i,3)*n(t_e(i,2),3);
x0 = ye(end,:);
x0(node_state(node_state(:,1)==t_e(i,2),2)) = n(t_e(i,2),3);
tstart = te(end);
if c6==1
y_total=Y;
t_total=t;
c6=2;
else
y_total=[y_total;Y(2:end,:)];
t_total=[t_total;t(2:end,:)];
end
end
tspan = [tstart tend];
[t,Y] = ode23s(@(t,y) ode_fun(t,y,N,p,n,dx,sg,d,F,T), tspan, x0);
y_total=[y_total;Y(2:end,:)];
t_total=[t_total;t(2:end,:)];
toc

function [value,isterminal,direction]=event(t,y,t_e)
%%the value that we want to be zero%%
if t>t_e(1)
    value=0;
else
    value=1;
end
isterminal=1;         %%halt integration when event occurs
direction=0;          %%the zero can be approached from either direction
end

end
