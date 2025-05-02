function [position,isterminal,direction] = Event_CompleteMelting(t,T,input)

position = T(end)-input.tol;  % stop when ice volume is met
isterminal = 1; % terminate the solver 
direction = 0; % both directions

end