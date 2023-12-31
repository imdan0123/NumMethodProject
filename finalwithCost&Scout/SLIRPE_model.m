%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SLIRPE function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIRPE_model (pronounced like the 7-11 slushy drink) combines all the 
% SLIR equations and the plant growth equations with an external source 
% equation in one system. The dependent variables are stored in y as:
% y(1) = B (amount of population, surface area, that is berries)
% y(2) = P (total population, surface area, including berries and leaves: P = S+L+I+R)
% y(3) = S (susceptible population)
% y(4) = L (latent population)
% y(5) = I (infectious population)
% y(6) = R (recovered/removed population)
% y(7) = E (amount of new infections from external sources)
% y(8) = F (size of the spreading population, e.g. sporulating for a fungus) 
%
% and the parameters in a cell array:
% p{1} = beta_max (max rate of colony growth/new infections)
% p{2} = mu_i (inverse length of the infectious period in days)
% p{3} = T (array of temperature in C)
% p{4} = day (array of times in units of days)
% p{5} = A (total plant surface area at reference time)
% p{6} = Windspd (windspeed)
% p{7} = Winddir (wind direction, currently not used here)
% p{8} = eta     (release fraction scale factor)
% p{9} = kappa   (release fraction scale factor)  
% p{10}= xi      (release fraction offset)
% p{11}= Gamma   (spore production multiple)
% p{12}= alpha   (spore production 2nd factor)
%
% and the time input (idx) should be an integer for the iteration number
% Note that E is not calculated by the function (only integrated in time)
function [dydt] = SLIRPE_model(idx,y,e,mu_L,p)
    %assign parameters
    beta_max = p{1};
    mu_I     = p{2}; 
    T        = p{3};
    day      = p{4};
    A        = p{5};
    Windspd  = p{6};
    Winddir  = p{7};
    eta      = p{8};
    kappa    = p{9};       
    xi       = p{10};
    Gamma    = p{11};
    alpha    = p{12};

    %assign variables
    B = y(1);
    P = y(2);
    S = y(3);
    L = y(4);
    I = y(5);
    R = y(6);
    E = y(7);
    F = y(8);

    %calculated parameters
    if(ceil(idx)==floor(idx)) %when we are at an interger step
        T_used = T(idx);
        day_used = day(idx);
        mu_L_used = mu_L(idx);
        m_used = Windspd(idx);
    else %when we are at a half step (for rk4)
        idx = floor(idx);
        T_used = 0.5*(T(idx)+T(idx+1));
        day_used = 0.5*(day(idx)+day(idx+1));
        mu_L_used = 0.5*(mu_L(idx)+mu_L(idx+1));
        m_used = 0.5*(Windspd(idx)+Windspd(idx+1));
    end
    beta = beta_max*Sall_temp_effect(T_used); %pathogen growth rate
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dydt(1) = Pb_model(vine(idx).P(t-1), T_used, F, A);
%     dydt(2) = dydt(1) + Pl_model(vine(idx).P(t-1), T_used, E, mu_L_used, dt, mu_I, I);
    T_e = -0.35968 + (T_used * 0.10789) - ( (T_used^2) * 0.00214); 
    Dpl = (T_e * day_used *1.33) ;
    dydt(1) = T_e * ((B * 0.1724) - (B^2 * 0.0000212));
    dydt(2) = (Dpl + dydt(1))/A;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %disease variables
    dydt(3) = -beta*S*I + dydt(2)/A; % change in S
    dydt(4) = beta*S*I - mu_L_used*L + e; % change in L
    dydt(5) = mu_L_used*L - mu_I*I; % change in I
    dydt(6) = mu_I*I; % change in R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % External source (E) and Spreading population (F)
dydt(7) = e; % External source
if I > 0
    dydt(8) = Gamma * exp(alpha * I * A) - F * R; % Implement F function
else
    dydt(8) = 0;
end
    if(I==0)%spore production shouldn't start before infection (quirk of exponential curve fit)
        dydt(8) = 0;
    else
        %YOUR CODE GOES HERE for our F function
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%
