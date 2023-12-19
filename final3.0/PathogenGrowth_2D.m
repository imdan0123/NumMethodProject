% function to compute the increase in icidence of a plant pathogen (or any
% pathogen with stationary hosts).  This code is 2-D because it predicts 
% the growth of the pathogen for a 2-D array of hosts each with a 
% population that grows in size with time.  Additionally, it calculates 
% airborne transmission of disease between hosts in the 2D array.  It is 
% a modified version of the PathogenGrowth_0D function used in lab 09.
%
% The primary equations governing this are compound interest equations for
% growth:
%
% dS/dt = -beta*S*I+dP/dt
%
% dL/dt = beta*S*I-mu_L*L+dE/dt
%
% dI/dt = mu_L*L*I-mu_I*I
%
% dP/dt = dB/dt + d(leaf)/dt 
%
% dR/dt = mu_I*I
%
% dE/dt = e
%
% dF/dt = Gamma*exp(alpha*I) - F*R_frac
%
% with:
% S    = susceptible fraction of population (susceptible tissue)
% beta = rate of infection growth for healthy population (fraction per day)
% L    = fraction of tissue infected and latent (e.g., dormant before infection)
% I    = fraction of tissue infected and producing inoculum
% R    = fraction of tissue recovered (or removed) from population
% P    = size of the total population (plant surface area)
% B    = surface area of berries
% E    = amount (fraction of population) of introduced disease from external sources
% F    = size of the spreading population, e.g. sporulating for a fungus 
% t    = time
% mu_L = rate of decrease in latent population (fraction per day)
% mu_I = rate of decrease in infectious population (fraction per day)
% e    = rate of import from external sources
%
% inputs: vine (structure containing the initial size of susceptible population); beta; mu;
% tspan (days to simulate array);
% output: S,L,I,R,P,E,time (vector of simulation times), and B

function [vine] = PathogenGrowth_2D(vine,beta_max,mu_L_target,mu_I,A,...
    eta,kappa,xi,Gamma,alpha,T,U,V,tspan)

%declare global variables
global NpX NpY Nsteps

% set parameters in a cell array
p{1} = beta_max; %(max rate of new infections)
p{2} = 1/mu_I;   %(inverse length of the infectious period in days)
p{3} = T;        %(temperature in C)
p{4} = tspan;    %(time in days array)
p{5} = A;        %(total plant surface area at reference time)
p{6} = sqrt(U.^2+V.^2);  %windspeed
p{7} = atand(V./U);      %wind direction
p{8} = eta;       %release fraction scale factor
p{9} = kappa;     %release fraction scale factor   
p{10}= xi;       %release fraction offset
p{11}= Gamma;    %spore production multiple
p{12}= alpha;   %spore production 2nd factor
% iterations = 0; 
% flag = 0; 
% declare function handles
odefun = @(t,y,e,g) SLIRPE_model(t,y,e,g,p);

% loop over timesteps (starting at 2)
for t=2:Nsteps

    disp(['day=',num2str(tspan(t),'%.2f'),' infected plants=',int2str(sum([vine.IsInfect]))])
    dt=tspan(t)-tspan(t-1); %timestep

    %update list of infected vines
    ActiveVines = find([vine.IsInfect]);
    %initialize deposition flux for timestep
    DepFlux_sum=zeros(1,NpX*NpY);

    for idx = ActiveVines
        %calculate positions of other vines relative to each infected plant
        Xplume = [vine.X]-vine(idx).X;
        Yplume = [vine.Y]-vine(idx).Y;
        
        DepFlux = GaussianPlumeDep(Xplume,Yplume,p{6}(t),p{7}(t),...
            vine(idx).S(t-1),vine(idx).F(t-1));
        NNaNInd = ~isnan(DepFlux);
        DepFlux_sum(NNaNInd) = DepFlux_sum(NNaNInd) + DepFlux(NNaNInd);
    end
    
    %now for all vines
    for i=1:NpX
        for j=1:NpY
            cnt=i+(j-1)*NpX; %index counter for vectorized vine structure
            
            %check if vines have just become latent, if so calculate mu_L 
            %and flip their LatentSwitch so that the calc only happens 1 time.
            if((vine(cnt).L(t-1) > 1e-8) && (vine(cnt).LatentSwitch == false))
                vine(cnt).mu_L = latentperiod(t,dt,Nsteps,mu_L_target,...
                    zeros(size(T)),T);
                vine(cnt).LatentSwitch = true;
            end
% Initialize plant area with variability
% for i = 1:NpX
%     for j = 1:NpY
%         cnt = i + (j-1) * NpX;
%         
%         % Calculate variability based on a normal distribution
%         variability = 0.2 * P_i_ave * randn;
%         
%         % Initialize plant area with variability
%         vine(cnt).P(t-1) = P_i_ave + variability;
% 
%         % Other initializations...
%     

            % set initial conditions for time integration (could use deal here)
            y0(1) = vine(cnt).B(t-1); %(amount of population, surface area, that is berries)
            y0(2) = vine(cnt).P(t-1); %(total population, surface area, including berries and leaves)
            y0(3) = vine(cnt).S(t-1); %(initial susceptible population fraction)
            y0(4) = vine(cnt).L(t-1); %(initial latent population fraction)
            y0(5) = vine(cnt).I(t-1); %(initial infectious population fraction)
            y0(6) = vine(cnt).R(t-1); %(initial recovered population fraction)
            y0(7) = vine(cnt).E(t-1); %(initial external population fraction)
            y0(8) = vine(cnt).F(t-1); %(size of the spreading population)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % integrate using a time integration method from class
            %%%% INSERT YOUR CODE HERE

%             % MATLAB Built in ODE 45 built in funtion method
%             opt = odeset('RelTol',1e-6,'AbsTol',1e-8);  % Adjust tolerances
%             [t, y] = ode45(@(t,y) odefun(t, y, DepFlux_sum(cnt), vine(cnt).mu_L), ...
%                [tspan(t-1), tspan(t)], y0, opt);
% 
%             % Extract the solution at the final time point
%             y = y(end, :);

             % Use of the Runge-Kutta 4 method for time intigration from
             % class 

%              function [y] = TimeInt(odefun, t, dt, y0, DepFlux, mu_L)
%              % Runge-Kutta 4th order time integration
%              k1 = dt * odefun(t, y0, DepFlux, mu_L);
%              k3 = dt * odefun(t + 0.5 * dt, y0 + 0.5 * k2, DepFlux, mu_L);
%              k4 = dt * odefun(t + dt, y0 + k3, DepFlux, mu_L);
% 
%              y = y0 + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
%              end




            % time int call for a single step advancement 
            %[y] = TimeInt(odefun, t, dt, y0, DepFlux_sum(cnt), vine(cnt).mu_L);

  [y] = TimeInt(odefun, t, dt, y0, DepFlux_sum(cnt), vine(cnt).mu_L);


            % 
            % for time integration, note for
            %%%% odefun to work as given above your call needs to look like:
            %
            %[y] = TimeInt(odefun,t,dt,y0,DepFlux_sum(cnt),vine(cnt).mu_L);
            %
            %NOTE: recognize that you are only integrating 1 time step!
            %your routine can be more general than that but recognize that
            %this point is in the middle of a time loop!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % set outputs
            vine(cnt).B(t) = y(1);
            vine(cnt).P(t) = y(2);
            vine(cnt).S(t) = y(3);
            vine(cnt).L(t) = y(4);
            vine(cnt).I(t) = y(5);
            vine(cnt).R(t) = y(6);
            vine(cnt).E(t) = y(7);
            vine(cnt).F(t) = y(8);

            %define a threshold for dispersal to start (equiv to 2.5mm diameter
            %sporulating colony (Calonnec et al)
            if(vine(cnt).I(t)>=0.25^2/4*pi/A)
                vine(cnt).IsInfect=true;
            end
            %turn off dispersal if we fall below the above size
            if(vine(cnt).I(t)<0.25^2/4*pi/A && vine(cnt).IsInfect==true)
                vine(cnt).IsInfect=false;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%        RECOMMENDED LOCATION FOR YOUR SCOUTING ROUTINE   
    
%        if 
%         scout_cnt = scout_cnt+1; %update our scouting storage counter
%         for sc = 1:Scout.Nsens %for each scout
%             switch lower(Scout.Type)
% 
%                 case 'random'
%                     PosToScout = randperm(Sout.Numpts); %pick random plants to scout
%                     %initialize time counter and position
%                     SX = 0; SY = 0; ScTime = 0;
% 
%                         for idx = 1:Scout.Numpts
%                             ScDistance = sqrt((Sx-vine(PosToScout(idx)).X)^2 + ...
%                                 (SY-vine(PosToScout(idx)).Y)^2);
%                             %update time and position
%                             ScTime = ScTime + ScDistance/Scout.Speed;
%                             SX = vine(PosToScout(idx)).X;
%                             SY = vine(posToScout(Idx)).Y;
% 
%                             if(ScTime <= 3600)%as  long as we still have time
%                                 InfectSize = vine(PosToScout(idx)).L(t)*A;
%                                 ScoutSize = max([(40*Scout.Speed)^2/4*pi pi]);
%                                 if(InfectSize >= ScoutSize)
%                                     Scout.Success(sc, scout_cnt) = true;
%                                     return %stop if success
%                                 end
%                             else
%                                 break
%                             end
%                         end
%             end
%         end
%     
           if tspan(t) <= (tspan(t-1) +1)
        scout_cnt = scout_cnt+1; %update our scouting storage counter
        for sc = 1:Scout.NSens %for each scout
            switch lower(Scout.Type)

                case 'random'
                    PosToScout = randperm(Scout.Numpts); %pick random plants to scout
                    %initialize time counter and position
                    SX = 0; SY = 0; ScTime = 0;

                        for idx = 1:Scout.Numpts
                            ScDistance = sqrt((SX-vine(PosToScout(idx)).X)^2 + ...
                                (SY-vine(PosToScout(idx)).Y)^2);
                            %update time and position
                            ScTime = ScTime + ScDistance/Scout.Speed;
                            SX = vine(PosToScout(idx)).X;
                            SY = vine(PosToScout(idx)).Y;

                            if(ScTime <= 3600)%as  long as we still have time
                                InfectSize = vine(PosToScout(idx)).L(t)*A;
                                ScoutSize = max([(40*Scout.Speed)^2/4*pi pi]);
                                if(InfectSize >= ScoutSize)
                                    Scout.Success(sc, scout_cnt) = true;
                                    [cost] = CostCalc(Scout);
                                    return %stop if success
                                end
                            else
                                break
                            end
                        end

            end
        end
       end
   if tspan(t) >= length(tspan) - 1
   % If the pathogen is not detected by the end of the simulation, calculate the cost for the last day
   [cost] = CostCalc(Scout);
end
end
    
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

end
% Runge kutta four method used for time intigration
 function [y] = TimeInt(odefun, t, dt, y0, DepFlux, mu_L)
             % Runge-Kutta 4th order time integration
             k1 = dt * odefun(t, y0, DepFlux, mu_L);
             k2 = dt * odefun(t + 0.5 * dt, y0 + 0.5 * k1, DepFlux, mu_L);
             k3 = dt * odefun(t + 0.5 * dt, y0 + 0.5 * k2, DepFlux, mu_L);
             k4 = dt * odefun(t + dt, y0 + k3, DepFlux, mu_L);

             y = y0 + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
           


 end

  function [cost] = CostCalc(Scout)
    penaltyday = 10; % 1st day to start counting

    %find first successful day of scouting
    day = find(sum(Scout.Success,1)>0,1);

    % if not found, cost is max
    if(isempty(day))
        day = 61;
    end
    if(day<penaltyday)
        cost = Scout.NSens * Scout.HrCost * day;
    else
        cost = Scout.NSens * Scout.HrCost * day + (day-penaltyday) * Scout.DayCost;
    end
    fprintf("Total Cost: $%0.4f\n",cost)
end

 %[y] = TimeInt(odefun, t, dt, y0, DepFlux_sum(cnt), vine(cnt).mu_L);
 %[y] = TimeInt(odefun, t, dt, y0, DepFlux_sum(cnt), vine(cnt).mu_L);