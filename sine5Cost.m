function cost = sine5Cost(sineParams,data)
    % Inputs:
        % sineParams = 3 sinewave parameters [1x3]
        % data = matrix of y values and corresponding circle positions [2xN]
            % Where first row is y values, and second row is trigger positions (1-5)

    % Sine function to fit (5 discrete points from 0pi to 2pi - lines up with data)
    sine5Out = @(sineParams) sineParams(1) + (sineParams(2)*sin(sineParams(3) + (linspace(0,2*pi,5)) ));
 
    % Get 5 sine points from 3 parameters
    y = sine5Out(sineParams);
    % Remove any NaNs from data
    data(:,isnan(data(2,:))) = [];
    % Sum of squared deviations between model and data
    cost = nansum((data(1,:) - y(data(2,:))).^2);
end