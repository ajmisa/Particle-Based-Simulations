% Lees het CSV-bestand in
filename = 'energies.csv';  % Vervang door de naam van je .csv-bestand
data = readtable(filename);

% Haal de 4e kolom op
fC = data{:, 4};  % De haakjes {} zorgen ervoor dat je de data als array krijgt

min_y = min(fC);
max_y = max(fC);

% Set y-axis limits 10% below the minimum and 10% above the maximum
ylim_lower = min_y - 0.0001 * abs(min_y);
ylim_upper = max_y + 0.0001 * abs(max_y);

for i = 2:length(fC)
    percent_diff = 100 * (fC(i) - fC(i-1)) / abs(fC(i-1));  % Percent difference
    
    % Annotate the plot with the percentage difference
    text(i-0.5, fC(i), sprintf('%.2f%%', percent_diff), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Plot de 4e kolom
figure;
plot(fC);  % '-o' voegt markers toe aan de punten op de plot
title('Plot van de 4e kolom');
ylim([ylim_lower ylim_upper]);
xlabel('Rijnummer');
ylabel('Waarde');
grid on;

