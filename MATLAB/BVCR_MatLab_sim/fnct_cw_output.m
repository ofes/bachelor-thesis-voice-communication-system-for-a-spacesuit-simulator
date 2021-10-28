%% FUNCTION: Display optimal tilting angles and alignment
% Displays the optimal tilting angle and alignment for a PV generator in at a given mission loaction on Earth.

% Organization:     OeWF (Austrian Space Forum)
% Author:           Omar Filip El Sendiouny
% Project:          Serenity BU-COMMs
% Date:             15.12.2020
% Version:          1 

function [] = fnct_cw_output(show_command_window_output, ...
    N_md, gamma_S, alpha_S, date_array, beta, t_S, W_G_el, SOC_curr)
    count = zeros(3, 1);
    align_str = {'FLAT', 'SOUTH', 'NORTH'};
    if (show_command_window_output)
        fprintf('*******************************************\n');
        fprintf('[OUTPUT] PV GENERATOR INSTALLATION:\n');
        fprintf('\n');
        fprintf('Mission day | Opt. tilt angle | Orientation\n');
        fprintf('-------------------------------------------\n');
        for n = 1:N_md
            
            fprintf(datestr(date_array(n)));
            if(gamma_S(n, t_S == 12) == 90) % Sun is at its zenith for t_S = 12h
                count(1) = count(1) + 1;
                fprintf(' |');
                fprintf('     %.2f deg    |     ', round(beta(n),2));
                fprintf(align_str{1});
                fprintf('\n');
                
            elseif(alpha_S(n, t_S == 12) == 0) % Sun is visible in the south for t_S = 12h
                count(2) = count(2) + 1;
                fprintf(' |');
                
                if(beta(n) < 10)
                    fprintf('     %.2f deg    |    ', round(beta(n),2));
                else
                    fprintf('    %.2f deg    |    ', round(beta(n),2));
                end
                
                fprintf(align_str{2});
                fprintf('\n');
                
            elseif(alpha_S(n, t_S == 12) == 180) % Sun is visible in the north for t_S = 12h
                count(3) = count(3) + 1;
                fprintf(' |');
                
                if(beta(n) < 10)
                    fprintf('     %.2f deg    |    ', round(beta(n),2));
                else
                    fprintf('    %.2f deg    |    ', round(beta(n),2));
                end
                
                fprintf(align_str{3});
                fprintf('\n');
            end
        end
        [~, count_idx] = max(count);
        fprintf('-------------------------------------------\n');
        
        fprintf('    Mean    |');
        
        if(mean(beta(n)) < 10)
            fprintf('     %.2f deg    |    ', round(mean(beta),2));
        else
            fprintf('    %.2f deg    |    ', round(mean(beta),2));
        end
        
        if(strcmp(align_str{count_idx},'FLAT'))
            fprintf(' ');
            fprintf(align_str{count_idx});
        else
            fprintf(align_str{count_idx});
        end
        fprintf('\n\n');
        
        fprintf('***************************************************\n');
        fprintf('[OUTPUT] PV GENERATOR ENERGY YIELD:\n');
        fprintf('\n');
        fprintf('Applies for mean tilting angle and orientation.\n');
        fprintf('The energy yield is the electrical energy yield.\n');
        fprintf('SOC full shows if the battery could be fully\n');
        fprintf('charged during the day.\n')
        fprintf('\n');
        fprintf('Mission day | Energy yield | SOC full\n');
        fprintf('-------------------------------------\n');
        for n = 1:N_md
            fprintf(datestr(date_array(n)));
            fprintf(' |');
            
            if(W_G_el(n) < 1000) % display energy yield
                fprintf('   %0.2f Wh  |',round(W_G_el(n), 2));
            elseif(W_G_el(n) < 100)
                fprintf('    %0.2f Wh  |',round(W_G_el(n), 2));
            elseif(W_G_el(n) < 10)
                fprintf('     %0.2f Wh  |',round(W_G_el(n), 2));
            else % > 10000
                fprintf('  %0.2f Wh  |',round(W_G_el(n), 2));
            end
            
            if(max(SOC_curr(n,:)) == 1)
                fprintf('   YES   ')
            else
                fprintf('    NO   ')
            end
            
          
           
            fprintf('\n');
        end
        fprintf('\n');
        fprintf('***************************************************\n');
    end      
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





