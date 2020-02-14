function playprint

fprintf('\nCounter: ')
for i=1:10
      if i>1
          for j=0:log10(i-1)
              fprintf('\b'); % delete previous counter display
          end
      end
      fprintf('%d', i);
      pause(.05); % allows time for display to update
end
fprintf('\n')
