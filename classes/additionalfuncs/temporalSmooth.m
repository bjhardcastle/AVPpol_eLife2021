%temporal smoothing

        function stack_in_raw_out = temporalSmooth(stack_in_raw, window_size,window_beginning)
        stack_in_raw = stack_in_raw(:,:,window_beginning:end);
            
            [y_dim, x_dim, slices] = size(stack_in_raw);
            
            %preallocate
            stack_in_raw_out = zeros(y_dim, x_dim, slices);
            
            for i = 1:size(stack_in_raw,3)-1
                
                %make sure you are not running out of bounds at end
                if i+window_size > size(stack_in_raw, 3)
                    
                    window_size = window_size-1;
                    
                end
                
                stack_in_raw_out(:,:,i) =  mean(stack_in_raw(:,:,i:(i+window_size)),3);
                
            end
            
        end