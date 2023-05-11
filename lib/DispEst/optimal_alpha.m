function alpha = optimal_alpha(img_param, met_param)
%OPTIMAL_ALPHA Look-up table for optimal alpha

switch(met_param.alg)
    case 'ack'
    switch(met_param.p)
        case 1
            switch(img_param.snr)
                case 5
                    alpha = 1225;
                case 20
                    alpha = 2050;
                case 60
                    alpha = 3000;
            end
        case 2
            switch(img_param.snr)
                case 5
                    alpha = 2;
                case 20
                    alpha = 6;
                case 60
                    alpha = 7.5;
            end
    end
   
    case 'ncc'
        switch(img_param.snr)
            case 5
                alpha = 0.9;
            case 20
                alpha = 0.9;
            case 60
                alpha = 1.5;
        end
end
