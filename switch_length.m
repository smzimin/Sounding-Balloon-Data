function [ left, right, fgr, num ] = switch_length( shirota )

switch shirota
     case 0
        left = -18;
        right = 18;
        fgr = 1;
        num = 1;
    case 1
        left = -18;
        right = -2;
        fgr = 1;
        num = 1;
    case 2
        left = -2;
        right = 0;
        fgr = 1;
        num = 2;
    case 3
        left = 0;
        right = 2;
        fgr = 2;
        num = 1;
    case 4
        left = 2;
        right = 4;
        fgr = 2;
        num = 2;
    case 5
        left = 4;
        right = 6;
        fgr = 3;
        num = 1;
    case 6
        left = 6;
        right = 9;
        fgr = 3;
        num = 2;
    otherwise
        
end

end

