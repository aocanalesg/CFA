function [loc] = locate(xi,nx,x_grid)

i_low = 0;
i_up = nx+1;

while (i_up-i_low > 1)
    i_mid = (i_up + i_low)/2;
        if (xi > x_grid(round(i_mid)))
          i_low = i_mid;
		else
		  i_up = i_mid;
		end
end

if (xi <= x_grid(1))
    loc = 1;
elseif (xi >= x_grid(nx)) 
    loc = nx-1; 	 
else
    loc = round(i_low);
end