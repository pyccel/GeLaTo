
function kernel_b(arr_jacobian, arr_Ni1_0, arr_Ni1_s, arr_Ni2_0, arr_Ni2_s, arr_wvol, arr_x, arr_y, n1, n2, n_rows)
setmetatable(_ENV, { __index=math })  
function_f = require 'function_f' 
setmetatable(_ENV, { __index=function_f }) 


local n1 = math.floor(n1)
local n2 = math.floor(n2)
local n_rows = math.floor(n_rows)
local contribution = 0.0
for g1 = 1, n1 do
for g2 = 1, n2 do
local gg = g2 + n2*(g1 - 1)
local x = arr_x[gg]
local y = arr_y[gg]
local wvol = arr_wvol[gg]
local Ni_0 = arr_Ni1_0[g1]*arr_Ni2_0[g2]
local Ni_u = arr_Ni1_s[g1]*arr_Ni2_0[g2]
local Ni_v = arr_Ni1_0[g1]*arr_Ni2_s[g2]
local ggg = 4*gg - 3
local j11 = arr_jacobian[ggg]
local ggg = 4*gg - 2
local j12 = arr_jacobian[ggg]
local ggg = 4*gg - 1
local j21 = arr_jacobian[ggg]
local ggg = 4*gg
local j22 = arr_jacobian[ggg]
local Ni_x = Ni_u*j11 + Ni_v*j12
local Ni_y = Ni_u*j21 + Ni_v*j22
contribution = Ni_0*wvol*f(x, y) + contribution
end
end
output = contribution
return output

end

function kernel_a(arr_jacobian, arr_Ni1_0, arr_Ni1_s, arr_Ni2_0, arr_Ni2_s, arr_Nj1_0, arr_Nj1_s, arr_Nj2_0, arr_Nj2_s, arr_wvol, arr_x, arr_y, n1, n2, n_cols, n_rows)
setmetatable(_ENV, { __index=math })  


local n1 = math.floor(n1)
local n2 = math.floor(n2)
local n_cols = math.floor(n_cols)
local n_rows = math.floor(n_rows)
local contribution = 0.0
for g1 = 1, n1 do
for g2 = 1, n2 do
local gg = g2 + n2*(g1 - 1)
local x = arr_x[gg]
local y = arr_y[gg]
local wvol = arr_wvol[gg]
local Ni_0 = arr_Ni1_0[g1]*arr_Ni2_0[g2]
local Ni_u = arr_Ni1_s[g1]*arr_Ni2_0[g2]
local Ni_v = arr_Ni1_0[g1]*arr_Ni2_s[g2]
local Nj_0 = arr_Nj1_0[g1]*arr_Nj2_0[g2]
local Nj_u = arr_Nj1_s[g1]*arr_Nj2_0[g2]
local Nj_v = arr_Nj1_0[g1]*arr_Nj2_s[g2]
local ggg = 4*gg - 3
local j11 = arr_jacobian[ggg]
local ggg = 4*gg - 2
local j12 = arr_jacobian[ggg]
local ggg = 4*gg - 1
local j21 = arr_jacobian[ggg]
local ggg = 4*gg
local j22 = arr_jacobian[ggg]
local Ni_x = Ni_u*j11 + Ni_v*j12
local Ni_y = Ni_u*j21 + Ni_v*j22
local Nj_x = Nj_u*j11 + Nj_v*j12
local Nj_y = Nj_u*j21 + Nj_v*j22
contribution = contribution + wvol*(Ni_x*Nj_x + Ni_y*Nj_y)
end
end
output = contribution
return output

end
