
function kernel_b(arr_jacobian, arr_Ni1_0, arr_Ni1_s, arr_wvol, arr_x, n1, n_rows)
setmetatable(_ENV, { __index=math })  
function_f = require 'function_f' 
setmetatable(_ENV, { __index=function_f }) 


local n1 = math.floor(n1)
local n_rows = math.floor(n_rows)
local contribution = 0.0
for g1 = 1, n1 do
local x = arr_x[g1]
local wvol = arr_wvol[g1]
local Ni_0 = arr_Ni1_0[g1]
local Ni_u = arr_Ni1_s[g1]
local ggg = g1
local j11 = arr_jacobian[ggg]
local Ni_x = Ni_u*j11
contribution = Ni_0*wvol*f(x) + contribution
end
output = contribution
return output

end

function kernel_a(arr_jacobian, arr_Ni1_0, arr_Ni1_s, arr_Nj1_0, arr_Nj1_s, arr_wvol, arr_x, n1, n_cols, n_rows)
setmetatable(_ENV, { __index=math })  


local n1 = math.floor(n1)
local n_cols = math.floor(n_cols)
local n_rows = math.floor(n_rows)
local contribution = 0.0
for g1 = 1, n1 do
local x = arr_x[g1]
local wvol = arr_wvol[g1]
local Ni_0 = arr_Ni1_0[g1]
local Ni_u = arr_Ni1_s[g1]
local Nj_0 = arr_Nj1_0[g1]
local Nj_u = arr_Nj1_s[g1]
local ggg = g1
local j11 = arr_jacobian[ggg]
local Ni_x = Ni_u*j11
local Nj_x = Nj_u*j11
contribution = Ni_x*Nj_x*wvol + contribution
end
output = contribution
return output

end
