using DelimitedFiles 

table = readdlm("data/constants/atomic_scattering_constants.csv", ','  , header= true)
vals = table[1]

elemental_constants = Dict{String, Array{Float64}}()
for row in eachrow(vals)
    elemental_constants[row[1]] = row[2:10]
end

function ApporoximateScatteringFactor(q::Float64, element::String)
   
    constants = elemental_constants[element]
    
    factor = constants[9] 
    i = 1 

    while i <= 4

        a = constants[2i-1] 
        b = constants[2i]    
      
        e_term = (-b * (q/(4 * pi))^2)
        factor = factor + (a*exp(e_term))
        
        i = i + 1
    end 

    return factor 
end