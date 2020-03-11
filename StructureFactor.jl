include("AtomicScattering.jl")

function CalcStrutureFactor(plane::Tuple{Int, Int, Int}, unit_cell::Array{Tuple{String, Tuple{Float64, Float64, Float64}}}, ScatteringVector::Float64)
    factor = 0 
    x,y,z = plane 
    sin_portion = 0 
    cos_portion = 0 
    for element_pos_pair in unit_cell
        h , k , l = element_pos_pair[2]

        scattering_factor = ApporoximateScatteringFactor(ScatteringVector, element_pos_pair[1])
    
        sine = scattering_factor * sin(2 * pi * (h*x + k*y+ l*z ))
        cose = scattering_factor * cos(2 * pi * (h*x + k*y+ l*z ))
       
        sin_portion += sine
        cos_portion += cose
       
    end
    
    return sin_portion^2 + cos_portion^2
end 

function CalcScatteringVector(ratio::Float64)
    return (4 * pi *  ratio)
end 