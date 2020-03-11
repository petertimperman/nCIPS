include("UnitCell.jl")
include("StructureFactor.jl")
include("BraggPeaks.jl")

using Pkg
Pkg.add("Plots")
Pkg.add("PyPlot")
using Plots
using Random 
NaCl_atoms = Array{Tuple{String, Tuple{Float64, Float64, Float64}}}(undef, 8)

NaCl_atoms[1]  = ("Na1+",(0.000000 , 0.000000  ,0.000000 ))
NaCl_atoms[ 2]  = ("Na1+",(0.000000 , 0.500000  ,0.500000 ))
NaCl_atoms[ 3]  = ("Na1+",(0.500000 , 0.000000  ,0.500000 ))
NaCl_atoms[ 4]  = ("Na1+",(0.500000 , 0.500000  ,0.000000 ))
NaCl_atoms[ 5]  = ("Cl1-",(0.000000 , 0.000000  ,0.500000 ))
NaCl_atoms[ 6]  = ("Cl1-",(0.000000 , 0.500000  ,0.000000 ))
NaCl_atoms[ 7]  = ("Cl1-",(0.500000 , 0.000000  ,0.000000 ))
NaCl_atoms[8]  = ("Cl1-",(.500000 , 0.500000  ,0.500000 ))

NaCl = UnitCell(5.69169356, 5.69169356,5.69169356, NaCl_atoms)
n_parts =10000
parts = Array{UnitCell}(undef,n_parts)

parts[1] = NaCl

for i in 2:n_parts
    scale = rand(-1000:1000)/10000
   
    parts[i] =  UnitCell(5.69169356+scale, 5.69169356,5.69169356, NaCl_atoms)
end


wave_length  =1.0 

step_size = .02 
start = 1 
steps = convert( Int, (89-1) / .02 ) 
diffracto_gram = Array{Float64}(undef , steps)
miller_plane = (1,1,1)
tetha = start
i = 1 
while i <= steps
    global i, step_size , diffracto_gram
    theta = tetha + ( i *step_size)
    println(i)
    diffracto_gram[i] = 0 
    for part in parts
        n = CalcBraggsLaw(part.d_a,wave_length,theta)
        
        if FilterBraggLaw(n) == 1 
            d = part.d_a/ sqrt(1^2 + 1 ^ 2 + 1 ^ 2) 
            dd = 1/(2d)
            
            f = CalcStrutureFactor(miller_plane, part.atoms, CalcScatteringVector(dd))
           
            diffracto_gram[i] = f + diffracto_gram[i]
        end
  
    end 
   
    # print(theta)
    # print(" ")
    # println(f)
    i = i + 1 
end

Plots.pyplot()
plot(diffracto_gram)