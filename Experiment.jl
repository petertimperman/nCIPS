include("StructureFactor.jl")
using Pkg
Pkg.add("Plots")
Pkg.add("PyPlot")
using Plots

NaCl = Array{Tuple{String, Tuple{Float64, Float64, Float64}}}(undef, 8)

NaCl[1]  = ("Na1+",(0.000000 , 0.000000  ,0.000000 ))
NaCl[ 2]  = ("Na1+",(0.000000 , 0.500000  ,0.500000 ))
NaCl[ 3]  = ("Na1+",(0.500000 , 0.000000  ,0.500000 ))
NaCl[ 4]  = ("Na1+",(0.500000 , 0.500000  ,0.000000 ))
NaCl[ 5]  = ("Cl1-",(0.000000 , 0.000000  ,0.000000 ))
NaCl[ 6]  = ("Cl1-",(0.000000 , 0.500000  ,0.500000 ))
NaCl[ 7]  = ("Cl1-",(0.500000 , 0.000000  ,0.500000 ))
NaCl[8]  = ("Cl1-",(.500000 , 0.500000  ,0.500000 ))

d = 5.638 / sqrt(1^2 + 1 ^ 2 + 1 ^ 2) 
dd = 1/(2d)
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
    
    f = CalcStrutureFactor(miller_plane, NaCl, CalcScatteringVector(dd))
    diffracto_gram[i] = f 
    print(theta)
    print(" ")
    println(f)
    i = i + 1 
end
Plots.pyplot()
plot(diffracto_gram)