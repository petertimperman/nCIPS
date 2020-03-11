struct UnitCell
    d_a::Float64
    d_b::Float64
    d_c::Float64
    atoms::Array{Tuple{String, Tuple{Float64, Float64, Float64}}}
end



# function UnitCellFromCif(cif_file_path::String)
#     if !occursin(r"\.cif$", cif_file_path)
#         d_a = 0.0
#         d_a = 0.0
#         d_a = 0.0 

#         for line in eachline(cif_file_path)
        
#         end 
#     end
# end