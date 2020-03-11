function CalcBraggsLaw(distance::Float64, wave_length::Float64 , theta::Float64)
(2 * distance * sin(deg2rad(theta))) / wave_length
end

function  FilterBraggLaw(n::Float64, toler::Float64 = .01 )
    if abs(round(Int64, n) - n ) < toler
        return 1
    end 
    return 0
end


