using AstroTime: Epoch, Period, unit, value
# using DataInterpolations: CubicSpline
using DataInterpolations

export TimeSeries

struct TimeSeries{S,TT,U,PT,T} <: AbstractArray{T,1}
    epoch::Epoch{S,TT}
    time::Vector{Period{U,PT}}
    data::Vector{T}
    # interp::BSplineInterpolation{Array{DType,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},true,DType}
    interp
    function TimeSeries(epoch::Epoch{S,TT}, time, data) where {S,TT}
        @show time_array = collect(time)
        @show U = unit(time_array[1])
        @show PT = eltype(time_array[1])
        @show T = eltype(data)
        # interp = BSplineInterpolation(data, float(value.(time)), 3, :ArcLen, :Average)
        interp = CubicSpline(data, float(value.(time)))
        @show typeof(epoch)
        @show typeof(time_array)
        @show typeof(data)
        new{S,TT,U,PT,T}(epoch, time_array, data, interp)
    end
end

function (ts::TimeSeries)(p::Period)
    p < ts.time[1] && throw(ArgumentError("`p` is too small"))
    p > ts.time[end] && throw(ArgumentError("`p` is too large"))
    return ts.interp(value(p))
end

function (ts::TimeSeries)(ep::Epoch)
    ep < ts.epoch && throw(ArgumentError("`ep` is too small"))
    ep > (ts.epoch + ts.time[end]) && throw(ArgumentError("`ep` is too large"))
    return ts.interp(value(unit(ts.time[1])(ep - ts.epoch)))
end

Base.getindex(ts::TimeSeries, idx) = ts.data[idx]
Base.size(ts::TimeSeries) = size(ts.data)
Base.eltype(ts::TimeSeries) = eltype(ts.data)
