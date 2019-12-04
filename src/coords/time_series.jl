using AstroTime: Epoch, Period, unit, value
using DataInterpolations: BSplineInterpolation

export TimeSeries

struct TimeSeries{Scale, TType, Unit, PType, DType} <: AbstractArray{DType, 1}
    epoch::Epoch{Scale, TType}
    time::Vector{Period{Unit, PType}}
    data::Vector{DType}
    interp::BSplineInterpolation{Array{DType,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},true,DType}
    function TimeSeries(epoch::Epoch{Scale, TType}, time, data) where {Scale, TType}
        time_array = collect(time)
        interp = BSplineInterpolation(data, float(value.(time)), 3, :ArcLen, :Average)
        new{Scale, TType, unit(time[1]), eltype(time[1]), eltype(data)}(epoch, time_array, data, interp)
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
