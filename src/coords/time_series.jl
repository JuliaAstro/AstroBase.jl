using AstroTime: Epoch, Period, unit, value
using DataInterpolations: BSplineInterpolation

export TimeSeries

struct TimeSeries{TType, PType, DType, Scale, Unit, IntType} <: AbstractArray{DType, 1}
    epoch::Epoch{Scale, TType}
    time::Vector{Period{Unit, PType}}
    data::Vector{DType}
    interp::IntType
end

function TimeSeries(epoch, time, data)
    time_array = collect(time)
    interp = BSplineInterpolation(data, float(value.(time)), 3, :ArcLen, :Average)
    TimeSeries(epoch, time_array, data, interp)
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
