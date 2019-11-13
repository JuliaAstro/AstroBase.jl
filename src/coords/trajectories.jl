using AstroTime: Period, value, unit
using DataInterpolations: AbstractInterpolation, BSplineInterpolation
using SmoothingSplines
using LinearAlgebra

import Base: getindex, lastindex, show
import AstroBase: state

export Trajectory, initial, final, state, events, times,
    LogEntry, count_id, id, epoch, detector,
    TimeSeries

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

struct LogEntry
    id::Int
    detector::Symbol
    time::Float64
    ep::Epoch
end

id(l::LogEntry) = l.id
detector(l::LogEntry) = l.detector
epoch(l::LogEntry) = l.epoch
count_id(idx, log) = count(x->id(x) == idx, log)

struct Trajectory
    initial
    final
    times::Vector
    splines::Vector{SmoothingSpline}
    # components::Vector{TrajectoryComponents}
    array::Matrix
    events::Vector{LogEntry}
end

initial(tra::Trajectory) = tra.initial
final(tra::Trajectory) = tra.final
events(tra::Trajectory) = tra.events
times(tra::Trajectory) = tra.times

function Trajectory(initial, final, events=LogEntry[])
    Trajectory(initial, final, Vector(undef, 0), SmoothingSpline[], Matrix(undef, 0, 0), events)
end

function Trajectory(initial, final, t, vectors, events=LogEntry[])
    n = length(vectors[1])
    splines = Array{SmoothingSpline}(undef, n)
    arr = permutedims(hcat(vectors...))
    m = size(arr)[2]
    for i = 1:m
        splines[i] = fit(SmoothingSpline, t, arr[:,i], 0.0)
    end
    Trajectory(initial, final, t, splines, arr, events)
end

function show(io::IO, tra::Trajectory)
    println(io, "Trajectory")
    println(io, " Start date: $(epoch(initial(tra)))")
    println(io, " End date:   $(epoch(final(tra)))")
end

function interpolate(tra::Trajectory, time)
    SmoothingSplines.predict.(tra.splines, float(time))
end

function state(tra::Trajectory, time)
    if isempty(tra.times)
        error("Trajectory does not contain data points.")
    end
    rv = interpolate(tra, time)
    r = rv[1:3]
    v = rv[4:6]
    ep1 = epoch(initial(tra)) + time * seconds
    f = frame(initial(tra))
    b = body(initial(tra))
    return State(ep1, r, v, f, b)
end

function state(tra::Trajectory, ep::Epoch)
    ep0 = epoch(initial(tra))
    time = typeof(ep0)(ep) - ep0
    tra(seconds(time))
end

state(tra::Trajectory, time::Period) = state(tra, get(seconds(time)))

(tra::Trajectory)(time) = state(tra, time)

getindex(tra::Trajectory, idx) = state(tra, tra.times[idx])
lastindex(tra::Trajectory) = length(tra.times)
