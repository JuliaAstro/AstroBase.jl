using AstroTime: Period, Epoch, TimeScale
using SmoothingSplines
using LinearAlgebra
using Roots: find_zero

import Base: getindex, lastindex, show
import AstroBase: state
import AstroTime: unit

export TypedTrajectory, Trajectory, initial, final, state, events, times,
    LogEntry, count_id, id, epoch, detector

abstract type AbstractTrajectory{Scale, Frame, Body} end

timescale(s::AbstractTrajectory{Scale}) where {Scale} = Scale
frame(s::AbstractTrajectory{_S, Frame}) where {_S, Frame} = Frame
body(s::AbstractTrajectory{_S, _F, Body}) where {_S, _F, Body} = Body

function transpose_aoa(input::Vector{Vector{T}}) where T
    n = length(input[1])
    output = Vector{Vector{T}}(undef, n)
    matrix = hcat(input...)
    for i = 1:n
        output[i] = matrix[i, :]
    end
    return output
end

struct TypedTrajectory{Scale, Frame, Body, T, TType, Unit, PType, NT} <: AbstractTrajectory{Scale, Frame, Body}
    series::NT
    events::Vector{Event{Scale, TType}}
    data::Vector{Vector{T}}
    time::Vector{Period{Unit, PType}}
    function TypedTrajectory(epoch::Epoch{Scale, TType},
                             time::Vector{Period{Unit, PType}},
                             data::Vector{Vector{T}};
                             frame::Frame=icrf,
                             body::Body=earth,
                             names::Vector{Symbol}=Symbol[]) where {Scale, Frame, Body, TType, Unit, PType, T}
        columns = transpose_aoa(data)
        if isempty(names)
            n = length(columns)
            if n == 3
                names = [:x, :y, :z]
            elseif n == 6
                names = [:x, :y, :z, :vx, :vy, :vz]
            else
                throw(ArgumentError("Explicit `names` need to provided if `data` does not have either 3 or 6 columns."))
            end
        end
        series = (; zip(names, [TimeSeries(epoch, time, c) for c in columns])...)
        events = Vector{Event{Scale, TType}}[]
        new{Scale::TimeScale, frame::AbstractFrame, body::CelestialBody,
            T, TType, Unit, PType, typeof(series)}(series, events, data, time)
    end
end

unit(tra::TypedTrajectory) = unit(tra.time[1])
Base.float(p::Period) = float(value(p))

function find_events!(tra::TypedTrajectory, detectors::Vector{Detector})
    for detector in detectors
        for (i, t) in enumerate(tra.time[1:end-1])
            t1 = tra.time[i+1]
            s0 = sign(detect(detector, t, tra))
            s1 = sign(detect(detector, t1, tra))
            detected = s1 > s0
            detected || continue

            t_event = find_zero(t -> detect(detector, t * unit(tra), tra), (value(t), value(t1)))
            evt = Event(epoch(tra) + t_event * unit(tra), nameof(typeof(detector)))
            evt in tra.events || push!(tra.events, evt)
        end
    end
    sort!(tra.events, lt=(a, b)->isless(a.epoch, b.epoch))
end

function (tra::TypedTrajectory)(t)
    typ = promote_type(eltype.(values(tra.series))...)
    out = Vector{typ}(undef, length(tra.series))
    for (i, ts) in enumerate(tra.series)
        out[i] = ts(t)
    end
    return out
end

epoch(tra::TypedTrajectory) = tra.series[1].epoch

function keplerian(tra::TypedTrajectory, t)
    rv = tra(t)[1:6]
    return keplerian(rv[1:3], rv[4:6], grav_param(body(tra)))
end

function State(tra::TypedTrajectory, t::Period)
    rv = tra(t)[1:6]
    return State(epoch(tra) + t, rv[1:3], rv[4:6],
                 scale=timescale(tra), frame=frame(tra), body=body(tra))
end

function State(tra::TypedTrajectory, ep::Epoch)
    rv = tra(t)[1:6]
    return State(ep, rv[1:3], rv[4:6],
                 scale=timescale(tra), frame=frame(tra), body=body(tra))
end

function KeplerianState(tra::TypedTrajectory, t::Period)
    ele = keplerian(tra, t)
    return KeplerianState(epoch(tra) + t, ele...,
                          scale=timescale(tra), frame=frame(tra), body=body(tra))
end

function KeplerianState(tra::TypedTrajectory, ep::Epoch)
    ele = keplerian(tra, ep)
    return KeplerianState(ep, ele...,
                          scale=timescale(tra), frame=frame(tra), body=body(tra))
end

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
