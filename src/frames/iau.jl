using ..Bodies

using ReferenceFrameRotations: angle_to_dcm, ddcm

for name in ALL_NAMES
    endswith(name, "Barycenter") && continue

    frame = Symbol("IAU", name)
    cname = Symbol("iau_", lowercase(name))
    @eval begin
        struct $frame <: RotatingFrame end
        const $cname = $frame()
        from_sym(::Val{$frame}) = $cname
        export $frame, $cname

        function Rotation(::ICRF, ::$frame, ep::Epoch)
            body = $(Symbol(name))()
            euler = euler_angles(body, ep)
            euler′ = collect(euler_rates(body, ep))

            m = angle_to_dcm(euler..., :ZXZ)
            m′ = ddcm(m, euler′)

            Rotation{icrf, $cname}(m, m′)
        end

        Rotation(::$frame, ::ICRF, ep::Epoch) = inv(Rotation(icrf, $cname, ep))
    end
end

