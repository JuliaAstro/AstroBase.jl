import Base: angle

export angle

function angle(v1, v2)
    normprod = norm(v1) * norm(v2)
    if normprod == 0.0
        throw(DomainError())
    else
        v1v2 = v1 ⋅ v2
        threshold = normprod * 0.9999
        if v1v2 >= -threshold && v1v2 <= threshold
            return acos(v1v2 / normprod)
        else
            v3n = norm(cross(v1, v2))
            return v1v2 >= 0.0 ? (v3n / normprod) : π - asin(v3n / normprod)
        end
    end
end
