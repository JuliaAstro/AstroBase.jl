using Unitful

import Unitful: km, s, kg, 𝐋, 𝐓, Length, Time

export km, s, kps

const kps = km/s

@derived_dimension Velocity 𝐋/𝐓

VectorKM = typeof([1.0]*km)
VectorKPS = typeof([1.0]*kps)
