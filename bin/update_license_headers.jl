using Dates

package = "AstroBase"
start_year = 2018
current_year = year(now())

header =
"""
#
# Copyright (c) $start_year-$current_year Helge Eichhorn and the $package.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
"""

for (root, dirs, files) in walkdir("src")
    for file in files
		splitext(file)[end] != ".jl" && continue
        fp = joinpath(root, file)
		str = open(f->read(f, String), fp, "r")
		if !startswith(str, "#\n# Copyright")
			println("Added header for: ", fp)
			open(fp, "w") do f
				write(f, header * str)
			end
		else
			n = length(header)
			str[1:n] == header && continue
			println("Updated header for: ", fp)
			open(fp, "w") do f
				write(f, header * str[n+1:end])
			end
		end
    end
end
