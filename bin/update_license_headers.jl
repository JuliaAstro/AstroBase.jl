using Dates

package = "AstroBase"
maintainer = "Helge Eichhorn"
start_year = 2018
current_year = year(now())

header =
"""
#
# Copyright (c) $start_year-$current_year $maintainer and the $package.jl contributors
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
            str1 = header * str[n+1:end]
            str1 = replace(str, r"#\n([a-z])"=>s"#\n\n\1"; count=1)
            str1 == str && continue
            println("Updated header for: ", fp)
            open(fp, "w") do f
                write(f, str1)
            end
        end
    end
end
