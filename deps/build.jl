#!/usr/bin/env julia
using BinDeps
import Compat
import Compat.Sys

@BinDeps.setup

imagemagick = library_dependency("imagemagick", aliases = ["ImageMagick"])

provides(AptGet, "imagemagick", imagemagick, os = :Linux)
provides(Yum, "ImageMagick", imagemagick, os = :Linux)
provides(Pacman, "imagemagick", imagemagick, os = :Linux)

if Compat.Sys.isapple()
    if Pkg.installed("Homebrew") === nothing
        error("Homebrew julia package not installed, " *
              "please run Pkg.add(\"Homebrew\")")
    end
    using Homebrew
    provides(Homebrew.HB, "imagemagick", imagemagick, os = :Darwin)

elseif Compat.Sys.iswindows()
    using WinRPM
    provides(WinRPM.RPM, "imagemagick", imagemagick, os = :Windows)
end

@BinDeps.install Dict(:imagemagick => :imagemagick)
