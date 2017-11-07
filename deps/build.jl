#!/usr/bin/env julia
using BinDeps
using Compat

@BinDeps.setup

imagemagick = library_dependency("imagemagick", aliases = ["ImageMagick"])

if is_linux()
    provides(AptGet, "imagemagick", imagemagick, os = :Linux)
    provides(Yum, "ImageMagick", imagemagick, os = :Linux)
    provides(Pacman, "imagemagick", imagemagick, os = :Linux)

elseif is_apple()
    if Pkg.installed("Homebrew") === nothing
        error("Homebrew julia package not installed, " *
              "please run Pkg.add(\"Homebrew\")")
    end
    using Homebrew
    provides(Homebrew.HB, "imagemagick", imagemagick, os = :Darwin)

elseif is_windows()
    using WinRPM
    provides(WinRPM.RPM, "ImageMagick", imagemagick, os = :Windows)
end

@BinDeps.install Dict(:imagemagick => :imagemagick)
