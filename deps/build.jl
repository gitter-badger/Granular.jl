#!/usr/bin/env julia
using BinDeps

@BinDeps.setup

imagemagick = library_dependency("imagemagick")

if is_apple()
    using Homebrew
    provides(Homebrew.HB, "imagemagick", imagemagick, os = :Darwin)
end

@BinDeps.install Dict(:imagemagick => :imagemagick)
