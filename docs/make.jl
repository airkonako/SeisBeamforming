using Documenter, SeisBeamforming

makedocs(;
    modules=[SeisBeamforming],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/kura-okubo/SeisBeamforming.jl/blob/{commit}{path}#L{line}",
    sitename="SeisBeamforming.jl",
    authors="kurama",
    assets=String[],
)

deploydocs(;
    repo="github.com/kura-okubo/SeisBeamforming.jl",
)
