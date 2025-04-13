# MySplines.jl

This repository contains a spline implementation aiming to be fast to construct and call.

It mainly exports `Spline1D`, with a fair amount of options, but the simplest
use is as follows:
```
using MySplines
x = 1.0:1:10
y = @. cos(x^2 / 9)
spl = Spline1D(x, y)
spl(4.56)
```

Note: If I were to redo this, I would probably start with
[`BasicInterpolators.jl`](https://github.com/markmbaum/BasicInterpolators.jl).
