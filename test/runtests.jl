using Splines
using Test

using Statistics
using Random
using Splines
using PyPlot
using Profile

using Dierckx
#using DierckxSpline
#Dierckx = DierckxSpline


function iterate_mine(x, y)
	for i=1:1000
		s = Splines.Spline1D(x,y)
	end
end

function iterate_mine(x, y, spline)
	for i=1:1000
		s = Splines.Spline1D(x, y, spline)
	end
end


function iterate_Dierckx(x, y)
	for i=1:1000
		s = Dierckx.Spline1D(x,y,bc="zero")
	end
end


function iterate_mine_integrate(spline, xlo, xhi; N=1000)
	for i=1:N
		I = Splines.integrate(spline, xlo, xhi)
	end
end


function iterate_Dierckx_integrate(spline, xlo, xhi; N=1000)
	for i=1:N
		I = Dierckx.integrate(spline, xlo, xhi)
	end
end



function speedtester(x,y)
	N = length(x)
	spline = Splines.Spline1D(x, y)
	iterate_mine(x,y,spline)
	@time iterate_mine(x,y,spline)
	iterate_mine(x,y)
	@time iterate_mine(x,y)
	iterate_Dierckx(x,y)
	@time iterate_Dierckx(x,y)
end

function speedtester_integrate(x, y, xlo, xhi; N=1000)
	s = Splines.Spline1D(x, y)
	spl = Dierckx.Spline1D(x, y, bc="zero")
	iterate_mine_integrate(s, xlo, xhi; N=N)
	@time iterate_mine_integrate(s, xlo, xhi; N=N)
	#Profile.clear()
	#@profile iterate_mine_integrate(s, xlo, xhi; N=N)
	#Profile.print()
	iterate_Dierckx_integrate(spl, xlo, xhi; N=N)
	@time iterate_Dierckx_integrate(spl, xlo, xhi; N=N)
end


function testminimal()
	println("testminimal:")
	N = 4
	x = Array(range(0, stop=1, length=N))
	y = rand(N)
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	speedtester(x, y)

	xx = range(-0.1, stop=1.1, length=1000)
	figure()
	scatter(x,y)
	plot(xx, s.(xx))
	plot(xx, spl.(xx), color="0.85")
	ylim(-0.1,1.1)
end


function testminimal_with_derivative()
	println("testminimal with derivative:")
	N = 10
	x = Array(range(0, stop=1, length=N))
	y = rand(N)
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	xx = range(-0.1, stop=1.1, length=1000)
	figure()
	scatter(x,y)
	plot(xx, s.(xx))
	plot(xx, spl.(xx), color="0.85")
	for i=1:4
		xderiv = rand()
		yderiv = Splines.derivative(s, xderiv)
		xxd = range(xderiv - 0.05, stop=xderiv + 0.05, length=50)
		yyd = @. yderiv * (xxd - xderiv) + s(xderiv)
		plot(xxd, yyd, color="b", ls="--")
	end
	ylim(-0.1,1.1)
end

function testminimal_loglog()
	println("testminimal log-log with derivative:")
	N = 10
	x = Array(range(0.1, stop=1, length=N))
	y = 0.1 .+ rand(N)
	s = Splines.Spline1Dloglog(x,y)
	s1 = Splines.Spline1Dtrans(x,y,log,identity)
	s2 = Splines.Spline1Dtrans(x,y,identity,log)
	spl = Splines.Spline1D(x,y)

	xx = range(0.05, stop=1.1, length=1000)
	figure()
	scatter(x,y)
	plot(xx, s.(xx), label="log-log")
	plot(xx, s1.(xx), label="log-lin")
	plot(xx, s2.(xx), label="lin-log")
	plot(xx, spl.(xx), color="0.85", label="linear")
	legend()
end


function testminimal_k0()
	println("testminimal k=0:")
	N = 10
	x = Array(range(0, stop=1, length=N+1))
	y = rand(N)
	s0 = Splines.Spline1D(x,y; k=0)
	#s3 = Splines.Spline1D(x,y)

	xx = range(-0.1, stop=1.1, length=1000)
	figure()
	scatter(middle.(x[1:end-1], x[2:end]), y)
	plot(xx, s0.(xx), label="0th-order")
	#plot(xx, s3.(xx), label="3rd-order")
	ylim(-0.1,1.1)
end


function testminimal_k0_loglog()
	println("testminimal k=0 with log-log:")
	N = 10
	#x = Array(range(0.1, stop=1, length=N+1))
	x = Array(10 .^ range(log10(0.1), stop=log10(1), length=N+1))
	y = 0.1 .+ rand(N)
	s0 = Splines.Spline1Dloglog(x,y,k=0)
	s1 = Splines.Spline1Dtrans(x,y,log,identity,k=0)
	s2 = Splines.Spline1Dtrans(x,y,identity,log,k=0)
	s3 = Splines.Spline1Dtrans(x,y,identity,identity,k=0)

	xx = range(0.05, stop=1.1, length=1000)
	figure()
	scatter(middle.(x[1:end-1], x[2:end]), y)
	plot(xx, s0.(xx), label="log-log")
	plot(xx, s1.(xx), label="log-lin")
	plot(xx, s2.(xx), label="lin-log")
	plot(xx, s3.(xx), label="lin-lin")
	legend()
end


function testminimal_k1()
	println("testminimal k=1:")
	N = 10
	x = Array(range(0, stop=1, length=N))
	y = rand(N)
	s1 = Splines.Spline1D(x,y; k=1)
	#s3 = Splines.Spline1D(x,y)

	xx = range(-0.1, stop=1.1, length=1000)
	Splines.derivative.(s1, xx)
	figure()
	scatter(x, y)
	plot(xx, s1.(xx), label="1st-order")
	plot(xx, Splines.derivative.(s1, xx), label="1st-order derivative")
	#plot(xx, s3.(xx), label="3rd-order")
	#ylim(-0.1,1.1)
	legend()
end


function testrand()
	println("testrand:")
	N = 20
	Random.seed!(4)
	x = Array(range(0, stop=1, length=N))
	y = rand(N)
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	speedtester(x, y)

	println(s(0.45))
	println(spl(0.45))

	xx = range(-0.1, stop=1.1, length=10000)
	yy = s.(xx)
	figure()
	scatter(x,y)
	plot(xx, yy)
	plot(xx, spl.(xx), color="0.85")
	ylim(-0.1,1.1)
end


function testoox()
	println("testoox:")
	N = 10
	x = Array{Float64}(undef, N)
	x[1:2] = [0.02, 0.03]
	x[3:N] = range(0.04, stop=1, length=N-2)
	y = 1 ./ x
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	speedtester(x, y)

	println(s(0.45))
	println(spl(0.45))

	xx = range(-0.1, stop=1.1, length=1000)
	yy = s.(xx)
	figure()
	scatter(x,y)
	plot(xx, yy)
	plot(xx, spl.(xx), color="0.85")
	plot(xx, 1 ./ xx, color="0.85",  linestyle="--")
	ylim(-1, 35)
end


function testgauss()
	println("testgauss:")
	N = 10
	x = Array(range(-1, stop=1, length=N))
	y = @. exp(-x^2 / 0.1^2)
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	speedtester(x, y)

	println(s(0.45))
	println(spl(0.45))

	xx = range(-1.1, stop=1.1, length=1000)
	yy = s.(xx)
	figure()
	scatter(x,y)
	plot(xx, yy)
	plot(xx, spl.(xx), color="0.85")
	plot(xx, exp.(-xx.^2 / 0.1^2), color="0.85",  linestyle="--")
end


function testlargearray()
	println("testlargearray:")
	N = 1000
	x = Array(range(-1, stop=1, length=N))
	y = @. exp(-x^2 / 0.1^2)
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	speedtester(x, y)

	println(s(0.45))
	println(spl(0.45))

	xx = range(-1.1, stop=1.1, length=1000)
	yy = s.(xx)
	figure()
	scatter(x,y)
	plot(xx, yy)
	plot(xx, spl.(xx), color="0.85")
	plot(xx, exp.(-xx.^2 / 0.1^2), color="0.85",  linestyle="--")
end


function testspeed()
	N = 500
	x = Array(range(-1, stop=1, length=N))
	y = @. exp(-x^2 / 0.3^2)
	speedtester(x, y)
end


function testintegrate()
	println("testintegrate:")
	N = 500
	x = Array{Float64,1}(range(-1, stop=1, length=N))
	y = @. exp(-x^2 / 0.3^2)
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	println(Splines.integrate(s, -0.9, 0.8))
	println(Dierckx.integrate(spl, -0.9, 0.8))

	speedtester_integrate(x, y, -0.9, 0.8; N=100000)
	speedtester_integrate(x, y, -0.9, -0.88; N=100000)

	#xx = range(-1.1, stop=1.1, length=1000)
	#yy = s(xx)
	#figure()
	#scatter(x,y)
	#plot(xx, yy)
	#plot(xx, spl.(xx), color="0.85")
end


function testintegrate_beyondbndry()
	println("testintegrate:")
	N = 10
	x = Array(range(-1, stop=1, length=N))
	y = rand(N)
	s = Splines.Spline1D(x,y)
	spl = Dierckx.Spline1D(x,y,bc="zero")

	I = Splines.integrate(s, -1.9, 1.4)
	I = Splines.integrate(s, 1.1, 1.4)
	I = Splines.integrate(s, -1.8, -1.4)

	xx = range(-1.1, stop=1.1, length=1000)
	figure()
	scatter(x, y)
	plot(xx, s.(xx))
	plot(xx, spl.(xx), color="0.85")
end


function testallocate()
	println("testallocate:")
        x = 10.0 .^ (-3:0.01:1)
	y = rand(length(x))

	#s = Splines.Spline1D(x, y; k=3)
	#s = Splines.Spline1D(3)
	#spl = Splines.Spline1D()
	#Splines.Spline1D(x, y, s)

	#s = Splines.Spline1Dloglog(x, y)
	s = Splines.Spline1D(x, y)

        xx = rand(10000)

        @time s(xx[1])
        @time s(xx[4])
        @time s(xx[100])

        @time @. s(xx)
        @time @. s(xx)
        @time @. s(xx)
end


function testextrapolations()
	println("testextrapolations:")
	x = 1:0.5:4
	y = randn(length(x))
	s_zero = Splines.Spline1D(x, y, extrapolation=Splines.zero)
	s_boundary = Splines.Spline1D(x, y, extrapolation=Splines.boundary)
	s_linear = Splines.Spline1D(x, y, extrapolation=Splines.linear)
	s_powerlaw = Splines.Spline1D(x, y, extrapolation=Splines.powerlaw)

	xx = 0:0.01:5
	figure()
	title("Extrapolations")
	scatter(x, y)
	plot(xx, s_zero.(xx), label="Zero")
	plot(xx, s_boundary.(xx), label="Boundary")
	plot(xx, s_linear.(xx), label="Linear")
	plot(xx, s_powerlaw.(xx), label="Power-law")
	ylim(minimum(y)-0.5, maximum(y)+0.5)
	legend()
end


function test_inverse()
    z = [0.11109999999999996, 0.36915874194714393, 0.7333961254531127, 1.2287736648072456]
    r = [324.3780626283167, 1008.6025862083494, 1815.1139372153434, 2665.022835526069]

    szr = Splines.Spline1D(z, r, extrapolation=Splines.linear)
    srz = Splines.Spline1D(r, z, extrapolation=Splines.linear)

    @show issorted(z) issorted(r)
    @assert issorted(z)
    @assert issorted(r)


    # First broken test case:
    z1 = 0.8
    r1 = szr(z1)
    r2 = 1944.404
    z2 = srz(r2)
    @show z1 z2 r1 r2

    # Are r1,r2 and z1,z2 sorted the same way?
    @test_broken sign(r1-r2) == sign(z1-z2)


    # bigger differences are OK?
    z3 = 0.801
    r3 = szr(z3)
    z4 = 0.799
    r4 = srz(z4)
    @show z3 r3
    @show z4 r4
    @test r3 > r1
    @test r3 > r2
    @test r4 < r1
    @test r4 < r2


    # Another broken case:
    z5 = 0.895
    r5 = sxy(z5)
    r6 = 2120.08
    z6 = syx(r6)
    @show z5 z6 r5 r6
    @test_broken sign(r6-r5) == sign(z6-z5)
end


@testset "Splines.jl" begin
    testminimal()
    testminimal_loglog()
    testminimal_with_derivative()
    testminimal_k0()
    testminimal_k0_loglog()
    testminimal_k1()
    testrand()
    testoox()
    testgauss()
    testlargearray()
    testspeed()
    testintegrate()
    testintegrate_beyondbndry()
    testallocate()
    testextrapolations()
    test_inverse()
end
