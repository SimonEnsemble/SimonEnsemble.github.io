### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 4193a928-a8f0-11ed-2167-fd4d9b5550bc
begin
	import Pkg; Pkg.activate()
	using OrdinaryDiffEq, CairoMakie
end

# ╔═╡ add74719-83e4-47f3-a397-d8911fe6f05b
update_theme!(fontsize=24)

# ╔═╡ 6182330e-3047-4608-95aa-6b0fc3fec484
begin
	# specify tank geometry
	H = 4.0 # tank height, m
	Lb = 5.0 # bottom base length, m
	Lt = 2.0 # top base length, m
	
	# specify resistance to flow
	c = 1.0 # awkward units
	
	# inlet flow rate (could be funciton of time)
	qᵢ = 1.5 # m³/s
	
	# initial liquid level
	h₀ = 0.0 
end

# ╔═╡ bf1da0a2-f66b-46fc-a982-ade30856c040
begin
	tspan = (0.0, 155.0) # solve for 0 s to 150 s
	
	# area from a helicopter view, m²
	A(h) = (h/H * Lt + (1 - h/H) * Lb) ^ 2
	
	# right-hand-side of ODE
	rhs(h, _, t) = (qᵢ - c * sqrt(h)) / A(h)
	
	# DifferentialEquations.jl syntax
	prob = ODEProblem(rhs, h₀, tspan)
	h = solve(prob, Tsit5())
end

# ╔═╡ c7e52d59-0084-4376-a8a3-8b29623aaa73
begin
	ts = range(0.0, stop=152, length=300)
	# easy as this to compute solution at an array of times!
	hs = h.(ts)

	fig = Figure()
	ax  = Axis(fig[1, 1], xlabel="time, t [s]", ylabel="liquid level, h(t) [m]")
	lines!(ts, hs, linewidth=4)
	hlines!(H, linestyle=:dash)
	xlims!(0, 151)
	save("../../static/blog/tank_prob/numerical_soln.png", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═4193a928-a8f0-11ed-2167-fd4d9b5550bc
# ╠═add74719-83e4-47f3-a397-d8911fe6f05b
# ╠═6182330e-3047-4608-95aa-6b0fc3fec484
# ╠═bf1da0a2-f66b-46fc-a982-ade30856c040
# ╠═c7e52d59-0084-4376-a8a3-8b29623aaa73
