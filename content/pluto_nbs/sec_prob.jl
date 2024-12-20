### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ dbbcbf84-be7d-11ef-3384-2b2a69580c7a
begin
	import Pkg; Pkg.activate()
	using CairoMakie, Graphs, GraphMakie, ColorSchemes, Printf, NetworkLayout, Test
end

# ╔═╡ d5c07c6f-0bbc-428d-a1b1-dde60d406cfe
n = 5

# ╔═╡ a0d12780-5ab2-48d0-9854-4fa73a268ae4
begin
	g = DiGraph(2 * n)
	for k = 1:n
		add_edge!(g, k, n+k)
		add_edge!(g, n+k, n+k)
		for ℓ = 1:n
			if ℓ > k
				add_edge!(g, k, ℓ)
			end
		end
	end
end

# ╔═╡ d4580bec-fe82-4283-981b-86df877ef5a1
function viz_g(
	g::DiGraph, 
	P::Union{Nothing, Matrix{Float64}}=nothing
)
	colormap = reverse(ColorSchemes.acton)
	fig = Figure()#size=(450, 450))
	ax = Axis(fig[1, 1])
	
	hidedecorations!(ax)
	hidespines!(ax)

	p_to_label(p::Float64) = p < 1e-8 ? "0" : @sprintf("%.2f", p)
	function node_to_label(i::Int)
		if i < nv(g) / 2 + 1
			return "$i"
		else
			return "$(i-n)*"
		end
	end

	if ! isnothing(P)
		elabels = [p_to_label(P[ed.src, ed.dst]) for ed in edges(g)]
		edge_color = [
			get(colormap, P[ed.src, ed.dst], (0.0, 0.5)) 
				for ed in edges(g)
		]
	else
		elabels = nothing
		edge_color = nothing
	end

	gp = graphplot!(
		g, 
		node_strokewidth=3,
		edge_color=edge_color,
		arrow_size=25,
		node_color="white",
		node_strokecolor="black",
		elabels=elabels,
		node_size=35,
		ilabels=[node_to_label(v) for v in vertices(g)],
		ilabels_fontsize=25,
		arrow_shift=:end,
		elabels_fontsize=12,
		curve_distance=1.0,
		edge_width=4,
		# layout=Stress()
		layout=Spring(C=10.0)
	)
	if ! isnothing(P)
		Colorbar(
			fig[1, 2], colormap=colormap, 
			label="probability", limits=(0, 0.5)
		)
	end
	# if ! isnothing(π)
	# 	Colorbar(
	# 		fig[1, 3], colormap=ColorSchemes.Blues, label="probability", limits=(0.0, maximum(π))
	# 	)
	# end
	fig
end


# ╔═╡ 478ae2a4-c55f-46ad-b597-5045e898fd37
begin
	P = zeros(2 * n, 2 * n)
	for k = 1:n # time at which k-th candidate observed
		# enter absorbing state, i.e. this the best candidate u ever see
		P[k, n + k] = k / n
		# absorbing state'
		P[n + k, n + k] = 1
		
		for ℓ = 1:n # time at which ℓ-th candidate observed
			if ℓ > k
				P[k, ℓ] = k / (ℓ * (ℓ - 1))
			end
		end
	end
	P
end

# ╔═╡ c1da3a78-fcfb-4d68-8c14-461cced050e7
@test all(P * ones(2*n) .≈ ones(2*n))

# ╔═╡ 2aa85c66-e0bb-4a11-b1e3-b4438bdb53de
@test all(P .≤ 1.0)

# ╔═╡ debeac52-6c15-430a-99a0-5ab7ce018292
viz_g(g, P)

# ╔═╡ 9e255ead-d995-4da5-b64d-615003cc55c8
begin
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		titlesize=25,
		xlabel="arrival time of candidate i, k",
		ylabel="arrival time of candidate i+1, ℓ",
		xticks=1:2*n, 
		yticks=1:2*n,
		aspect=DataAspect()
	)

	# hack to make p = 0.0 map to white
	colormap = vcat(ColorSchemes.grays[1.0], ColorSchemes.viridis[0.0:0.01:1.0])
	hm = heatmap!(P, colorrange=(0, 1), colormap=colormap)
	Colorbar(fig[1, 2], hm, label="transition probability", labelsize=20)
	fig
end

# ╔═╡ 2b685bbc-a6f7-4693-98c0-7481e2cfe571
π₀ = [i == 1 ? 1.0 : 0.0 for i = 1:2*n]

# ╔═╡ 9e3251d3-d44d-4df8-ae14-2550167fc3cb
π₀' * P^2

# ╔═╡ a4d206a5-46aa-4527-b713-f939f3f02688
barplot((π₀' * P * P * P)[:])

# ╔═╡ c99e89a9-1480-4e35-8621-a5b07cc336fe


# ╔═╡ Cell order:
# ╠═dbbcbf84-be7d-11ef-3384-2b2a69580c7a
# ╠═d5c07c6f-0bbc-428d-a1b1-dde60d406cfe
# ╠═a0d12780-5ab2-48d0-9854-4fa73a268ae4
# ╠═d4580bec-fe82-4283-981b-86df877ef5a1
# ╠═478ae2a4-c55f-46ad-b597-5045e898fd37
# ╠═c1da3a78-fcfb-4d68-8c14-461cced050e7
# ╠═2aa85c66-e0bb-4a11-b1e3-b4438bdb53de
# ╠═debeac52-6c15-430a-99a0-5ab7ce018292
# ╠═9e255ead-d995-4da5-b64d-615003cc55c8
# ╠═2b685bbc-a6f7-4693-98c0-7481e2cfe571
# ╠═9e3251d3-d44d-4df8-ae14-2550167fc3cb
# ╠═a4d206a5-46aa-4527-b713-f939f3f02688
# ╠═c99e89a9-1480-4e35-8621-a5b07cc336fe
