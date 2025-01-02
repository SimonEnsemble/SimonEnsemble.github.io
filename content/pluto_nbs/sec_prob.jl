### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ dbbcbf84-be7d-11ef-3384-2b2a69580c7a
begin
	import Pkg; Pkg.activate()
	using CairoMakie, Graphs, GraphMakie, ColorSchemes, Printf, NetworkLayout, Test, MakieThemes

	set_theme!(theme_ggthemr(:lilac))
	update_theme!(fontsize=20, linewidth=3, markersize=20)
end

# ╔═╡ be97b1d7-fe22-4112-bb20-19e5bede769b
colors = MakieThemes.GGThemr.ColorTheme[:lilac][:swatch]

# ╔═╡ 4d03b774-6114-4632-88ed-d698e3a54342
bkg_color = MakieThemes.GGThemr.ColorTheme[:lilac][:background]

# ╔═╡ aa08c613-192c-4e89-9bab-6b42fe64918f
probability_cmap = reverse(ColorSchemes.acton)

# ╔═╡ cfffd689-2318-4653-99b8-152c4bc2fdaf
md"## policy independent MC"

# ╔═╡ d5c07c6f-0bbc-428d-a1b1-dde60d406cfe
n = 5

# ╔═╡ a0d12780-5ab2-48d0-9854-4fa73a268ae4
function MC_graph(n::Int)
	g = DiGraph(n + 1)
	# state n: end (absorbing)
	A = n + 1
	add_edge!(g, A, A)
	for t = 1:n
		add_edge!(g, t, A) # could be last
		for t_next = t+1:n
			add_edge!(g, t, t_next)
		end
	end
	return g
end

# ╔═╡ d4580bec-fe82-4283-981b-86df877ef5a1
function viz_g(
	g::DiGraph,
	mode::String="policy",
	P::Union{Nothing, Matrix{Float64}}=nothing;
	t★::Union{Nothing, Int}=nothing
)
	@assert mode in ["policy", "non-policy"]

	# number of candidates
	n = mode == "policy" ? nv(g) - 2 : nv(g) - 1
	
	fig = Figure()#size=(450, 450))
	ax = Axis(fig[1, 1])
	
	hidedecorations!(ax)
	hidespines!(ax)

	p_to_label(p::Float64) = p < 1e-8 ? "0" : @sprintf("%.2f", p)
	function node_to_label(i::Int)
		if i ≤ n
			return "$i"
		end
		if mode == "policy"
			if i == n + 1
				return "L"
			elseif i == n + 2
				return "W"
			end
		end
		if mode == "non-policy" && i == n + 1
			return "E"
		end
	end

	node_to_color(i) = i > n ? colors[6] : colors[5]

	if ! isnothing(P)
		elabels = [p_to_label(P[ed.src, ed.dst]) for ed in edges(g)]
		edge_color = [
			get(probability_cmap, P[ed.src, ed.dst], (0.0, 0.5)) 
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
		node_color=[node_to_color(v) for v in vertices(g)],
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
			fig[1, 2], colormap=probability_cmap, 
			label="transition probability", limits=(0, 1.0)
		)
	end
	if ! isnothing(t★)
		@assert mode == "policy"
		Label(
			fig[1, 1], "t*=$t★", tellwidth=false, tellheight=false,
			halign=0.8, valign=0.8
		)
	end
	save("g_$mode.pdf", fig)
	fig
end

# ╔═╡ 478ae2a4-c55f-46ad-b597-5045e898fd37
function MC_transition_matrix(n::Int)
	P = zeros(n+1, n+1)
	A = n+1 # end, absorbing state
	P[A, A] = 1
	for t = 1:n # time at which i-th candidate observed
		# enter absorbing state, i.e. this the best candidate u ever see
		P[t, A] = t / n
		
		for t_next = t+1:n # time at which ℓ-th candidate observed
			P[t, t_next] = t / (t_next * (t_next - 1))
		end
	end
	P
end

# ╔═╡ 5b2837f4-dfce-48c4-a9cf-36565767581b
function gimme_tick_labels(n::Int, mode::String)
	return vcat(
		["$i" for i = 1:n],
		mode == "policy" ? ["L", "W"] : ["E"]
	)
end

# ╔═╡ 9e255ead-d995-4da5-b64d-615003cc55c8
function viz_P(P::Matrix{Float64}, mode::String)
	@assert mode in ["policy", "non-policy"]
	
	# number of candidates
	n = mode == "policy" ? size(P)[1] - 2 : size(P)[1] - 1

	ticklabels = gimme_tick_labels(n, mode)
	
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		titlesize=25,
		xlabel="current state",
		ylabel="next state",
		xticks=(1:size(P)[1], ticklabels),
		yticks=(1:size(P)[1], ticklabels),
		aspect=DataAspect()
	)

	# hack to make p = 0.0 map to white
	colormap = vcat(bkg_color, probability_cmap[0.0:0.01:1.0])
	#ColorSchemes.viridis[0.0:0.01:1.0])
	hm = heatmap!(P, colorrange=(0, 1), colormap=colormap)
	Colorbar(fig[1, 2], hm, label="transition probability", ticks=[i/4 for i = 0:4])
	save("P_$mode.pdf", fig)
	fig
end

# ╔═╡ a96c1505-5851-4e48-a52b-98a02576a8c3
function viz_π(π::Vector{Float64}, mode::String, t::Int)
	n = mode == "policy" ? length(π) - 2 : length(π) - 1
	
	ticklabels = gimme_tick_labels(n, mode)

	fig = Figure()
	ax = Axis(
		fig[1, 1], xlabel="stage, t", ylabel=rich("ℙ(T", subscript("$t"), "=t)"),
		xticks=(1:length(π), ticklabels)
	)
	ylims!(0, 1)
	barplot!(π, color=[get(probability_cmap, πᵢ) for πᵢ in π])
	fig
end

# ╔═╡ b64569d1-3825-4569-a648-f33b278ff002
function gimmie_π₀(n::Int, mode::String)
	if mode == "non-policy"
		π₀ = zeros(n+1)
	elseif mode == "policy"
		π₀ = zeros(n+2)
	end
	π₀[1] = 1.0
	return π₀
end

# ╔═╡ a45403a5-fc90-44c3-aff8-37d5d3b90421
function viz_πs(P::Matrix{Float64}, mode::String, t_max::Int)
	n = mode == "policy" ? size(P)[1] - 2 : size(P)[1] - 1
	
	π₀ = gimmie_π₀(n, mode)
	
	ticklabels = gimme_tick_labels(n, mode)

	fig = Figure()
	axs = [Axis(
			fig[t+1, 1], 
			ylabel=rich("ℙ(T", subscript("$t"), "=s)")
		) 
		for t = 0:t_max
	]
	linkaxes!(axs...)
	linkyaxes!(axs...)
	axs[end].xlabel = "state, s"
	
	for t = 0:t_max
		ylims!(axs[t+1], 0, 1)
		
		if t == t_max
			axs[t+1].xticks = (1:size(P)[1], ticklabels)
		else
			axs[t+1].xticks = (1:size(P)[1], ["" for _ in ticklabels])
		end
		
		π = π₀' * P ^ t
		
		barplot!(axs[t+1], π[:], color=[get(probability_cmap, πᵢ) for πᵢ in π[:]])
	end
	save("pi_over_time.pdf", fig)
	fig
end

# ╔═╡ 2b685bbc-a6f7-4693-98c0-7481e2cfe571
π₀ = gimmie_π₀(n, "non-policy")

# ╔═╡ e24e6587-f41d-4f05-bb32-260463c9b4dd
t = 2

# ╔═╡ c99e89a9-1480-4e35-8621-a5b07cc336fe
md"## policy-dependent MC
 $t^*$: after that or equal to that stage select the candidate.
"

# ╔═╡ e98898b5-4437-4fe5-92d6-f1fb035331a8
function MC_graph(n::Int, t★::Int)
	L = n + 1
	W = n + 2
	
	g = DiGraph(n+2)
	# self-loops
	add_edge!(g, L, L) 
	add_edge!(g, W, W)
	
	for t = 1:n
		if t ≤ t★
			# reject candidate. can only lose right after.
			add_edge!(g, t, L)
		else
			# accept candidatee could loose or win
			add_edge!(g, t, L)
			add_edge!(g, t, W)
		end
		# ah, can't lose if see candidate at end
		rem_edge!(g, n, L)
		for t_next = t+1:n
			if t ≤ t★
				add_edge!(g, t, t_next)
			end
		end
	end
	return g
end

# ╔═╡ 8c38dcc7-b1f5-4e75-b105-1c61090901b8
g = MC_graph(n)

# ╔═╡ ed545ef2-53d9-4edb-92f2-f672241cf390
function MC_transition_matrix(n::Int, t★::Int)
	L = n + 1
	W = n + 2
	
	P = zeros(n + 2, n + 2)
	# state t: candidate observed
	# state n + 1: lose
	# state n + 2: win
	for t = 1:n # time candidate i observed
		if t ≤ t★
			# possibility we lose, by this candidate, who we reject, being #1
			P[t, L] = t / n
		else
			# possibility we lose, by this candidate, who we accept, not being #1
			P[t, L] = 1 - t / n
			# possibility we win, by this candidate, who we accept, being #1
			P[t, W] = t / n
		end
		
		for t_next = t+1:n # time at which candidate i+1 observed
			if t ≤ t★ # only continue observing candidates under this condition
				P[t, t_next] = t / (t_next * (t_next - 1))
			end
		end
	end
	# absorbing states
	P[L, L] = P[W, W] = 1.0
	return P
end

# ╔═╡ eeae5e01-1b04-4855-a889-3ed05a82f529
P = MC_transition_matrix(n)

# ╔═╡ c1da3a78-fcfb-4d68-8c14-461cced050e7
@test all(P * ones(n+1) .≈ ones(n+1))

# ╔═╡ 2aa85c66-e0bb-4a11-b1e3-b4438bdb53de
@test all(P .≤ 1.0)

# ╔═╡ debeac52-6c15-430a-99a0-5ab7ce018292
viz_g(g, "non-policy", P)

# ╔═╡ c13431c0-1a8a-437d-b947-9b3d8c783fe4
viz_P(P, "non-policy")

# ╔═╡ 9e3251d3-d44d-4df8-ae14-2550167fc3cb
π = π₀' * P ^ t

# ╔═╡ a4d206a5-46aa-4527-b713-f939f3f02688
viz_π(π[:], "non-policy", t)

# ╔═╡ 1e8d615d-25e8-4c2b-be16-4a372670ec0c
t★ = 2 # policy

# ╔═╡ 02f29171-3438-4620-a6f2-8f95624fca7b
g_policy = MC_graph(n, t★)

# ╔═╡ f612ba1a-a3c3-49c1-96cd-41be7b1cfdfe
P_policy = MC_transition_matrix(n, t★)

# ╔═╡ 1eec5172-b76c-4469-b0db-8db3aeac8902
viz_P(P_policy, "policy")

# ╔═╡ 1d909b88-9b4a-45cd-b883-6ec8ab57e09f
@test P_policy * ones(n+2) == ones(n+2)

# ╔═╡ 5f7eb3a5-1263-4054-b618-0d9973b41ada
viz_g(g_policy, "policy", P_policy, t★=t★)

# ╔═╡ cbf80fcc-155b-4cc1-bbd3-8c0024f59fac
π_policy = gimmie_π₀(n, "policy")' * P_policy ^ t

# ╔═╡ 1d7f9db9-98d7-4be0-b90b-959d1ec9be1e
viz_π(π_policy[:], "policy", t)

# ╔═╡ 06484fb8-5ddf-4e9f-9178-1fee8aa33add
viz_πs(P_policy, "policy", t★+1)

# ╔═╡ f5982842-2310-47d4-93dc-3a223f439781
function p_win(n::Int, t★::Int)
	P = MC_transition_matrix(n, t★)
	π_t★ = gimmie_π₀(n, "policy")' * P ^ (t★+1)
	return π_t★[end]
end

# ╔═╡ 37753f25-c165-4e5d-9959-5f5238b94a1a
p_wins = [p_win(n, t★) for t★=0:n-1]

# ╔═╡ 4331d0ea-d709-4798-a4a5-01c8b20525bb
function viz_p_wins(n::Int)
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel=rich("length of observation phase, t", superscript("*")),
		ylabel=rich("p", subscript("W"), superscript("(t*=$t★)"))
	)

	# computations
	p_wins = [p_win(n, t★) for t★=0:n-1]
	t★_opt = argmax(p_wins)
	p_opt = maximum(p_wins)
	@show t★_opt-1, p_opt
	
	scatterlines!(0:n-1, p_wins)
	lines!(
		[t★_opt-1, t★_opt-1], [0, p_opt], 
		linestyle=:dash, color=colors[1], linewidth=2
	)
	ylims!(0, 0.5)
	save("p_win_$n.pdf", fig)
	fig
end

# ╔═╡ 74d2c390-1ee5-45db-8bee-65e98dbb84e0
viz_p_wins(n)

# ╔═╡ ace6fd4f-9d0a-47c4-b5e2-bbdb4e9866e8
md"# many states"

# ╔═╡ c7630eb9-0e98-428f-a3cc-0080ed5918e5
n_big = 80

# ╔═╡ f6dea05b-53c8-45cb-9e99-840728bb9feb
viz_P(MC_transition_matrix(n_big), "non-policy")

# ╔═╡ 1af4f256-9631-41ae-8241-a28f5dae3400
viz_p_wins(80)

# ╔═╡ Cell order:
# ╠═dbbcbf84-be7d-11ef-3384-2b2a69580c7a
# ╠═be97b1d7-fe22-4112-bb20-19e5bede769b
# ╠═4d03b774-6114-4632-88ed-d698e3a54342
# ╠═aa08c613-192c-4e89-9bab-6b42fe64918f
# ╟─cfffd689-2318-4653-99b8-152c4bc2fdaf
# ╠═d5c07c6f-0bbc-428d-a1b1-dde60d406cfe
# ╠═a0d12780-5ab2-48d0-9854-4fa73a268ae4
# ╠═d4580bec-fe82-4283-981b-86df877ef5a1
# ╠═478ae2a4-c55f-46ad-b597-5045e898fd37
# ╠═eeae5e01-1b04-4855-a889-3ed05a82f529
# ╠═c1da3a78-fcfb-4d68-8c14-461cced050e7
# ╠═2aa85c66-e0bb-4a11-b1e3-b4438bdb53de
# ╠═8c38dcc7-b1f5-4e75-b105-1c61090901b8
# ╠═debeac52-6c15-430a-99a0-5ab7ce018292
# ╠═5b2837f4-dfce-48c4-a9cf-36565767581b
# ╠═9e255ead-d995-4da5-b64d-615003cc55c8
# ╠═c13431c0-1a8a-437d-b947-9b3d8c783fe4
# ╠═a96c1505-5851-4e48-a52b-98a02576a8c3
# ╠═a45403a5-fc90-44c3-aff8-37d5d3b90421
# ╠═2b685bbc-a6f7-4693-98c0-7481e2cfe571
# ╠═b64569d1-3825-4569-a648-f33b278ff002
# ╠═e24e6587-f41d-4f05-bb32-260463c9b4dd
# ╠═9e3251d3-d44d-4df8-ae14-2550167fc3cb
# ╠═a4d206a5-46aa-4527-b713-f939f3f02688
# ╟─c99e89a9-1480-4e35-8621-a5b07cc336fe
# ╠═e98898b5-4437-4fe5-92d6-f1fb035331a8
# ╠═ed545ef2-53d9-4edb-92f2-f672241cf390
# ╠═1e8d615d-25e8-4c2b-be16-4a372670ec0c
# ╠═02f29171-3438-4620-a6f2-8f95624fca7b
# ╠═f612ba1a-a3c3-49c1-96cd-41be7b1cfdfe
# ╠═1eec5172-b76c-4469-b0db-8db3aeac8902
# ╠═1d909b88-9b4a-45cd-b883-6ec8ab57e09f
# ╠═5f7eb3a5-1263-4054-b618-0d9973b41ada
# ╠═cbf80fcc-155b-4cc1-bbd3-8c0024f59fac
# ╠═1d7f9db9-98d7-4be0-b90b-959d1ec9be1e
# ╠═06484fb8-5ddf-4e9f-9178-1fee8aa33add
# ╠═f5982842-2310-47d4-93dc-3a223f439781
# ╠═37753f25-c165-4e5d-9959-5f5238b94a1a
# ╠═4331d0ea-d709-4798-a4a5-01c8b20525bb
# ╠═74d2c390-1ee5-45db-8bee-65e98dbb84e0
# ╟─ace6fd4f-9d0a-47c4-b5e2-bbdb4e9866e8
# ╠═c7630eb9-0e98-428f-a3cc-0080ed5918e5
# ╠═f6dea05b-53c8-45cb-9e99-840728bb9feb
# ╠═1af4f256-9631-41ae-8241-a28f5dae3400
