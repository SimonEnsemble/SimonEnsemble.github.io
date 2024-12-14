### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8b8dce70-81a8-11ef-0864-6f738a464431
begin
	import Pkg; Pkg.activate()
	using  PlutoTeachingTools, DataFrames, CairoMakie, LinearAlgebra, PlutoUI, Random, Distributions, JuMP, Graphs, GraphMakie, MetaGraphsNext, Printf, BenchmarkTools, Test, Ipopt, ColorSchemes, StatsBase, AlgebraOfGraphics

	AlgebraOfGraphics.set_aog_theme!(
			fonts=[AlgebraOfGraphics.firasans("Light"), 
				   AlgebraOfGraphics.firasans("Light")]
		)
	
	update_theme!(
		# fontsize=28, 
		# linewidth=4,
		# markersize=14,
		titlefont=AlgebraOfGraphics.firasans("Light"),
		# size=size
	)
end

# ╔═╡ 6a4f5731-d009-4da2-ac26-229b2840fede
TableOfContents()

# ╔═╡ cdf36e50-0fc6-44a4-a3e3-9e6f61bb9932
md"""
# a Markov chain-based patrol strategy for a surveillance robot

by Cory M. Simon and Paul Morris.

🤖 suppose we wish for a mobile surveillance robot equipped with a [e.g. chemical or radiation] sensor to patrol, i.e. continuously monitor, a set of key locations in an environment. 
we wish for the schedule by which the patrolling robot continually/repeatedly visits locations in the environment to be (i) _efficient_, so that the duration each location goes unvisited tends to be small, and (ii) _unpredictable_ (stochastic), so it cannot be easily exploited by an adversary.
both efficiency and unpredictability in the patrol strategy deter an adversary from attacking in the first place. we explain how to construct such a patrol strategy constituted by a Markov chain. 

🔭 we review _Problem 1_ in [_Duan and Bullo_ Annu. Rev. Control Robot. Auton. Syst. 4(1) (2020)]. 
specifically, we (1) model the environment as a directed graph, (2) allow a Markov chain with this graph's nodes as the state space to constitute the stochastic patrol strategy, then (3) assign the patrol strategy as a Markov chain with (a) a desired stationary distribution and (b) minimal mean hitting time. we write `JuMP.jl` code in Julia to numerically solve the optimization problem on a grid graph, reproducing _Duan and Bullo_'s Fig. 3. 

!!! note \"main references\"
	X Duan and F Bullo. \"Markov Chain–Based Stochastic Strategies for Robotic Surveillance\". _Annual Review of Control, Robotics, and Autonomous Systems_. (2020) [link](https://www.annualreviews.org/content/journals/10.1146/annurev-control-071520-120123)
	
	Robert P. Dobrow. Introduction to Stochastic Processes with R. (2016) [link](https://onlinelibrary.wiley.com/doi/book/10.1002/9781118740712). Chs. 2 and 3.

	Kemeny and Snell. Finite Markov chains. 1960. Chs. 2 and 4.

	Andreas Klappenecker. notes on Markov chains. [link](https://people.engr.tamu.edu/andreas-klappenecker/csce658-s18/markov_chains.pdf).

	Maria José Serna Iglesias. notes on Markov Chains. [link](https://www.cs.upc.edu/~mjserna/docencia/ra-miri/2020/15-RA-MIRI-MC-Stationary.pdf)

	J Hunter. The computation of the mean first passage times for Markov chains. _Linear Algebra and its Applications_ (2018)

## model of the environment

we model the environment as a directed graph $G=(\mathcal{V}, \mathcal{E})$ where $\mathcal{V}=\{1, ..., N\}$ is the set of nodes and $\mathcal{E}$ is the set of arcs. each node represents a key location in the environment. each arc $(u, v)\in\mathcal{E}$ reflects the possibility to directly travel from node $u\in\mathcal{V}$ to node $v\in\mathcal{V}$ (without passing through another node first). self-loops $(u, u) \in \mathcal{E}$ are allowed.

💡 we abstract the robot traversing the environment in 2D/3D continuous space, between key locations, as walking on the graph $G$ from node to node, across arcs. we assume it takes one unit of time for the robot to traverse each arc (alternatively, the arcs could be weighted with travel times; see _Duan and Bullo_). we assume $G$ is _strongly connected_, meaning it's possible for the robot to start at any node $u\in\mathcal{V}$ then take a directed path to any other node $v\in\mathcal{V}$. (if the graph were not strongly connected, we'd need more than one robot to patrol all nodes in $G$ recurrently.)

!!! example
	below, we generate a directed grid graph as a spatial model of some environment.
"""

# ╔═╡ a4d6c9e3-5820-45a7-add7-756bd4c9a918
function generate_grid_graph(; n_row::Int64=3, n_col::Int64=3)
	g = DiGraph(n_row * n_col)
	
	# self-loops
	for u = 1:nv(g)
		add_edge!(g, u, u)
	end

	# construct maps: position (i, j) on grid <---> node ID
	pos_to_node = Dict{Tuple{Int, Int}, Int}()
	v = 0 # node ID
	for i = 1:n_row
		for j = 1:n_col
			v += 1
			pos_to_node[(i, j)] = v
		end
	end
	node_to_pos = Dict(v => pos for (pos, v) in pos_to_node)

	# join neighbors
	for i = 1:n_row
		for j = 1:n_col
			# node in position (i, j)
			u = pos_to_node[(i, j)]
			
			# join up neighbor
			if i != 1
				add_edge!(g, u, pos_to_node[(i - 1, j)])
			end
			# join down neighbor
			if i != n_row
				add_edge!(g, u, pos_to_node[(i + 1, j)])
			end
			# join left neighbor
			if j != 1
				add_edge!(g, u, pos_to_node[(i, j - 1)])
			end
			# join right neighbor
			if j != n_col
				add_edge!(g, u, pos_to_node[(i, j + 1)])
			end
		end
	end

	return g
end

# ╔═╡ af1b3c48-44dd-4464-8721-322a3a47bd77
function viz_g(
	g::DiGraph, 
	P::Union{Nothing, Matrix{Float64}}=nothing;
	v_hl::Union{Nothing, Int}=nothing,
	ed_hl::Tuple{Int, Int}=(-1, -1),
	π::Union{Nothing, Vector{Float64}}=nothing,
	margin = (0.1, 0.1, 0.1, 0.1)
)
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	hidedecorations!(ax)
	hidespines!(ax)

	p_to_label(p::Float64) = p < 1e-8 ? "0" : @sprintf("%.2f", p)

	if ! isnothing(P)
		elabels = [p_to_label(P[ed.src, ed.dst]) for ed in edges(g)]
		edge_color = [
			get(ColorSchemes.Greens, P[ed.src, ed.dst]) 
				for ed in edges(g)
		]
	else
		elabels = nothing
		edge_color = nothing
	end

	edge_width = [(ed.src, ed.dst) == ed_hl ? 8 : 3 for ed in edges(g)]
	arrow_size = [(ed.src, ed.dst) == ed_hl ? 40 : 20 for ed in edges(g)]

	if isnothing(π)
		node_color = ["white" for v = 1:nv(g)]
	else
		node_color = [
			get(ColorSchemes.Greens, π_i, (0.0, 1.0))
				for π_i in π
		]
	end
	if ! isnothing(v_hl)
		node_color[v_hl] = "yellow"
	end
	
	gp = graphplot!(
		g, 
		edge_width=edge_width, 
		node_strokewidth=3,
		arrow_size=arrow_size,
		node_color=node_color,
		node_strokecolor="black",
		elabels=elabels,
		edge_color=edge_color,
		node_size=30,
		ilabels=["$v" for v in vertices(g)],
		ilabels_fontsize=25,
		arrow_shift=:end,
		elabels_fontsize=10,
		curve_distance = -0.2
	)
	if ! isnothing(P)
		Colorbar(
			fig[1, 2], colormap=ColorSchemes.Greens, label="probability"
		)
	end
	# if ! isnothing(π)
	# 	Colorbar(
	# 		fig[1, 3], colormap=ColorSchemes.Blues, label="probability", limits=(0.0, maximum(π))
	# 	)
	# end
	fig
end

# ╔═╡ 4805e7ef-3fd2-43bd-954a-b96998d17e11
g = generate_grid_graph(n_row=3, n_col=3)

# ╔═╡ 8083336f-1e75-4efb-8bb9-48061ac1f303
N = nv(g) # number of states = number of nodes

# ╔═╡ abf83b73-3e76-41ce-8986-08416de867c6
@test is_strongly_connected(g)

# ╔═╡ b5a744db-4d7a-44f0-93d4-b808be4e8f12
viz_g(g)

# ╔═╡ 2ac4a30c-e321-43fb-a43d-f2c216ff6430
md"## a Markov-chain based patrol schedule

💡 we consider a family of stochastic robot patrol strategies on the graph $G$ constituted by a discrete-time, finite, time-homogeneous, irreducible Markov chain with (i) a state space equal to the node set $\mathcal{V}$ and (ii) transition probabilities $p_{ij}=0$ if $(i, j)\notin \mathcal{E}$. the patrol strategy is then specified by the transition matrix $P$ of the Markov chain.

📗 a **Markov chain (MC)** is a sequence of random variables $(X_0, X_1, ...)$ where the random variable $X_k$ represents the state of the robot at time step $k\in\{0, 1, ...\}$. we assign the [discrete] _state space_ of the Markov chain to be the node set $\mathcal{V}$, meaning the realization of each random variable $X_k$ is a node in the graph $G$. (so, the Markov chain is _finite_.) 
if $X_k=x$, then, the surveillance robot resides at node $x\in\mathcal{V}$ at time $k$. in practice, at each time step, the robot samples the next node in the chain to determine the node to visit next; a realization of the MC gives the node visitation schedule.

the distributions of the $X_k$'s follow the **Markov property**:
```math
\begin{align}
	\mathbb{P}[X_{k+1} = x_{k+1} \mid X_0=x_0, ..., X_{k-1}=x_{k-1}, X_k = x_k] = \mathbb{P}[X_{k+1} = x_{k+1} \mid X_k = x_k] \\
	\text{for all } k \geq 0 \text{ and } x_0, ..., x_k \in \mathcal{V}.
\end{align}
```
i.e. if the robot currently resides at some node $i$, the probability the robot next hops to node $j$ depends _only_ on the node $i$ where it currently resides; the history of the walk the robot took to get to its current location $i$ does not matter.

**time-homogeneity** implies the transition probabilities are static over time.
```math
\begin{align}
	\mathbb{P}[X_{k+1} = j \mid X_k = i] = \mathbb{P}[X_{1} = j \mid X_0 = i] \\
	\text{for all } k \geq 0.
\end{align}
```

the properties of a time-homogenous, finite MC are fully characterized by its (1) cardinality-$N$ state space $\mathcal{V}$ and (2) an $N \times N$ (row) stochastic **transition matrix** $P$, whose elements provide the transition probabilities at every time step:
```math
\begin{align}
	 \mathbb{P}[X_{k+1} = j \mid X_k = i] = p_{ij} \text{ for all } k \geq 0 \text{ and } i, j\in \mathcal{V}.
\end{align}
```
since the entries of $P$ are probabilities, $0 \leq p_{ij} \leq 1$ for $i, j\in\mathcal{V}$. more, the sum along each row of $P$ is 1, i.e. $\sum_{j \in \mathcal{V}} p_{ij}=1$ for $i\in\mathcal{V}$, because the robot has to next transition to _some_ node $j \in \mathcal{V}$ when residing in node $i$. in matrix form, we have $P \mathbb{1} = \mathbb{1}$ with $\mathbb{1}$ a vector of ones.

for the MC to be a **feasible** robot patrol strategy, we constrain $p_{ij}=0$ for all $(i, j) \notin \mathcal{E}$ since it's impossible for the robot to directly travel from node $i$ to node $j$ in the real environment if arc $(i, j)$ isn't in the graph $G$ spatially modeling the environment.

a **probability [column] vector** $\pi \in \mathbb{R}^N$ satisfies $0 \leq \pi_i \leq 1$ with $\sum_{i=1}^N \pi_i=1$ and gives a probability distribution over the node set $\mathcal{V}$; $\pi_i$ assigns a probability to node $i\in\mathcal{V}$. 

given a probability distribution over states of the robot at time $k$, $\pi_k$,
the transition matrix $P$ gives the **probability distribution of the robot over next states**, at time $k+1$, $\pi_{k+1}$. 
looking at it a sequence of events: (1) sample the current state, then (2) sample the next state, we have $\pi_{k+1}^\intercal = \pi_k^\intercal P$.  for example, if $\pi_k$ is a one-hot vector with a $1$ in entry $i$, we recover the definition of $p_{ij}$ since $\pi_k^\intercal P$ just picks out row $i$ of $P$. as a consequence of the Markov property and time-homogeneity, $\pi_k^\intercal = \pi_0^\intercal P^k$.

!!! example
	below, we generate a random feasible Markovian robot patrol strategy on our grid graph and visualize its transition matrix $P$ and transition graph.
"

# ╔═╡ a335f814-10ca-44ba-9c94-2cd418a55fb9
function random_P(g::DiGraph)
	@assert is_connected(g)

	# samples from Uniform([0, 1])
	P = rand(nv(g), nv(g))

	# set pᵢⱼ = 0 for (i, j) ∉ edge list
	for i = 1:nv(g)
		for j = 1:nv(g)
			if ! has_edge(g, i, j)
				P[i, j] = 0.0
			end
		end
		# ensure row sum is 1.0
		P[i, :] = P[i, :] / sum(P[i, :])
	end

	return P
end

# ╔═╡ 99e88ed2-3ca3-4705-8114-88e23c35a320
P = random_P(g)

# ╔═╡ e9cff9b8-a99d-4493-b6ad-08d24ed13322
# check each entry is a feasible probability
@test all(P .≤ 1.0) && all(P .≥ 0.0)

# ╔═╡ 2ebba248-d212-4c6b-a74a-6554448e8503
# check each row sums to 1.0
@test all(P * ones(N) .≈ ones(N))

# ╔═╡ 1ddfd1fa-63e2-4d9a-a462-ad0096966710
viz_g(g, P)

# ╔═╡ 5416b864-daea-433d-b7cc-d609a9c7bb64
function viz_heatmap(P::Matrix{Float64}, label::String)
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		# title="Markovian patrol strategy", 
		titlesize=25, 
		xlabel=label == "transition probability" ? 
			rich("state, x", CairoMakie.subscript("k")) : "state i"
		, 
		xlabelsize=25, 
		ylabelsize=25, 
		ylabel=label == "transition probability" ? 
			rich("state, x", CairoMakie.subscript("k+1")) : "state j", 
		xticks=[i for i=1:size(P, 1)], 
		yticks=[i for i=1:size(P, 2)], 
		aspect=DataAspect()
	)

	# hack to make p = 0.0 map to white
	if label == "transition probability"
		colormap = vcat(ColorSchemes.grays[1.0], ColorSchemes.algae[0.0:0.01:1.0])
		hm = heatmap!(P, colorrange=(0, 1), colormap=colormap)
	else
		hm = heatmap!(P)
	end
	Colorbar(fig[1, 2], hm, label=label, labelsize=20)

	return fig
end

# ╔═╡ 1a08aadd-6aa7-4db4-b0ba-e89089857c91
viz_heatmap(P, "transition probability")

# ╔═╡ fcfc675f-ba5a-4391-a71c-6f9be854752d
md"
### communicating states, recurrence, and irreducible MCs

💡 we constrain the MC constituting our robot patrol strategy to be irreducible to make all states recurrent. 

the **transition graph** of the MC is constructed with node set $\mathcal{V}$ (same as the state space of the MC) and edge set $\mathcal{E}^\prime \subseteq \mathcal{E}$ where $(i, j) \in \mathcal{E}^\prime \Leftrightarrow p_{ij}>0$. it describes the possible transitions in the MC.

a state $v\in\mathcal{V}$ is **accessible** from state $u\in\mathcal{V}$, i.e. $u\rightarrow v$, in a MC with transition matrix $P$ iff $(P^n)_{uv}>0$ for some integer $n \geq 0$. i.e., starting in state $u$ there is a nonzero probability of reaching $v$ after a sufficient number ($n$) of steps. $u \rightarrow v$ iff there is a directed path from $u$ to $v$ in the transition graph of the MC.

two states $u\in\mathcal{V}$ and $v\in\mathcal{V}$ **communicate**, i.e. $u \leftrightarrow v$, if $u \rightarrow v$ and $v \rightarrow u$.

an MC is **irreducible** iff all ordered pairs of states $(u, v) \in \mathcal{V} \times \mathcal{V}$ communicate. the transition graph of an MC is strongly connected iff the MC is irreducible. 

within a **communication class** $\mathcal{C}\subseteq \mathcal{V}$, all states communicate with each other (i.e., $u \leftrightarrow v$ for $u,v\in\mathcal{C}$), and no state in the class communicates with any state outside of it (i.e., $\nexists u \in \mathcal{C}, v\in \mathcal{V}\setminus \mathcal{C} : u \leftrightarrow v$). an irreducible Markov chain has exactly one communication class.

a state $v\in\mathcal{V}$ is **recurrent** in an MC if, once the MC begins in state $v$, it will eventually return to state $v$ with certainty. if a state is not recurrent, it is **transient**. beginning in a transient state $u$, there is a nonzero probability that the MC will _never_ return to $u$. 👀 while a recurrent state is visited infinitely many times over the course of the infinite Markov chain, a transient state is visited finite times. 

📔 _all_ states in a finite, irreducible MC are recurrent. so, a robot patrolling according to an irreducible MC continually explores all states in its state space ✔.

a stronger condition than irreducibility is regularity: for a **regular** MC with transition matrix $P$, $P^n>0$ (i.e. all entries of $P^n$ are positive) for some integer $n\geq0$. this tells us it's possible for the MC to be in *any* state after a sufficient number ($n$) of time steps, regardless of the initial state. 
* a regular MC is certainly irreducible, as it's possible to travel from any particular state to any other state. 
* however, an irreducible MC may not be regular, particularly if it is **periodic**. 
the **period** of a _state_ $v\in\mathcal{V}$ is defined as the greatest common divisor of the set $\{ k\geq 1 : (P^k)_{v, v}>0 \}$ giving possible return times.
if the period of a state is one, that state is aperiodic---and, otherwise, periodic. 
a MC is aperiodic if all of its states are aperiodic and, otherwise, periodic. (since the states of a communication class all have the same period, each state in an irreducible MC has the same period.)

note, an irreducible MC with $p_{vv}=0$ for some state $v\in\mathcal{V}$ (giving a self-loop in the transition graph) is necessarily aperiodic.

📔 the finite MC is irreducible and aperiodic $\iff$ the finite MC is regular.

!!! example
	first, we check our randomly-generated MC is irreducible and regular. second, we provide a simple, illustrative example of a non-irreducible MC and an irreducible but non-regular MC.
"

# ╔═╡ a222abfb-70ee-4771-bdb6-2a9266c621c6
# check the Markov chain is irreducible
@test is_strongly_connected(DiGraph(P .≥ 0))

# ╔═╡ acb88a77-c4f8-4310-bab4-f7ff1e38b1ce
# check the Markov chain is regular
@test all(P ^ nv(g) .> 0)

# ╔═╡ a41a8584-c79d-49be-951a-761c5e53e21a
md"**e.g.** a reducible MC.


not irreducible because, while $1\leftrightarrow 2$, $1\nleftrightarrow 3$ and $2\nleftrightarrow 3$. the directed transition graph is not strongly connected, since we can't travel e.g. from node 3 to node 2.
"

# ╔═╡ 897698ab-e557-455a-8482-b281d6039ee2
begin
	# construct directed graph
	local g = DiGraph(3)
	for (u, v) in [(1, 2), (2, 1), (2, 2), (1, 3), (3, 3)]
		add_edge!(g, u, v)
	end

	# assign transition probabilities
	local P = [
		0.0 0.7 0.3; 
		0.8 0.2 0.0;
		0  0  1
	]

	# check for strongly connected transition graph
	if ! is_strongly_connected(DiGraph(P .≥ 0))
		print("transition graph not strongly connected.")
	end
	
	local fig = viz_g(g, P)
	local ax = current_axis(fig)
	fig
end

# ╔═╡ 1564cf3e-13d4-41a2-abcc-44ad3a43c9bb
md"**e.g.** an irreducible but periodic MC.

the directed transition graph is strongly connected, as we can traverse the graph from any node to any other node. 

the set of possible return times to node 2 are {3, 6, 9, 12, 15, ...}, which has a greatest common divisor of 3. therefore state 2 is periodic, and the MC is periodic.
(note, all states have the same period.)
"

# ╔═╡ d50be50a-eac0-4851-97bd-bf65ce96da3b
begin
	# constructed directed graph
	local g = DiGraph(8)
	for (u, v) in [
		(2, 8), (8, 4), (4, 5), (5, 6), (6, 7), (7, 2), (2, 3), (3, 1), (1, 2)
	]
		add_edge!(g, u, v)
	end

	# asign transtion probabilities
	local P = random_P(g)
	
	if is_strongly_connected(g)
		println("transition graph is strongly connected.")
	end

	print("P⁹⁷ > 0?")
	if ! all(P ^ 97 .> 0)
		print(" no.")
	end

	# viz
	local fig = viz_g(g, P)
	local ax = current_axis(fig)
	fig
end

# ╔═╡ 7cd6488b-e829-4e20-b638-e756485188c5
md"**e.g.** an irreducible, aperiodic and thus regular MC.

(a tricky one---I thought this would be periodic...)

the directed transition graph is strongly connected.

the set of possible return times to node 2 are {3, 6, 7, 10, 13, 14, 16, ...}, which has a greatest common divisor of 1. therefore state 2 is aperiodic, and the MC is aperiodic.

so after sufficient steps, the MC can reach _any_ state regardless of the initial state. 
e.g. consider starting at node 2. suppose we wish to return to node 2 in $z$ steps. then we must express: $z=3x+7y$ where integers $x$ and $y$ are the number of times the short and long loop are taken, respectively. easiest place to start is when $z$ is a multiple of three and can be written as $z=3k$ for some integer $k$. then $x:=k$ and $y:=0$ gives us our desired return time $z=3k\geq 3$. now, we need to show that we can increment this return time by 1 and 2. 
```math
z := 3 k + 1 = 3k + 7(1) - 3(2) = 3(k-2) + 7(1).
```
so we achieve a return time of $z=3k+1 \geq 7$ for $x:=k-2$ and $y:=1$.
finally, 
```math
z := 3 k + 2 = 3k + 2[7(1) - 3(2)] = 3(k-4) + 7(2).
```
and we achieve a return time $z=3k+2\geq 14$ for $x:=k-4$ and $y:=2$.
(_credit_: Axel Saenz Rodriguez for this argument.)
"

# ╔═╡ 633c0722-7da1-4326-afd4-38da9df98cfc
begin
	# constructed directed graph
	local g = DiGraph(9)
	for (u, v) in [
		(2, 9), (9, 8), (8, 4), (4, 5), (5, 6), (6, 7), (7, 2), (2, 3), (3, 1), (1, 2)
	]
		add_edge!(g, u, v)
	end

	# asign transtion probabilities
	local P = random_P(g)
	
	if is_strongly_connected(g)
		println("transition graph is strongly connected.")
	end

	@show all(P ^ 97 .> 0)

	# viz
	local fig = viz_g(g, P)
	local ax = current_axis(fig)
	fig
end

# ╔═╡ 58f52271-9e0f-4697-ae39-f1713481b989
md"### stationary distribution of a Markov chain

💡 the **stationary distribution** $\pi^*$ of a finite, discrete-time, homogenous MC with transition matrix $P$ is a special probability distribution over its state space (so, $\pi^*$ is a probability vector), satsifying $\pi^{* \intercal} =\pi^{*\intercal} P$. a _regular_ MC has a unique stationary distribution with positive entries i.e. $\pi^* > 0$.

🤔 **interpreting the stationary distribution of an MC**:

1. if the probability distribution over _current_ states of the robot following the MC is $\pi^*$, then the probability distribution over the _next_ states of the robot is _still_ $\pi^*$. if the chain starts in the stationary distribution, it stays in the stationary distribution. reflecting this, the stationary distribution is sometimes called the _invariant_ distribution.

2. the stationary distribution of a finite, irreducible MC gives the fraction of time the Markov chain (and, thus, robot) spends in each state (node) over the long run: $$\pi_i^* = \lim_{t \rightarrow \infty} \frac{1}{t+1} \sum_{k=0}^t \mathcal{I}[X_k=i]$$. the sum with the indicator function $\mathcal{I}[\cdot]$ counts the number of visits to node $i$ up until time $t$ by a robot following the MC. [Thm. 3.6 of _Dobrow_]

for a regular MC, the stationary distribution is the same as the limiting distribution, i.e.:
```math
\lim_{k \rightarrow \infty} P^k = \mathbb{1}\pi^{*\intercal}
```
the matrix $\mathbb{1}\pi^{*\intercal}$ has the stationary distribution $\pi^{*\intercal}$ in its rows. elementwise, this means
```math
\pi_j^* = \lim_{k \rightarrow \infty} (P^k)_{ij} \text{ for all } i,j \in\mathcal{V}.
```
so regardless of the starting state $i$, the probability the MC lies in state $j$ after a long time is given by the stationary distribution.

💡 the stationary distribution $\pi^*$ associated with the MC constituting a robot patrol strategy is important because it describes the long-run frequency that the robot visits each node in the graph. for robot patrol, we wish to _specify_ the stationary distribution of our robot to ensure e.g. (a) each node is visited with the same frequency over the long-run or, alternatively, (b) a subset of important, high-value locations are visisted _more_ frequently over the long-run than less important locations. 
so, here we make the stationary distribution an _input_ to the algorithm that finds the optimal stochastic patrol strategy. i.e. we _constrain_ the MC patrol strategy to have the stationary distribution $\pi^*$ _we_ specify.

💻 how can we **compute the stationary distribution** of a Markov chain with transition matrix $P$? well, we have a system of equations for the unknown $\pi^*$:
```math
\pi^{*\intercal}=\pi^{*\intercal}P \implies (P^\intercal-I)\pi^*=\mathbb{0},
```
with one of the equations redundant owing to each row of $P$ summing to one, but compensated for by the constraint that $\pi^*$ is a probability vector:
```math
\sum_{i=1}^N \pi^*_i=1 \implies \mathbb{1}^\intercal \pi^*=1.
```
we write this system of equations as a linear matrix equation
```math
\begin{bmatrix}
P^\intercal-I \\ \mathbb{1}^\intercal
\end{bmatrix} \pi^* = 
\begin{bmatrix}
\mathbb{0} \\ 1
\end{bmatrix}
```
then numerically solve for $\pi^*$ (via Gaussian elimination).

!!! example
	we solve for and visualize the stationary distribution of our random Markovian patrol strategy below.
"

# ╔═╡ 286e3d55-c127-43a5-94ca-eb58afbe4e63
function stationary_distn(P::Matrix{Float64})
	N = size(P)[1]
	π★ = vcat(P' - I(N), ones(N)') \ vcat(zeros(N), 1)
	return π★
end

# ╔═╡ 9abd6864-3915-49b8-ae4e-4487e0454fdb
π★ = stationary_distn(P)

# ╔═╡ adb9b1e6-ff7b-4eb9-abbd-af432fee03dd
@test all(π★ .≥ 0) && all(π★ .≤ 1.0)

# ╔═╡ dd502b60-2e9e-4d63-9a04-76a7b80dbeb3
@test sum(π★) .≈ 1.0

# ╔═╡ 6173f7df-921b-4ab9-89ad-f682cf49718f
@test π★' * P ≈ π★'

# ╔═╡ b5cdc231-3896-4e54-aba6-41a143b37df2
begin
	fig = Figure()
	ax = Axis(
		fig[1, 1], xlabel="node", ylabel="probability", xticks=1:N, 
		title="stationary dist'n"
	)
	barplot!(1:N, π★, color="green")
	ylims!(0, nothing)
	fig
end

# ╔═╡ 56ab7f79-dbd5-45f4-8b7b-a3f15c107e5a
viz_g(g, P, π=π★)

# ╔═╡ c5c5cb68-2c8b-4a63-83c6-b2f1cf45de39
# test that stationary dist'n = the limiting distribution
@test P ^ 2500 ≈ ones(N) * π★'

# ╔═╡ 48c4e067-bd23-4a72-a301-9738976d994e
md"
### mean first hitting time from state $i$ to state $j$ 

💬 Markov chain slang for \"visit state $i$\" is \"hit up state $i$\".

📗 the MC's **first hitting time from state** $i$ **to state** $j$ is a random variable representing the time, measured by the number of transitions i.e. robot arc-hops, for a robot starting at node $i$ to first reach node $j$:
```math
T_{ij} := \min \{ k \geq 1 : X_k=j  \mid X_0=i\} \text{ for } i, j \in \mathcal{V}.
```
for the special case $i=j$, $T_{ii}$ is the **recurrence time of state** $i$---the time for the robot starting at node $i$ to first return to node $i$.

for a finite, irreducible MC, the hitting time $T_{ij}$ is finite for all pairs of states $i,j\in\mathcal{V}$. 

via the **strong Markov property**, the MC **regenerates** itself (starts anew) once it travels from some state $i^\prime$ to another state $i$. so, we could also write:
```math
T_{ij} = \min \{ k \geq 1 : X_{s+k}=j  \mid X_{s}=i\} \text{ for } i, j \in \mathcal{V}
```
where the MC starts at state $i$ at an arbitrary time $s\geq0$.

the **expected first hitting time from state $i$ to state $j$**, $\mathbb{E}[T_{ij}]$, is the expected number of steps it takes for the robot to first visit node $j$ when starting at node $i$. and $\mathbb{E}[T_{ii}]$ is the **mean recurrence time of state $i$**. with $m_{ij} := \mathbb{E}[T_{ij}]$, we store the mean first hitting times between each ordered pair of states in the $N\times N$ **mean first hitting time matrix** $M$ whose entry $(i, j)$ is $m_{ij}$. 

#### a relation for $M$

the first hitting time matrix $M$ follows from the transition matrix $P$ of the MC.

we can write equations relating the $m_{ij}$'s through a thought experiment. suppose at time $k=0$, the robot resides at state $X_0=i$ and we're interested in $m_{ij}$.
next, either (1) the robot hops across arc $(i, j)$ to our node $j$ of interest, i.e. $X_1=j$, with probability $p_{ij}$, or (2) the robot hops to another node $\ell \neq j$, i.e. $X_1=\ell$, with probability $p_{i\ell}$ each. 
in the former case, the realized first hitting time from $i$ to $j$ is one. in the latter case, the Markov chain has regenerated itself, so the hitting time from $i$ to $j$ is one (the step to get to $\ell$) plus the random variable $T_{\ell j}$ giving the hitting time from the intermediate state $\ell$ to the state of interest $j$. consequently, for all $(i, j)\in\mathcal{E}$:
```math
\begin{align}
m_{ij} &= p_{ij} 1 + \sum_{\substack{\ell=1 \\ \ell \neq j}}^N p_{i\ell} (1 + m_{\ell j}) \\
       &= p_{ij}  + \sum_{\substack{\ell=1 \\ \ell \neq j}}^N p_{i\ell} + \sum_{\substack{\ell=1 \\ \ell \neq j}}^N p_{i\ell} m_{\ell j} \\
	&= 1 + \sum_{\ell=1}^N p_{i\ell} m_{\ell j} - p_{ij} m_{jj}
\end{align}
```
giving coupled equations for the $m_{ij}$'s. the massaging in the last two lines is to help show we can write these equations in matrix form as
```math
M= \mathbb{1}\mathbb{1}^\intercal + P(M - \text{diag}(M)).
```
the diagonal matrix $\text{diag}(M)$ has $m_{ii}$'s down its diagonal.

#### the mean recurrence times

if we multiply the above eqn. on the _left_ by the stationary distribution $\pi^{*\intercal}$ we get:
```math
\begin{align}
\pi^{* \intercal} M & = \pi^{* \intercal}\mathbb{1}\mathbb{1}^\intercal + \pi^{* \intercal}P(M - \text{diag}(M)) \\
\pi^{* \intercal} M & = \mathbb{1}^\intercal + \pi^{* \intercal}(M - \text{diag}(M)) \\
\pi^{* \intercal} \text{diag}(M) & = \mathbb{1}^\intercal \\
& \implies m_{ii} = 1/\pi^*_i \text{ for } i \in \{1, ..., N\}.
\end{align}
```
where the first simplification owes to $\pi^*$ being a probability vector and the definition of the stationary distribution.
so, the mean recurrence time of a state $i$ is the inverse of the probability of that state in the stationary distribution of the MC. makes total sense, as the long-run frequency of visiting state $i$ is $\pi^*_i$, so the mean recurrence time of state $i$ is the inverse of that frequency. 

an alternative derivation of $m_{ii}=1/\pi^*_i$ is provided in Theorem 3.6 of Dobrow.

#### computing the mean first hitting time matrix

key to computing the mean first hitting time matrix $M$ of a regular MC is the **fundamental matrix of a regular MC** with transition matrix $P$:
```math
Z:=(I-(P-\mathbb{1}\pi^{*\intercal}))^{-1},
```
which exists and, from the geometric series for matrices, the binomial theorem for matrices, and the fact $(\mathbb{1}\pi^{*\intercal})^k=\mathbb{1}\pi^{*\intercal}$, is equal to:
```math
Z=I - \sum_{k=1}^\infty (P^k-\mathbb{1}\pi^{*\intercal}).
```
see Theorem 4.3.1 of Kemeny and Snell. vaguely, we see $Z$ characterizes how fast the regular MC converges to its limiting distribution, since we're looking at a difference between $P^k$ and $\mathbb{1}\pi^{*\intercal}$.

the mean first hitting time matrix can be computed from the fundamental matrix:
```math
M = (I-Z+\mathbb{1}\mathbb{1}^\intercal \text{diag}(Z)) \text{diagm}([1/\pi_1^*, ..., 1/\pi_N^*]).
```
one can verify by showing this expression satisfies our relation for $M$ above. see Theorem 4.4.7 of Kemeny and Snell.

!!! example
	we compute and visualize the mean first hitting times of our random MC below.
"

# ╔═╡ d64ec721-e50c-4cc7-9320-6035519b6f88
function fundamental_matrix(
	# transition matrix
	P::Matrix{Float64}, 
	# stationary distribution
	π★::Union{Vector{Float64}, Nothing}=nothing
)
	N = size(P)[1] # number of states
	
	# compute stationary distribution if we haven't already
	if isnothing(π★)
		π★ = stationary_distn(P)
	end
	
	# build the fundamental matrix
	Z = inv(I(N) - (P - ones(N) * π★'))
	
	return Z
end

# ╔═╡ be9b316d-f8e6-4e18-91ef-0fb2f9c255fb
# the fundamental matrix
Z = fundamental_matrix(P)

# ╔═╡ 26b8bbd6-3010-4d25-bace-202fd05fda45
function mean_first_hitting_time_matrix(P::Matrix{Float64})
	N = size(P)[1]
	π★ = stationary_distn(P)
	Z = fundamental_matrix(P, π★)
	M = (I(N) - Z + ones(N, N) * diagm(diag(Z))) * diagm(1 ./ π★)
	return M
end

# ╔═╡ a00a16e7-b3e0-4f59-bc67-c3200370c886
M = mean_first_hitting_time_matrix(P)

# ╔═╡ dea567a7-28a0-4107-a318-43af2796ddbf
begin
	local fig = hist(M[:], axis=(;xlabel="mean first hitting time mᵢⱼ", ylabel="# pairs of states"))
	xlims!(0, nothing)
	ylims!(0, nothing)
	fig
end

# ╔═╡ 7fec50e2-e20b-45f0-b7ab-a40382dfad63
viz_heatmap(M, "mean first hitting time")

# ╔═╡ 09181a10-34b0-4191-8129-65d008d3321b
# test that our implicit equation for M is satisfied.
@test all(M .≈ ones(N, N) + P * (M - diagm(diag(M))))

# ╔═╡ 36847cdd-64e8-4d06-8e77-6355336b5142
md"### the Kemeny constant of an MC

💡 a key quantity for scoring the efficiency of a stochastic robot patrol strategy based on a regular MC is the Kemeny constant of that MC.

the **Kemeny constant**, or **mean first hitting time**, of the finite, homogenous, regular MC with transition matrix $P$ is:
```math
	\mathcal{M}(P):= \sum_{j=1}^n \pi_j^* m_{i j}
```
for whatever state $i \in \mathcal{V}$ we arbitrarily choose since the quantity is the same for all states. 

❗ though not explicitly indicated here, keep in mind that _both_ the stationary distribution $\pi^*$ and the mean first hitting times in $M$ are a function of the transition matrix $P$.

to see that $\mathcal{M}(P)$ is independent of $i$: 
use $m_{ii}=1/\pi^*_i$ then multiply our implicit eqn. for $M$ on the _right_ by $\pi^*$, giving:
```math
\begin{align}
M \pi^* & = \mathbb{1}\mathbb{1}^\intercal \pi^* + P(M - \text{diagm}([1/\pi_1^*, ..., 1/\pi_N^*])) \pi^* \\
M \pi^* & = \mathbb{1} + PM \pi^* - P \mathbb{1} \\
M \pi^* & = PM \pi^* \\
& \implies M\pi^* \text{ is an eigenvector of } P \text{ with eigenvalue } 1.
\end{align}
```
now, we already know $P\mathbb{1}=\mathbb{1}$ since each row of $P$ sums to one. Lemma 3.15 of Dobrow tells us that $\mathbb{1}$ is a unique eigenvector (up to scaling) of $P$ as a consequence of $P$ being regular. hence, we know $M\pi^*=c\mathbb{1}$ for some constant $c$. this useful because it tells us the rows of $M\pi^*$ are equal [to the constant $c$, whatever that is...]:
```math
m_{i, 1}\pi_1^* + \cdots + m_{i, N} \pi_N^* = m_{i^\prime, 1}\pi_1^* + \cdots + m_{i^\prime, N} \pi_N^* = c \text{ for } i, i^\prime \in \mathcal{V} .
```
therefore, the quantity
```math
\sum_{j=1}^n \pi_j^* m_{i j}
```
is the same for all $i$, so we call it the _Kemeny constant_ of the [entire---not just a state of the] MC.

🤔 **interpreting the Kemeny constant**:

multiplying $\mathcal{M}(P)$ by one in a fancy way:
```math
	\mathcal{M}(P)= \left(\sum_{i=1}^n \pi_i^* \right) \sum_{j=1}^n \pi_j^* m_{i, j}= \sum_{i=1}^n \sum_{j=1}^n \pi_i^* \pi_j^* m_{ij} = \pi^{*\intercal} M \pi^*.
```
shows that the Kemeny constant is the expected first hitting time from some state sampled from the stationary distribution to another state sampled from the stationary distribution.

another interpretation of the Kemeny constant comes from [Levene & Loizou. \"Kemeny's Constant and the Random Surfer.\" _The American Mathematical Monthly_. (2002)], which I paraphrase: suppose the robot following the MC's transition matrix initially had a mission to reach a destination node $j\in\mathcal{V}$, but (i) forgot its destination and (ii) is currently lost (i.e. doesn't know its current state). Kemeny's constant $\mathcal{M}(P)$ is the expected number of hops the robot must take to reach its forgotten destination. 

**the importance of the Kemeny constant for MC-based robot patrol strategies**:

💡  _Duan and F Bullo_ state, the Kemeny constant \"$\mathcal{M}(P)$ measures the expected time it takes for the surveillance agent to transition between two locations randomly selected according to the stationary distribution\". 
we wish for the mean hitting time to be small, so that the robot doesn't leave a big gap of time between its visits between pairs of nodes, giving efficiency in its patrol schedule.

#### computing the Kemeny constant

we can compute the Kemeny constant by computing the mean first hitting time matrix $M$ as above, then using $\pi^{*\intercal}M\pi^*$. 

more efficiently, we can write:

```math
\begin{align}
\mathcal{M}(P) & = \pi^{*\intercal}M\pi^* \\
 & = \pi^{*\intercal}(I-Z+\mathbb{1}\mathbb{1}^\intercal \text{diag}(Z)) \text{diagm}([1/\pi_1^*, ..., 1/\pi_N^*]) \pi^* \\
 & = \pi^{*\intercal}(I-Z+\mathbb{1}\mathbb{1}^\intercal \text{diag}(Z)) \mathbb{1} \\
 & = \pi^{*\intercal}I \mathbb{1} -\pi^{*\intercal} Z\mathbb{1} +\pi^{*\intercal}\mathbb{1}\mathbb{1}^\intercal \text{diag}(Z)\mathbb{1} \\
 & = 1 -\pi^{*\intercal} Z\mathbb{1} +\text{tr}(Z) \\
\end{align}
```
where the trace of the fundamental matrix $Z$ pops out.
now,
```math
\begin{align}
\pi^{*\intercal} Z\mathbb{1} & = \pi^{*\intercal} \left(I - \sum_{k=1}^\infty (P^k-\mathbb{1}\pi^{*\intercal})\right) \mathbb{1} \\
& = \pi^{*\intercal} I \mathbb{1} - \sum_{k=1}^\infty \pi^{*\intercal} (P^k-\mathbb{1}\pi^{*\intercal})\mathbb{1} \\
& = 1 - \sum_{k=1}^\infty (\pi^{*\intercal} P^k\mathbb{1}-\pi^{*\intercal} \mathbb{1}\pi^{*\intercal}\mathbb{1}) \\
& = 1 - \sum_{k=1}^\infty (\pi^{*\intercal} \mathbb{1}-\pi^{*\intercal} \mathbb{1}) \\
& = 1.
\end{align}
```

combining the last two equations, we have:
```math
\mathcal{M}(P)=\text{tr}(Z).
```
beautiful! 🌄

!!! example
	we compute the Kemeny constant of our random MC below---in two different ways.
"

# ╔═╡ 44062100-2699-4727-a594-4c6b100f163f
# Kemeny constant. method 1
ℳ = π★' * M * π★

# ╔═╡ dd7ea312-c369-48fc-a54e-7b5f53b11d36
# Kemeny constant. method 2
tr(Z)

# ╔═╡ 03e86033-187a-49dd-a695-4ba076b4006d
md"""
## the optimization problem to give the MC-based patrol schedule

💡 
given a directed graph $G=(\mathcal{V}, \mathcal{E})$ modeling our environment, our design principle for a Markov chain serving as a stochastic patrol strategy for a robot is to specify its stationary distribution $\pi^*$ then tune the transition matrix $P$ (under arc-hopping feasibility constraints) to minimize the mean hitting time/ Kemeny constant of the MC. our desire is for the time it takes for the robot to travel between states to tend to be small.

the optimization problem we must solve is then the nonlinear program:

```math
\begin{align}
	&\min_{P \in \mathbb{R}^{N\times N}} \mathcal{M}(P) \\
	&p_{ij} \geq 0 \text{ for } (i, j) \in \mathcal{E}\\
	&p_{ij} = 0 \text{ for } (i, j) \notin \mathcal{E} \\
	&P \mathbb{1} = \mathbb{1} \\
	&\pi^{* \intercal} P = \pi^* 
\end{align}
```

the constraints ensure, respectively: (1) the transition probabilities are non-negative for arcs present in the directed graph modeling the environment; (2) the transition probabilities are zero for arcs not present in the directed graph modeling the environment; (3) each row of the transition matrix sums to 1; and (4) the chain has the desired stationary distribution. (note, $p_{ij}\leq 1$ is omitted because it is implied by constraints (1), (2), and (3).)

💻 below is our Julia code in `JuMP.jl` to solve this nonlinear program. since the problem is not convex, we solve it many times over different initial guesses for the transition matrix, then return the best solution we find.
"""

# ╔═╡ b29f9ee7-3f11-44b5-a95d-4784624ab827
function opt_patrol_strategy(
	# graph model of the environment
	g::DiGraph, 
	# desired stationary distribution
	π★::Vector{Float64}; 
	# initial guess for optimal transition matrix
	P₀::Union{Nothing, Matrix{Float64}}=nothing,
	JuMP_verbose::Bool=true
)
	n = nv(g) # size of MC's state space = number of nodes in graph
	@assert n == length(π★)

	model = Model(Ipopt.Optimizer)
	if ! JuMP_verbose
		set_silent(model)
	end

	#=
    decision variables
        (the transition matrix)
    =#
    @variable(model, P[1:n, 1:n] .≥ 0)

    #=
    constraints
    =#
    # if edge (i, j) not in the graph, then pᵢⱼ = 0.
    for u = 1:n
        for v = 1:n
            if ! has_edge(g, u, v)
                @constraint(model, P[u, v] == 0)
            end
        end
    end

    # each row of stochastic matrix P sums to one.
    @constraint(model, P * ones(n) == ones(n))

    # the desired stationary distribution is satisfied.
    @constraint(model, P' * π★ == π★)

    #=
    objective function
        (minimize Kemeney constant)
		see https://github.com/jump-dev/JuMP.jl/issues/2060 
			for why it's written this awkward way...
    =#
    function kemeny_constant(P_flat...)
        # rebuild transition matrix
        P = reshape(collect(P_flat), n, n)

        # fundamenatal matrix
        Z = inv(I(n) - (P - ones(n) * π★'))

        # Kemeny constant
        ℳ = tr(Z)

        return ℳ
    end

    @operator(model, ℳ, n ^ 2, kemeny_constant)
    @objective(model, Min, ℳ(P...))

    #=
    initialize solution
    =#
	if ! isnothing(P₀)
		set_start_value.(P, P₀)
	end
	
	#=
    solve the optimization problem
       --> extract optimizer and optimum.
    =#
    JuMP.optimize!(model)

    ℳ_opt = objective_value(model)

    P_opt = zeros(n, n)
    for i = 1:n
        for j = 1:n
            P_opt[i, j] = value.(P)[i, j].value
        end
    end

	# TODO return the objective too.
	return P_opt, ℳ_opt
end

# ╔═╡ f8a98fff-c19d-480e-8146-66a8c9925382
function opt_patrol_strategy(g::DiGraph, π★::Vector{Float64}, n_runs::Int64)
	# initialize mean hitting time and optimum strategy
	ℳ_opt = Inf
	P_opt = zeros(length(π), length(π))

	for r = 1:n_runs
		# choose random starting transition matrix
		P₀ = random_P(g)
		
		# solve optimization problem
		P, ℳ = opt_patrol_strategy(g, π★, P₀=P₀, JuMP_verbose=false)
		
		# if improved sol'n found, store it!
		if ℳ < ℳ_opt
			ℳ_opt = ℳ
			P_opt = P
		end
	end

	# return best sol'n
	return P_opt, ℳ_opt
end

# ╔═╡ 54d82a51-ae9b-4075-a965-9c15ac2fc1a9
md"
!!! example
	we design an optimal patrol strategy for our grid graph, supposing the desired stationary distribution is uniform.

	🚀
"

# ╔═╡ 12c187c6-9ca8-4b0a-8d01-0380772a6a12
md"our desired stationary distribution, suppose, is uniform. we want the robot to visit each node with an equal frequency over the long-run."

# ╔═╡ 62891a43-1705-433c-b511-8f200ac436e4
desired_π★ = [1/nv(g) for i = 1:nv(g)]

# ╔═╡ 1e25852b-c158-44c6-bfe1-7213e9467eb7
md"solving the nonlinear program for the transition matrix giving the minimal Kemeny constant, under the constraint of our desired stationary distribution..."

# ╔═╡ 214b84ae-b0b1-47ed-8e33-8efdcb4433ca
P_opt, ℳ_opt = opt_patrol_strategy(g, desired_π★, 25)

# ╔═╡ ea1616dd-cb46-4269-831a-2adaa57b25d3
print("minimal Kemeny constant: ", round(ℳ_opt, digits=2))

# ╔═╡ 338581b4-6cf1-4359-9d16-1e8e68c37114
# got the desired stationary distribution?
@test all(desired_π★ .≈ stationary_distn(P_opt))

# ╔═╡ eb9bf19f-b975-4821-aa16-cc5bfe71f231
# a proper transition matrix?
@test all(P_opt * ones(N) .≈ ones(N))

# ╔═╡ 0dc9800b-95b5-460a-b01e-33ce99216969
md"looking at the mean hitting times for the optimal patrol strategy.

🚀 woah, much smaller than our randomly chosen MC!"

# ╔═╡ e5528470-640a-4811-ba19-2d3cee199156
M_opt = mean_first_hitting_time_matrix(P_opt)

# ╔═╡ 7837a664-85f9-4df1-ae5d-efbce50ce2fd
# Kemeny constant consistent w transition matrix and computed hitting times?
@test ℳ_opt ≈ desired_π★' * M_opt * desired_π★

# ╔═╡ 681201db-dbe2-4a50-8da4-c6cf910e7c83
begin
	local fig = hist(M_opt[:], 
		axis=(;xlabel="mean first hitting time mᵢⱼ", ylabel="# pairs of states")
	)
	vlines!(ℳ_opt, color="black")
	xlims!(0, nothing)
	ylims!(0, nothing)
	fig
end

# ╔═╡ 7cc8566b-7f03-40d9-8d2e-461eb15612df
viz_heatmap(M_opt, "mean first hitting time")

# ╔═╡ df16675c-bb14-48b8-92ba-766dd50ac37b
# is the MC regular?
@test all(P_opt ^ (2*N) .> 0.)

# ╔═╡ 8fb9f5bf-727b-4fee-9c0f-fcd613124674
viz_g(g, P_opt, π=desired_π★)

# ╔═╡ aca9d535-c3d3-4173-9111-5fac35dfe364
md"### 🤖 simulating the MC patrol schedule


equipped with an optimal MC-based patrol strategy on our graph $G$ having our desired stationary distribution---described by the transition matrix $P_{opt}$ giving the minimal Kemeny constant 😎---the schedule by which we instruct the robot to continually visit locations in the environment is specified by a realization of this MC. 
below, we sample a realization of the Markov chain and make a video of the robot patrolling the environment. very practical!
"

# ╔═╡ 2efad111-9f3e-42fa-ba38-5a06763faa49
# given initial state x₀, sample MC described by transition matrix P for T steps.
function sample_markov_chain(x₀::Int, P::Matrix{Float64}, T::Int)
	# number of states
	n = size(P, 1)
	
	# store realization of the MC here.
	x = zeros(Int, T)
	x[1] = x₀ # initial state
	for t = 2:T
		x[t] = sample(
			# candidate next states
			1:n, 
			# probabilities of each next state
			ProbabilityWeights(P[x[t-1], :][:])
		)
	end
	# return realization of MC
	return x
end

# ╔═╡ 84559da2-ca22-4719-b209-1a7b98b6b120
x = sample_markov_chain(1, P_opt, 1000)

# ╔═╡ e3f75f0b-4c17-40d5-91d5-421c3d02fee6
@bind clock_time Clock(0.25, true)

# ╔═╡ 63785e7d-a244-4e67-9d00-fc393581b1a1
mc_time = ceil(Int, clock_time / 2)

# ╔═╡ 94a86a38-6c58-449c-b3e6-626b7f1a21fb
if clock_time % 2 == 1
	viz_g(g, P_opt, v_hl=x[mc_time])
else
	viz_g(g, P_opt, v_hl=x[mc_time], ed_hl=(x[mc_time], x[mc_time+1]))
end

# ╔═╡ Cell order:
# ╟─8b8dce70-81a8-11ef-0864-6f738a464431
# ╟─6a4f5731-d009-4da2-ac26-229b2840fede
# ╟─cdf36e50-0fc6-44a4-a3e3-9e6f61bb9932
# ╟─a4d6c9e3-5820-45a7-add7-756bd4c9a918
# ╟─af1b3c48-44dd-4464-8721-322a3a47bd77
# ╠═4805e7ef-3fd2-43bd-954a-b96998d17e11
# ╠═8083336f-1e75-4efb-8bb9-48061ac1f303
# ╠═abf83b73-3e76-41ce-8986-08416de867c6
# ╟─b5a744db-4d7a-44f0-93d4-b808be4e8f12
# ╟─2ac4a30c-e321-43fb-a43d-f2c216ff6430
# ╠═a335f814-10ca-44ba-9c94-2cd418a55fb9
# ╠═99e88ed2-3ca3-4705-8114-88e23c35a320
# ╟─1a08aadd-6aa7-4db4-b0ba-e89089857c91
# ╠═e9cff9b8-a99d-4493-b6ad-08d24ed13322
# ╠═2ebba248-d212-4c6b-a74a-6554448e8503
# ╟─1ddfd1fa-63e2-4d9a-a462-ad0096966710
# ╟─5416b864-daea-433d-b7cc-d609a9c7bb64
# ╟─fcfc675f-ba5a-4391-a71c-6f9be854752d
# ╠═a222abfb-70ee-4771-bdb6-2a9266c621c6
# ╠═acb88a77-c4f8-4310-bab4-f7ff1e38b1ce
# ╟─a41a8584-c79d-49be-951a-761c5e53e21a
# ╟─897698ab-e557-455a-8482-b281d6039ee2
# ╟─1564cf3e-13d4-41a2-abcc-44ad3a43c9bb
# ╟─d50be50a-eac0-4851-97bd-bf65ce96da3b
# ╟─7cd6488b-e829-4e20-b638-e756485188c5
# ╟─633c0722-7da1-4326-afd4-38da9df98cfc
# ╟─58f52271-9e0f-4697-ae39-f1713481b989
# ╠═286e3d55-c127-43a5-94ca-eb58afbe4e63
# ╠═9abd6864-3915-49b8-ae4e-4487e0454fdb
# ╠═adb9b1e6-ff7b-4eb9-abbd-af432fee03dd
# ╠═dd502b60-2e9e-4d63-9a04-76a7b80dbeb3
# ╠═6173f7df-921b-4ab9-89ad-f682cf49718f
# ╟─b5cdc231-3896-4e54-aba6-41a143b37df2
# ╟─56ab7f79-dbd5-45f4-8b7b-a3f15c107e5a
# ╠═c5c5cb68-2c8b-4a63-83c6-b2f1cf45de39
# ╟─48c4e067-bd23-4a72-a301-9738976d994e
# ╠═d64ec721-e50c-4cc7-9320-6035519b6f88
# ╠═be9b316d-f8e6-4e18-91ef-0fb2f9c255fb
# ╠═26b8bbd6-3010-4d25-bace-202fd05fda45
# ╠═a00a16e7-b3e0-4f59-bc67-c3200370c886
# ╟─dea567a7-28a0-4107-a318-43af2796ddbf
# ╟─7fec50e2-e20b-45f0-b7ab-a40382dfad63
# ╠═09181a10-34b0-4191-8129-65d008d3321b
# ╟─36847cdd-64e8-4d06-8e77-6355336b5142
# ╠═44062100-2699-4727-a594-4c6b100f163f
# ╠═dd7ea312-c369-48fc-a54e-7b5f53b11d36
# ╟─03e86033-187a-49dd-a695-4ba076b4006d
# ╠═b29f9ee7-3f11-44b5-a95d-4784624ab827
# ╠═f8a98fff-c19d-480e-8146-66a8c9925382
# ╟─54d82a51-ae9b-4075-a965-9c15ac2fc1a9
# ╟─12c187c6-9ca8-4b0a-8d01-0380772a6a12
# ╠═62891a43-1705-433c-b511-8f200ac436e4
# ╟─1e25852b-c158-44c6-bfe1-7213e9467eb7
# ╠═214b84ae-b0b1-47ed-8e33-8efdcb4433ca
# ╟─ea1616dd-cb46-4269-831a-2adaa57b25d3
# ╠═338581b4-6cf1-4359-9d16-1e8e68c37114
# ╠═eb9bf19f-b975-4821-aa16-cc5bfe71f231
# ╟─0dc9800b-95b5-460a-b01e-33ce99216969
# ╠═e5528470-640a-4811-ba19-2d3cee199156
# ╠═7837a664-85f9-4df1-ae5d-efbce50ce2fd
# ╟─681201db-dbe2-4a50-8da4-c6cf910e7c83
# ╟─7cc8566b-7f03-40d9-8d2e-461eb15612df
# ╠═df16675c-bb14-48b8-92ba-766dd50ac37b
# ╟─8fb9f5bf-727b-4fee-9c0f-fcd613124674
# ╟─aca9d535-c3d3-4173-9111-5fac35dfe364
# ╠═2efad111-9f3e-42fa-ba38-5a06763faa49
# ╠═84559da2-ca22-4719-b209-1a7b98b6b120
# ╟─e3f75f0b-4c17-40d5-91d5-421c3d02fee6
# ╟─63785e7d-a244-4e67-9d00-fc393581b1a1
# ╟─94a86a38-6c58-449c-b3e6-626b7f1a21fb
