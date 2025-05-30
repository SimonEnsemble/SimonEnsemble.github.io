### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 7a733b92-35a6-11f0-0323-f9de5de91771
begin
    import Pkg; Pkg.activate()
    using CairoMakie, PlutoUI, JuMP, DataFrames
    import HiGHS

    update_theme!(fontsize=18, linewidth=3)
end

# ╔═╡ 83d420f4-3c94-467f-8c5a-879b6201e74c
md"# solving a linear program for optimal wine blending in Julia

!!! note \"reference\"
	this is a modified version of the fun, clever, and practical problem presented in:
	> D. Putler. \"Getting to Know Optimization: Linear Programming\". Alteryx Blog. 2020. [link](https://community.alteryx.com/t5/Data-Science/Getting-to-Know-Optimization-Linear-Programming/ba-p/513793)
"

# ╔═╡ 837796eb-0f64-4e3d-8e45-fdf5ff3a358b
html"<img src=\"https://raw.githubusercontent.com/SimonEnsemble/CHE361_W2025/refs/heads/main/drawings/wine_blending.png\" width=400>"

# ╔═╡ 600116f7-21c0-463f-a7ce-8edcace8542f
md"
🍷 suppose we own a winery and have produced two large vats of wine---one, pure Syrah; the other, pure Grenache. the following table presents data regarding our cost to produce each variety of wine, the alcohol and acid content of each variety, and our inventory of each variety.

| wine variety | cost   | alcohol content   | acid content | inventory |
|--------------|--------|-------------------|--------------|--------------|
| Syrah        | \$4.5/L | 11% by volume     | 0.38 g/L     | 3000 L |
| Grenache     | \$6.5/L | 13% by volume     | 0.31 g/L     | 4500 L |

a customer will pay us \$8/L for up to 5,000 L of a Syrah/Grenache blend meeting two specifications relating to perceived quality. first, the blend should contain at least 11.5% alcohol by volume. second, the blend should contain less than 0.36 g/L acid.

to maximize our profit, what volume of Grenache and what volume of Syrah should we blend together for the customer?

💻 let's formulate and solve this optimization problem---a linear program---for our optimal wine blend using [`JuMP.jl`](https://jump.dev/JuMP.jl/stable/) in Julia. 

💡 in this optimization problem, the:
* _decision variables_ are the volume of Syrah and volume of Grenache to blend together.
* _objective function_ [of the decision variables] we wish to maximize is profit.
* [inequality] _constraints_ specify a minimal alcohol and maximal acid content of the wine blend, a maximal total volume of wine blend the customer will purchase from us, and a maximal volume of each variety the blend can contain based on our inventory. (not to mention, the non-negative constraints that are obvious to us, but not to a computer.)
"

# ╔═╡ f030c170-4f8c-4e70-86b1-0944935d550b
md"🍷 first, we put the data regarding the wine into a dictionary for convenience."

# ╔═╡ 711f0fb2-3a27-470e-ba4a-66b9f0ead816
data = Dict(
	"cost"      => Dict("Grenache" => 6.5,    "Syrah" => 4.5),    # $/L
	"acid"      => Dict("Grenache" => 0.31,   "Syrah" => 0.38),   # g/L
	"alcohol"   => Dict("Grenache" => 13.0,   "Syrah" => 11.0),   # vol %
	"inventory" => Dict("Grenache" => 4500.0, "Syrah" => 3000.0), # L
	"V"         => 5000.0                                         # max volume, L
)

# ╔═╡ a1ca6688-ed64-4758-9828-c9dc8734f9a4
md"
🍷 next, we write `JuMP.jl` code to formulate then numerically solve the linear program.
"

# ╔═╡ 643a18d3-69fc-4560-8d17-51f861b6e9fc
begin
	# initialize the optimization model
	model = Model(HiGHS.Optimizer)
	set_silent(model)

	# attach decision variables
	wines = ["Syrah", "Grenache"]
	@variable(model, V[wines] .≥ 0.0) # volumes of each variety

	V_total = sum(V[wine] for wine in wines) # total volume of blend

	# attach objective function
	price = 8.0 # $
	@objective(
		model,
		Max,
		price * V_total - sum(data["cost"][wine] * V[wine] for wine in wines)
	) # profit

	#=
	attach constraints
	=#
	@constraint(
		model, 
		sum(V[wine] for wine in wines) ≤ data["V"]
	) # total volume

	min_alcohol = 11.5 # % by volume
	@constraint(
		model, 
		sum(V[wine] * data["alcohol"][wine] for wine in wines) ≥ min_alcohol * V_total
	) # alcohol

	max_acid = 0.36 # g/L
	@constraint(
		model, 
		sum(V[wine] * data["acid"][wine] for wine in wines) ≤ max_acid * V_total
	) # acid
	
	for wine in wines
		@constraint(model, V[wine] ≤ data["inventory"][wine]) # inventory
	end

	# numerically solve the linear program
	optimize!(model)

	# print the linear program
	latex_formulation(model)
end

# ╔═╡ 604fe9ed-cfd0-4475-b3a5-e8095ea2f9af
is_solved_and_feasible(model) # better be true!

# ╔═╡ 5c8d6ff0-f2e9-4848-8a85-a77c47c0e9fb
md"🍷 we print the optimal volume of Syrah and Grenache to include in the blend, the profit we reap by selling it to the customer, and the acid & alcohol content of the blend."

# ╔═╡ 46e53dad-15f4-46b4-850c-9db6a75390c7
for wine in wines
	println(wine, ": ", round(value(V[wine]), digits=1), " L")
end

# ╔═╡ 79b98c60-7dcd-4871-95ed-6ac0f727a028
println("total volume: ", round(value(V_total), digits=1))

# ╔═╡ 4c5f4cc8-d768-4dc2-b09b-1200246019e5
begin
	profit = round(objective_value(model), digits=2) # $
	println("profit: \$", profit)
end

# ╔═╡ 3687fd3d-497b-4c90-9761-58f24e492d29
begin
	println("specs of the optimal blend:")
	
	total_acid = sum(value(V[wine]) * data["acid"][wine] for wine in wines) # g
	println("\tacid [g/L]: ", total_acid / value(V_total))

	total_alcohol = sum(value(V[wine]) * data["alcohol"][wine] for wine in wines) # %⋅L
	println("\talcohol [% by volume]: ", total_alcohol / value(V_total))
end

# ╔═╡ f8d24f48-f907-4d7c-9166-9b97831cec30
md"🍷 we visualize the optimal wine blend with a bar plot."

# ╔═╡ 94bb6e34-6e2c-48b1-8d71-7748a732c824
barplot(
	1:length(wines), [value(V[wine]) for wine in wines], color="darkred",
	axis=(; 
		  xticks=(1:length(wines), wines), 
		  ylabel="volume [L]", 
		  title="optimal blend.\nprofit=\$$profit", 
		  limits=(nothing, nothing, 0, nothing),
		  titlefont=:regular
	)
)

# ╔═╡ 4141b7f2-734d-4cf3-a21b-215ad0cb62ab
md"🍷 which constraints were active? which were inactive?

the _active_ constraints are the maximal amount of blend the customer will accept and our inventory of Syrah, the cheaper wine variety.

the acid and alcohol specifications constituted _inactive_ constraints.
the Grenache inventory is also an inactive constraint.
"

# ╔═╡ e12baa22-e32b-4fda-bd9a-3df2afb402ef
md"🍷 since our design space is two dimensional, we can visualize this linear program. we draw the 2D design space, the constraint lines, the feasible region, some profit contours, and the optimal solution.

the _feasible region_ of wine blend designs constitutes the shaded gray region, where
* the blend has sufficient alcohol content
* the blend exhibits sufficiently low acid content
* the customer's capacity to buy the wine blend is not exceeded
* our inventory of each wine is respected
the feasible region is a polygon.

some _profit_ contours are shown as the dashed lines. blends belonging to a profit contour all give us the same profit.

the _optimal solution_ lies at a vertex of the feasible region. at this point, we cannot slide the profit contour to the upper right any further because all of the blends would be infeasible.
here, we clearly see that the only two active constraints are the customer demand and Syrah inventory. 

thinking about changing a constraint, we could make more profit if we had more Syrah in our inventory! with more Syrah, eventually, the acid constraint would come into play and be active.
"

# ╔═╡ 9a0bd383-3358-4100-850b-5b2511310bc8
begin
	V_max = 6000 # span of design space to show
	
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel=wines[1] * " [L]", 
		ylabel=wines[2] * " [L]",
		limits=(0, V_max, 0, V_max), 
		aspect=DataAspect()
	)

	#=
	constraints
	=#
	V₁ = range(0.0, V_max, length=500)
	
	# alcohol constraint line
	alc_slope = -(data["alcohol"][wines[1]] - min_alcohol) / (data["alcohol"][wines[2]] - min_alcohol)
	lines!(V₁, alc_slope * V₁, label="alcohol constraint")
	
	# acid constraint line
	aci_slope = -(data["acid"][wines[1]] - max_acid) / (data["acid"][wines[2]] - max_acid)
	lines!(V₁, aci_slope * V₁, label="acid constraint")

	# total volume constraint line
	total_volume = data["V"] .- V₁
	lines!(V₁, total_volume, label="customer demand")

	# inventory constraint line
	vlines!(data["inventory"][wines[1]], color=Cycled(4), label="inventory")
	hlines!(data["inventory"][wines[2]], color=Cycled(4))
	
	# feasible solution
	ids = V₁ .< data["inventory"][wines[1]]
	band!(
		V₁[ids], 
		aci_slope * V₁[ids], 
		min.(total_volume[ids], data["inventory"][wines[2]]),
		color=("black", 0.1), label="feasible region"
	)

	# profit contours
	δs = [price - data["cost"][wine] for wine in wines]
	profit_slope = -δs[1] / δs[2]
	for p in [1.0, 0.75, 0.5]
		lines!(
			V₁, p * profit / δs[2] .+ profit_slope * V₁, 
			linestyle=:dash, 
			color=("gray", p), 
			label=p == 1.0 ? "profit contour" : nothing
		)
	end
	
	# optimal solution
	scatter!(
		[value(V[wines[1]])], [value(V[wines[2]])], 
		marker=:star5, color="black", markersize=25,
		label="optimal solution"
	)

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ a1495447-bf20-46f1-b971-bb45c3e310a6
md"🍷 we gain a lot of insight by visualizing a linear program in 2D, as we did above. we can make sense of the solution and what constraints are active vs. inactive. 

for blends of more than three varieties, we lose this interpretability since we cannot visualize a design space that is greater than 4D.

interestingly, last month, I was drinking this Oregon wine that is composed of a blend of _nine_ different grapes❗ hard to imagine a nine-dimensional polytope comprising the feasible region.
"

# ╔═╡ 5dc5fa6c-1f52-447e-96fe-495c4dc336ad
html"<img src=\"https://raw.githubusercontent.com/SimonEnsemble/CHE361_W2025/refs/heads/main/drawings/wine_blend.jpeg\" width=400>"

# ╔═╡ dcf361c6-6cb2-4204-ad4e-934aaf824cc2
md"
❓ I wonder if any wineries actually take such a sophisticated approach to wine blending.

!!! note \"post-fermentation wine blends vs. co-fermentation\"
	we considered a _post_-fermentation wine blend. by contrast, some wine is made by blending the distinct varieties of grapes or grape juices together in the fermenter, then fermenting them together. in that case, I imagine it is much more difficult to predict what the specifications of the wine will be at the end of the fermentation.
"

# ╔═╡ Cell order:
# ╟─83d420f4-3c94-467f-8c5a-879b6201e74c
# ╟─837796eb-0f64-4e3d-8e45-fdf5ff3a358b
# ╟─600116f7-21c0-463f-a7ce-8edcace8542f
# ╠═7a733b92-35a6-11f0-0323-f9de5de91771
# ╟─f030c170-4f8c-4e70-86b1-0944935d550b
# ╠═711f0fb2-3a27-470e-ba4a-66b9f0ead816
# ╟─a1ca6688-ed64-4758-9828-c9dc8734f9a4
# ╠═643a18d3-69fc-4560-8d17-51f861b6e9fc
# ╠═604fe9ed-cfd0-4475-b3a5-e8095ea2f9af
# ╟─5c8d6ff0-f2e9-4848-8a85-a77c47c0e9fb
# ╠═46e53dad-15f4-46b4-850c-9db6a75390c7
# ╠═79b98c60-7dcd-4871-95ed-6ac0f727a028
# ╠═4c5f4cc8-d768-4dc2-b09b-1200246019e5
# ╠═3687fd3d-497b-4c90-9761-58f24e492d29
# ╟─f8d24f48-f907-4d7c-9166-9b97831cec30
# ╠═94bb6e34-6e2c-48b1-8d71-7748a732c824
# ╟─4141b7f2-734d-4cf3-a21b-215ad0cb62ab
# ╟─e12baa22-e32b-4fda-bd9a-3df2afb402ef
# ╠═9a0bd383-3358-4100-850b-5b2511310bc8
# ╟─a1495447-bf20-46f1-b971-bb45c3e310a6
# ╟─5dc5fa6c-1f52-447e-96fe-495c4dc336ad
# ╟─dcf361c6-6cb2-4204-ad4e-934aaf824cc2
