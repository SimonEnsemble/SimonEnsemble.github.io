### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ d0ea1620-1bc5-11ef-0988-0394597d4f32
begin
	import Pkg; Pkg.activate()
	using CairoMakie, JuMP, ColorSchemes, PlutoTeachingTools, PlutoUI
	import HiGHS;
end

# ╔═╡ 5b53b7a7-af29-45f0-ada1-a4f94a664fe0
color = Dict(zip(["hops", "malt", "corn"], ColorSchemes.Accent_3))

# ╔═╡ 37cc8195-e3ba-47b5-9c76-e2dc8d4b2ac0
set_theme!(theme_minimal()); update_theme!(fontsize=20, linewidth=3)

# ╔═╡ fbd25cce-79de-4dca-8895-27bd2d3a63db
TableOfContents()

# ╔═╡ 56385e01-63db-47b3-b630-c7a05cfb888c
md"
# numerically solving a linear program in Julia

[learn by example!]

!!! note
	to learn about the theory behind linear programming, I recommend the textbook (Ch. 7):
	> S. Dasgupta, C. H. Papadimitriou, and U. V. Vazirani. Algorithms. (2006)

## 🍺 problem setup

we are a small brewery capable of producing both ale and lager.

we wish to determine the number of barrels of both ale and lager to produce in order to maximize our profit. producing then selling a barrel of ale gives us \$39 profit, while a barrel of lager gives us \$69 profit.

however, ale and lager require different resources to produce, and we are constrained in our resources. a barrel of ale requires 5 kg corn, 4 kg hops, and 35 kg of barley malt. a barrel of lager requires 15 kg corn, 4 kg hops, and 20 kg of malt. meanwhile, our brewery is in possession of only 480 kg corn, 160 kg hops, and 1190 kg of malt.
"

# ╔═╡ 633ae1e9-d65b-4561-bed5-8dc5eb87039c
Foldable(
	"references",
	md"this linear programming problem originates from \"Robert Bland, Allocation of Resources by Linear Programming, _Scientific American_, Vol. 244, No. 6, June 1981\".
	it was covered in class notes [here](https://dmi.unibas.ch/fileadmin/user_upload/dmi/Studium/Computer_Science/Vorlesung_HS20/Main_Lecture_Scientific_Computing/CompSci20-Chapter03.pdf) and [here](https://www.cs.princeton.edu/courses/archive/spring03/cs226/lectures/lp-4up.pdf).
	"
)

# ╔═╡ 69ed26e0-4c02-4d6a-85e3-54af0e3bc9f3
md"
## storing the data in dictionaries

🍺 let's define some dictionaries that contain the data in the problem setup for later use. (also, a list of the resources and brews to loop over)."

# ╔═╡ d95c0541-047f-4ece-bb12-43cc825d6a3d
# list of resources
resources = ["corn", "hops", "malt"]

# ╔═╡ 0f41fca8-fffc-41a2-8bf6-f07894a9c00c
# list of brews
brews = ["ale", "lager"]

# ╔═╡ b35a803a-e611-43e2-a46e-153c0a3af8fd
# amount of resources we have
our_resources = Dict(
	"corn" => 480.0,  # kg
	"hops" => 160.0,  # kg
	"malt" => 1190.0  # kg
)

# ╔═╡ 4768889e-12f1-4506-a8fb-9206d110fd68
# resource requirements for producing ale and lager
resource_requirements = Dict(
	"ale" => Dict(
		"corn" => 5.0,  # kg/barrel
		"hops" => 4.0,  # kg/barrel
		"malt" => 35.0  # kg/barrel
	),
	"lager" => Dict(
		"corn" => 15.0,  # kg/barrel
		"hops" => 4.0,   # kg/barrel
		"malt" => 20.0   # kg/barrel
	)
)

# ╔═╡ 6b59850d-01ea-44f2-86bc-578866d5efc2
# profits from producing ale and lager
profit = Dict(
	"ale"   => 39.0, # $
	"lager" => 69.0, # $
)

# ╔═╡ 48d41db7-bf31-4dd5-860a-0862b4282227
md"## two extremes: all lagers or all ales"

# ╔═╡ f14b1ccf-4ca8-4444-b0c5-f1181d26173e
md"🍺 if our brewery were to devote all of its resources to ale production, how much profit would we make? how much corn, hops, and malt would we have leftover?

...what if we devoted all of our resources to lager production?
"

# ╔═╡ 0ea056a7-6e71-40b0-9584-1ec5104be3aa
for brew in brews
	# number of brews we could produce
	#   compute how many beers we can brew for each resource.
	#   take the min, corresponding to the limiting resource
	nb_brews = minimum(
		[our_resources[res] / resource_requirements[brew][res] for res in resources]
	)

	# compute the profit we'd make
	our_profit = nb_brews * profit[brew]

	# compute resources we'd have left
	resources_left = Dict(
		res => our_resources[res] - nb_brews * resource_requirements[brew][res]
		for res in resources
	)

	println(
		"💡 devoting all of our resources to $brew, we'd produce $nb_brews barrels of $brew and make $our_profit profit. we'd have $(resources_left["corn"]) kg corn, $(resources_left["hops"]) kg hops, $(resources_left["malt"]) kg malt leftover.\n"
	)
end

# ╔═╡ 358daf8d-bb21-4d19-a98a-404ea7c764e5
md"## feasibility of solutions"

# ╔═╡ 13f00426-1ed3-4919-bf5a-0f32c367245e
md"🍺 one employee uses their 'gut feeling' to propose we produce 10 ales and 40 lagers to make a profit of \$3150. explain the reason(s) why this proposed solution is _infeasible_.
"

# ╔═╡ 50abaf4b-1cbf-48ae-9650-6d4cf700c157
# the infeasible solution
infeasible_soln = Dict("ale" => 10, "lager" => 40)

# ╔═╡ 7a9e35f5-ea86-42a8-b4a1-c074dcf0055a
# the profit we'd make, if this were feasible
sum(profit[b] * infeasible_soln[b] for b in brews)

# ╔═╡ 31f21065-6d83-4180-8382-b26af53a50c3
begin
	println(
		"to produce $(infeasible_soln["ale"]) barrels of ale and $(infeasible_soln["lager"]) barrels of lager..."
	)
	# loop over resources
	for res in resources
		# calculate the amount of this resource req'd for this solution
		res_required = sum(
			resource_requirements[b][res] * infeasible_soln[b] 
				for b in brews
		)
		
		# check feasibility
		if res_required > our_resources[res]
			println(
				"\t$res_required kg of $res are required, but we only have $(our_resources[res]) kg at our brewery!"
			)
		end
	end
end

# ╔═╡ 2b345f5b-9b3b-46df-a220-bfa2e64bf2fb
md"## numerically solve the linear program

let's use `JuMP.jl`, an optimization library in Julia, to numerically solve this linear program.

!!! note \"JuMP.jl docs on linear programs\"
	see the `JuMP.jl` docs on linear programming with several examples [here](https://jump.dev/JuMP.jl/stable/tutorials/linear/introduction/).
"

# ╔═╡ 9537f6c9-5c28-48d5-81c5-3d908172e4d1
md"🍺 create a `JuMP.jl` optimization model with the `HiGHS` optimizer (open-source software that can solve linear programs)."

# ╔═╡ d37cb850-0eb5-4c36-8832-4d05c6960eb7
model = Model(HiGHS.Optimizer)

# ╔═╡ 90a1c7db-aaf0-43b8-a185-bcef96a2e69b
md"🍺 tell `JuMP.jl` about your decision variables."

# ╔═╡ 2eaaee47-aaab-4a4f-8967-232e2a3d52ca
@variable(model, x[brews] .≥ 0.0) # could specify as an Int here

# ╔═╡ 37430457-1186-4518-a4dd-950277fb4ca6
x["ale"], x["lager"] # there ya go

# ╔═╡ 2b74698b-5202-4ee8-a74a-9d66b778b014
md"🍺 tell `JuMP.jl` about your objective function."

# ╔═╡ a9aa7b20-6f63-4cb8-a0a5-8c9672e7bbce
@objective(model, Max, sum(x[brew] * profit[brew] for brew in brews))

# ╔═╡ ebf3cc2a-2a43-4985-822e-4a164db20513
md"🍺 tell `JuMP.jl` about your constraints."

# ╔═╡ 12a9f5b9-ab35-4453-b00f-fc7d7585294e
for res in resources
	@constraint(model,
		sum(resource_requirements[brew][res] * x[brew] for brew in brews) ≤ our_resources[res]
	)
end

# ╔═╡ dc6725f0-a7f4-48b1-bda8-bc61d866c5cf
md"🍺 ask `JuMP.jl` to display your linear program in mathematical form using `latex_formulation`."

# ╔═╡ 49819512-018c-4b06-8e2d-118c05f0b9ee
latex_formulation(model)

# ╔═╡ 17ab7aff-a6a9-4914-82ae-301266a5fe53
md"🍺 ask `JuMP.jl` to search for and find the optimal solution."

# ╔═╡ b372aa05-5a36-4c17-bf2e-9ba8bc2f4f28
optimize!(model)

# ╔═╡ 323c4baa-6afb-4002-bd05-edbd97eb947b
md"🍺 check with `JuMP.jl` that your linear program is actually feasible and that it has been successfully solved."

# ╔═╡ 32a4a7eb-d46e-480f-8c83-5969c2888cd0
is_solved_and_feasible(model)

# ╔═╡ 33ce18b0-ab90-43d5-9592-788f849f36eb
md"🍺 for the optimal solution,
* how many barrels of ale and lager should we produce to maximize profit?
* what is the profit we will make? [check that this profit improves upon the all-ale and all-lager strategy]
* how much corn, hops, and barley will be leftover?
"

# ╔═╡ e5874d3e-f89e-4ca6-b97b-9391adef8615
opt_soln = Dict(brew => value(x[brew]) for brew in brews)

# ╔═╡ 984ca36e-ae08-492f-9e62-257a00cdeafc
opt_profit = solution_summary(model).objective_value

# ╔═╡ 797a4e4e-84bf-450c-9653-149395cda551
barplot(
	1:2, [opt_soln[brew] for brew in brews],
	axis=(; 
		ylabel="# barrels", xticks=(1:2, brews), title="optimal solution", limits=(nothing, nothing, 0, nothing)
	)
)

# ╔═╡ e1dc7caa-9ad5-4ca3-a90e-80d90e54fe82
resources_left = Dict(
	res => our_resources[res] - sum(
		opt_soln[brew] * resource_requirements[brew][res]
		for brew in brews
	)
	for res in resources
)

# ╔═╡ fb129a16-3390-42cd-898c-56e75d18ad20
barplot(
	1:3, [resources_left[res] for res in resources], color=Cycled(2),
	axis=(; 
		ylabel="kg", xticks=(1:3, resources), title="resources left", limits=(nothing, nothing, 0, nothing)
	)
)

# ╔═╡ 8cd2c9c9-473a-4d53-9696-cfc9de8edeb6
md"🎆 the optimal solution is to brew $(round(Int, opt_soln[\"ale\"])) barrels of ale and $(round(Int, opt_soln[\"lager\"])) barrels of lager. this gives us a profit of \$$opt_profit. we'd have $(resources_left[\"corn\"]) kg corn, $(round(abs(resources_left[\"hops\"]))) kg hops, and $(round(abs(resources_left[\"malt\"]))) kg malt leftover."

# ╔═╡ 4ac7625d-9cda-48ed-bf46-c94d453ac01f
md"## visualization of the feasible region and optimal solution
since this linear program involves only two decision variables, we can visualize the constraints, the feasible region, and the solution in decision variable space.

🍺 visualize!
"

# ╔═╡ 929802a9-d2a0-4df1-8368-02118961378e
begin
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="ales [# barrels]", 
		ylabel="lagers [# barrels]",
		aspect=DataAspect()
	)
	
	# for plotting constraints and profit lines
	ale_barrels = collect(range(0.0, 35.0, length=400))
	
	# plot constraints
	l_lo = zeros(length(ale_barrels))
	l_hi = [Inf for _ = 1:length(ale_barrels)]
	for res in resources
		# constraint line: l = m * a + b
		m = - resource_requirements["ale"][res] / resource_requirements["lager"][res]
		b = our_resources[res] / resource_requirements["lager"][res]
		l_hi_local = m * ale_barrels .+ b
		
		lines!(ale_barrels, l_hi_local, label=res * " constraint", color=color[res])
		
		l_hi = min.(l_hi, l_hi_local)
	end
	band!(ale_barrels, l_lo, l_hi, color=("gray", 0.1), label="feasible region")
	
	# plot profit contour l = m * a + b
	b = opt_profit / profit["lager"]
	m = -profit["ale"] / profit["lager"]
	lines!(
		ale_barrels, m * ale_barrels .+ b, 
		linestyle=:dash, color="black", label="opt profit contour"
	)
	
	xlims!(0, nothing)
	ylims!(0, nothing)
	
	# opt soln
	scatter!(
		[opt_soln["ale"]], [opt_soln["lager"]], 
		color="red", marker=:star6, markersize=28, label="opt sol'n"
	)
	
	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 29d219d9-5ddd-496b-a23d-e5ab9f61c3aa
md"## sensitivity analysis

the profit we'll make from producing a barrel of ale and lager is subject to uncertainty because the prices fluctuate in the market. hence, our objective function is subject to uncertainty.

🍺 by how much can the profit we'd make on a barrel of lager _increase_ and _decrease_ (keeping the profit on ale fixed) before our optimal solution would change? how about for ale?

!!! note
	the `lp_sensitivity_report` in `JuMP.jl` gives us this information. see [the docs](https://jump.dev/JuMP.jl/stable/tutorials/linear/lp_sensitivity/).
"

# ╔═╡ 790f2d5d-d5d5-4c26-81a1-4dee94c7cf11
sensitivity = lp_sensitivity_report(model)

# ╔═╡ 6527bd28-6fd5-4a68-8877-6b9ef3995e31
for brew in brews
	println(
		"💡 holding the profit we make on a barrel of $(filter(b -> b != brew, brews)[1]) fixed, the profit we make on a barrel of $brew could decrease by <\$$(-round(Int, sensitivity[x[brew]][1])) or increase by <\$$(round(Int, sensitivity[x[brew]][2])) without changing the optimal solution.\n"
	)
end

# ╔═╡ Cell order:
# ╠═d0ea1620-1bc5-11ef-0988-0394597d4f32
# ╠═5b53b7a7-af29-45f0-ada1-a4f94a664fe0
# ╠═37cc8195-e3ba-47b5-9c76-e2dc8d4b2ac0
# ╠═fbd25cce-79de-4dca-8895-27bd2d3a63db
# ╟─56385e01-63db-47b3-b630-c7a05cfb888c
# ╟─633ae1e9-d65b-4561-bed5-8dc5eb87039c
# ╟─69ed26e0-4c02-4d6a-85e3-54af0e3bc9f3
# ╠═d95c0541-047f-4ece-bb12-43cc825d6a3d
# ╠═0f41fca8-fffc-41a2-8bf6-f07894a9c00c
# ╠═b35a803a-e611-43e2-a46e-153c0a3af8fd
# ╠═4768889e-12f1-4506-a8fb-9206d110fd68
# ╠═6b59850d-01ea-44f2-86bc-578866d5efc2
# ╟─48d41db7-bf31-4dd5-860a-0862b4282227
# ╟─f14b1ccf-4ca8-4444-b0c5-f1181d26173e
# ╠═0ea056a7-6e71-40b0-9584-1ec5104be3aa
# ╟─358daf8d-bb21-4d19-a98a-404ea7c764e5
# ╟─13f00426-1ed3-4919-bf5a-0f32c367245e
# ╠═50abaf4b-1cbf-48ae-9650-6d4cf700c157
# ╠═7a9e35f5-ea86-42a8-b4a1-c074dcf0055a
# ╠═31f21065-6d83-4180-8382-b26af53a50c3
# ╟─2b345f5b-9b3b-46df-a220-bfa2e64bf2fb
# ╟─9537f6c9-5c28-48d5-81c5-3d908172e4d1
# ╠═d37cb850-0eb5-4c36-8832-4d05c6960eb7
# ╟─90a1c7db-aaf0-43b8-a185-bcef96a2e69b
# ╠═2eaaee47-aaab-4a4f-8967-232e2a3d52ca
# ╠═37430457-1186-4518-a4dd-950277fb4ca6
# ╟─2b74698b-5202-4ee8-a74a-9d66b778b014
# ╠═a9aa7b20-6f63-4cb8-a0a5-8c9672e7bbce
# ╟─ebf3cc2a-2a43-4985-822e-4a164db20513
# ╠═12a9f5b9-ab35-4453-b00f-fc7d7585294e
# ╟─dc6725f0-a7f4-48b1-bda8-bc61d866c5cf
# ╠═49819512-018c-4b06-8e2d-118c05f0b9ee
# ╟─17ab7aff-a6a9-4914-82ae-301266a5fe53
# ╠═b372aa05-5a36-4c17-bf2e-9ba8bc2f4f28
# ╟─323c4baa-6afb-4002-bd05-edbd97eb947b
# ╠═32a4a7eb-d46e-480f-8c83-5969c2888cd0
# ╟─33ce18b0-ab90-43d5-9592-788f849f36eb
# ╠═e5874d3e-f89e-4ca6-b97b-9391adef8615
# ╠═984ca36e-ae08-492f-9e62-257a00cdeafc
# ╠═797a4e4e-84bf-450c-9653-149395cda551
# ╠═e1dc7caa-9ad5-4ca3-a90e-80d90e54fe82
# ╠═fb129a16-3390-42cd-898c-56e75d18ad20
# ╟─8cd2c9c9-473a-4d53-9696-cfc9de8edeb6
# ╟─4ac7625d-9cda-48ed-bf46-c94d453ac01f
# ╠═929802a9-d2a0-4df1-8368-02118961378e
# ╟─29d219d9-5ddd-496b-a23d-e5ab9f61c3aa
# ╠═790f2d5d-d5d5-4c26-81a1-4dee94c7cf11
# ╠═6527bd28-6fd5-4a68-8877-6b9ef3995e31
