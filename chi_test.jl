### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 85db4222-a7d1-11f1-8861-9b8b1c1400b0
begin
	import Pkg; Pkg.activate()
	using UnicodePlots, StatsBase, CairoMakie
	update_theme!(theme_minimal())
	update_theme!(fontsize=20, linewidth=3)
end

# ╔═╡ e8b6036c-91f1-4c75-91d1-c1a5f6cbf9b6
md"
## customers pick drinks at a bar
just opened a bar! we serve three drinks. each individual customer lines up and orders a single drink.
"

# ╔═╡ b42b91fc-e910-47b6-bd83-7fcbd679633d
drinks = ["🍷", "🍺", "🍸"]

# ╔═╡ 80df9ccd-1e9b-4088-80a0-dcf36672cff9
md"
## null hypothesis
we have a null hypothesis: all drinks are equally preferable. under this hypothesis, a customer selects a drink according to the uniform distribution over drinks.
"

# ╔═╡ af287958-b529-48ce-bc68-3ba4e3fd36fe
p_null = Dict("🍷" => 1/3, "🍺" => 1/3, "🍸" => 1/3) # drink-selection probabilities

# ╔═╡ 29296657-243f-488a-b942-55f5b7587903
md"
under the null hypothesis, we can simulate outcomes of customers purchasing drinks. 

we assume each customer chooses a drink independently of the others. (not true if, e.g., customers behind me see me choose 🍷 then think \"the cool ppl are choosing🍷, I'll order that too.\")
"

# ╔═╡ 0e2bc51b-35bc-44d8-8934-52a8fa952d7f
function simulate_null(n_customers)
	data = Dict(drink => 0 for drink in drinks)

	ps = ProbabilityWeights([p_null[drink] for drink in drinks])
	for c = 1:n_customers
		drink = sample(drinks, ps)
		data[drink] += 1
	end
	
	return data
end

# ╔═╡ 7a2f8f51-53e6-44a9-abd3-5fa756b49f41
md"simulation \#1"

# ╔═╡ d805d8af-baf8-43f3-8547-13f4e8de51fa
sim_data = simulate_null(30)

# ╔═╡ 516584c9-2071-4a15-a218-6792ca6c8f27
UnicodePlots.barplot(drinks, [sim_data[drink] for drink in drinks])

# ╔═╡ 5300bd90-7801-4f2b-afc3-22f271be0151
md"simulation \#2"

# ╔═╡ a042d091-7aa4-47d0-b1c4-8548347a84e7
sim_data2 = simulate_null(30)

# ╔═╡ ccc5e09d-da8a-47df-a3cd-521676597ab2
UnicodePlots.barplot(drinks, [sim_data2[drink] for drink in drinks])

# ╔═╡ b7e8e6da-a6dd-4a86-b0f0-82628f4063b1
md"the outcome of the simulation is stochastic, and usually the drink counts are not the same." 

# ╔═╡ 638a7b16-d482-4c80-8e35-efce14408ba0
md"
## data
on our first day, 30 customers visited our bar. each purchased a drink. here is the empirical distribution."

# ╔═╡ 6f6f0451-751c-433a-8c5a-5d2048a671cb
data = Dict("🍷" => 14, "🍺" => 7, "🍸" => 9)

# ╔═╡ 46266490-6da4-48c6-a722-64b94fb6ea81
n_customers = sum(values(data))

# ╔═╡ 8cafed20-2c06-4812-b2e7-e7441ee69ab6
UnicodePlots.barplot(drinks, [data[drink] for drink in drinks])

# ╔═╡ 196a379d-cae2-4aca-a1fb-84a2690aeb84
md"## the statistical question

does the observed distribution among the drinks differ in a statistically signficant way from the expected distribution under the null hypothesis? sure, our data does not have equal counts of drinks, but neither do the simualted data sets under the null hypothesis!
"

# ╔═╡ e14e0adc-872e-4304-bb14-452fea9759c4
md"## test statistic

for a total of 30 customers, define the test statistic:
 
 $\chi:=$ |🍷(observed) - 10|+ |🍺(observed) - 10|+ |🍸(observed) - 10|

to measure the difference between an observed distribution and the expected distribution under the null hypothesis.
"

# ╔═╡ 7c5d8b53-a0aa-41d0-9f97-69cd45af8305
χ_obs = sum(abs(data[drink] - n_customers / length(drinks)) for drink in drinks)

# ╔═╡ 6cb646c7-ff52-420d-9da5-32624caef9b3
md"## distribution of test statistic under null hypothesis

what does the distribution of the test statistic look like under the null hypothesis?
how does our observed test statistic compare?
"

# ╔═╡ 3a6ec4ff-b6d2-41ed-a390-588b3341f9e0
n_sims = 10_000

# ╔═╡ 100a41a0-4b68-4c47-8be8-661a637867f3
χs_null = [
	begin
		local sim_data = simulate_null(n_customers)
		χ = sum(
			abs(sim_data[drink] - n_customers / length(drinks)) 
			for drink in drinks
		)
	end
	for s = 1:n_sims
]

# ╔═╡ 1b93b0ce-4ed5-463d-9743-b443a2426bc0
begin
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="χ", ylabel="# outcomes", 
		limits=(0, nothing, 0, nothing)
	)
	hist!(χs_null)
	vlines!([χ_obs], color=Cycled(2), label=rich("χ", subscript("obs")))
	axislegend()
	fig
end

# ╔═╡ c97eab42-ce38-425c-824e-0c476d273d11
md"the observed test statistic falls pretty much in the middle of the distribution of the test statistic under the null hypothesis. that means our observed test statistic is consistent with the null hypothesis, actually. so we lack strong evidence to reject the null hypothesis based on our drink picks data with these 30 customers.

## p-value

the p-value is the probability that we'd get a test statistic equal to or more extreme than our observed test statistic if the null hypothesis were true. since the p-value is quite big, we cannot reject the null hypothesis.

if we collected way more data (i.e. more customers), we might have the power to reject the null hypothesis. not with this small data set.
"

# ╔═╡ be972bad-8c88-4aa5-a5a2-50a829aa18c0
p_val = sum(χs_null .≥ χ_obs) / n_sims

# ╔═╡ 773db6a7-6f75-4a29-a49b-b7426e5b8d91
md"## caveats

* we assumed independent decisions among the customers. probably the drinks that earlier customers chose influence the choices of later customers. 
* a single customer might order several drinks over the night. maybe the customer prefers variety or sticks to the same drink.
* drink preferences might change based on time of day, weather, season, day of the week.
* introducing a new drink might shift preferences. check out Luce's choice axiom.
"

# ╔═╡ Cell order:
# ╠═85db4222-a7d1-11f1-8861-9b8b1c1400b0
# ╟─e8b6036c-91f1-4c75-91d1-c1a5f6cbf9b6
# ╠═b42b91fc-e910-47b6-bd83-7fcbd679633d
# ╟─80df9ccd-1e9b-4088-80a0-dcf36672cff9
# ╠═af287958-b529-48ce-bc68-3ba4e3fd36fe
# ╟─29296657-243f-488a-b942-55f5b7587903
# ╠═0e2bc51b-35bc-44d8-8934-52a8fa952d7f
# ╟─7a2f8f51-53e6-44a9-abd3-5fa756b49f41
# ╠═d805d8af-baf8-43f3-8547-13f4e8de51fa
# ╠═516584c9-2071-4a15-a218-6792ca6c8f27
# ╟─5300bd90-7801-4f2b-afc3-22f271be0151
# ╠═a042d091-7aa4-47d0-b1c4-8548347a84e7
# ╠═ccc5e09d-da8a-47df-a3cd-521676597ab2
# ╟─b7e8e6da-a6dd-4a86-b0f0-82628f4063b1
# ╟─638a7b16-d482-4c80-8e35-efce14408ba0
# ╠═6f6f0451-751c-433a-8c5a-5d2048a671cb
# ╠═46266490-6da4-48c6-a722-64b94fb6ea81
# ╠═8cafed20-2c06-4812-b2e7-e7441ee69ab6
# ╟─196a379d-cae2-4aca-a1fb-84a2690aeb84
# ╟─e14e0adc-872e-4304-bb14-452fea9759c4
# ╠═7c5d8b53-a0aa-41d0-9f97-69cd45af8305
# ╟─6cb646c7-ff52-420d-9da5-32624caef9b3
# ╠═3a6ec4ff-b6d2-41ed-a390-588b3341f9e0
# ╠═100a41a0-4b68-4c47-8be8-661a637867f3
# ╠═1b93b0ce-4ed5-463d-9743-b443a2426bc0
# ╟─c97eab42-ce38-425c-824e-0c476d273d11
# ╠═be972bad-8c88-4aa5-a5a2-50a829aa18c0
# ╟─773db6a7-6f75-4a29-a49b-b7426e5b8d91
