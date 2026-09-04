### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 85db4222-a7d1-11f1-8861-9b8b1c1400b0
begin
	import Pkg; Pkg.activate()
	using UnicodePlots, StatsBase, CairoMakie, Random
	Random.seed!(1337)
	set_theme!(theme_minimal())
	update_theme!(fontsize=20, linewidth=3)
end

# ╔═╡ e8b6036c-91f1-4c75-91d1-c1a5f6cbf9b6
md"

# does an observed distribution over categories match an expected distribution?

💡 topic: a simulation-based approach for statistically testing the hypothesis that an observed distribution over categories differs from an expected distribution. 

!!! note \"Reference\"
	I learned this from the textbook \"Understanding Data\" by Alan Garfinkel and Yina Guo (Ch. 8).


## customers pick drinks at a bar
say we just opened a bar, and we serve three drinks. each customer orders a single drink.
"

# ╔═╡ b42b91fc-e910-47b6-bd83-7fcbd679633d
drinks = ["🍷", "🍺", "🍸"]

# ╔═╡ 80df9ccd-1e9b-4088-80a0-dcf36672cff9
md"
## null hypothesis
our null hypothesis is: all drinks are equally preferable; i.e. a customer picks a drink at random.
"

# ╔═╡ af287958-b529-48ce-bc68-3ba4e3fd36fe
# null hypothesis: uniform distribution over the drinks
p_null = Dict("🍷" => 1/3, "🍺" => 1/3, "🍸" => 1/3) 

# ╔═╡ 29296657-243f-488a-b942-55f5b7587903
md"
under the null hypothesis, we can simulate customers picking drinks. we assume each customer chooses a drink independently of the others.
"

# ╔═╡ 0e2bc51b-35bc-44d8-8934-52a8fa952d7f
function simulate_null(n_customers::Int, p_null::Dict{String, Float64})
	data = Dict(drink => 0 for drink in drinks)

	ps = ProbabilityWeights([p_null[drink] for drink in drinks])
	for c = 1:n_customers
		drink = sample(drinks, ps)
		data[drink] += 1
	end
	
	return data
end

# ╔═╡ 7a2f8f51-53e6-44a9-abd3-5fa756b49f41
md"below, we conduct a few simulations of drink-picking by 30 customers under the null hypothesis and inspect the resulting distribution over drinks.

simulation \#1
"

# ╔═╡ d805d8af-baf8-43f3-8547-13f4e8de51fa
sim_data = simulate_null(30, p_null)

# ╔═╡ 516584c9-2071-4a15-a218-6792ca6c8f27
UnicodePlots.barplot(drinks, [sim_data[drink] for drink in drinks])

# ╔═╡ 5300bd90-7801-4f2b-afc3-22f271be0151
md"simulation \#2"

# ╔═╡ ccc5e09d-da8a-47df-a3cd-521676597ab2
begin
	sim_data2 = simulate_null(30, p_null)
	UnicodePlots.barplot(drinks, [sim_data2[drink] for drink in drinks])
end

# ╔═╡ 88ee974b-c023-4d04-ad44-861b307d9d26
md"simulation \#3"

# ╔═╡ 9ca6382b-033f-4fe1-9539-8be73b781dd9
begin
	sim_data3 = simulate_null(30, p_null)
	UnicodePlots.barplot(drinks, [sim_data3[drink] for drink in drinks])
end

# ╔═╡ b7e8e6da-a6dd-4a86-b0f0-82628f4063b1
md"the outcome of the simulation is stochastic, and usually the drink counts are not the same despite the assumed uniform distribution." 

# ╔═╡ 638a7b16-d482-4c80-8e35-efce14408ba0
md"
## data
on our first day, 30 customers visited our bar. each purchased a drink. here is the observed distribution over drinks."

# ╔═╡ 6f6f0451-751c-433a-8c5a-5d2048a671cb
data = Dict("🍷" => 14, "🍺" => 7, "🍸" => 9)

# ╔═╡ 46266490-6da4-48c6-a722-64b94fb6ea81
n_customers = sum(values(data))

# ╔═╡ 8cafed20-2c06-4812-b2e7-e7441ee69ab6
UnicodePlots.barplot(drinks, [data[drink] for drink in drinks])

# ╔═╡ 196a379d-cae2-4aca-a1fb-84a2690aeb84
md"## the statistical question

does the observed distribution among the drinks differ in a statistically significant way from the expected distribution under the null hypothesis? sure, our data does not have equal counts of drinks, but neither do the simulated data sets under the null hypothesis!
"

# ╔═╡ e14e0adc-872e-4304-bb14-452fea9759c4
md"## test statistic

for a total of 30 customers, define the test statistic:
 
 $\chi:=$ |🍷(observed) - 10|+ |🍺(observed) - 10|+ |🍸(observed) - 10|

to measure the difference between an observed distribution and the expected distribution under the null hypothesis.
"

# ╔═╡ 0f79145e-3209-4003-96b0-e6ce6568b45d
function calculate_χ(data::Dict{String, Int}, p_null::Dict{String, Float64})
	n_customers = sum(values(data))
	return sum(
		abs(data[drink] - p_null[drink] * n_customers) 
		for drink in keys(data)
	)
end

# ╔═╡ 7c5d8b53-a0aa-41d0-9f97-69cd45af8305
χ_obs = calculate_χ(data, p_null)

# ╔═╡ 6cb646c7-ff52-420d-9da5-32624caef9b3
md"## distribution of test statistic under null hypothesis

what does the distribution of the test statistic look like under the null hypothesis?
how does our observed test statistic compare?
"

# ╔═╡ 3a6ec4ff-b6d2-41ed-a390-588b3341f9e0
n_sims = 10_000

# ╔═╡ 100a41a0-4b68-4c47-8be8-661a637867f3
χs_null = [
	calculate_χ(simulate_null(n_customers, p_null), p_null)
	for s = 1:n_sims
]

# ╔═╡ 1b93b0ce-4ed5-463d-9743-b443a2426bc0
function viz_χ_distn(χs_null::Vector{Float64}, χ_obs::Float64)
	χ_counts = countmap(χs_null)
	χ_vals = sort(collect(keys(χ_counts)))
	χ_freq = [χ_counts[χ] for χ in χ_vals]

	p = sum(χs_null .≥ χ_obs) / length(χs_null)

	fig = Figure(size=(500, 400))
	ax = Axis(
		fig[1, 1],
		xlabel=L"\chi",
		ylabel="# simulated outcomes\nunder the null hypothesis",
		xticks=0:4:maximum(χ_vals),
		limits=(-1.5, nothing, 0, nothing),
		title="p-value: $(round(p, digits=3))"
	)
	CairoMakie.barplot!(ax, χ_vals, χ_freq, width=1, color="green")
	vlines!(ax, [χ_obs], color=Cycled(2), label=L"\chi_{obs}")
	axislegend(ax)
	fig
end

# ╔═╡ d379cf62-27c6-4908-9271-b32bfe80d0c7
viz_χ_distn(χs_null, χ_obs)

# ╔═╡ c97eab42-ce38-425c-824e-0c476d273d11
md"the observed test statistic falls inside the bulk of the distribution of the test statistic under the null hypothesis. that means our observed test statistic is consistent with the null hypothesis. so we lack strong evidence to reject the null hypothesis based on our drink-selection data pertaining to these 30 customers.

## p-value

the p-value is the probability that we'd get a test statistic equal to or more extreme than our observed test statistic if the null hypothesis were true.
"

# ╔═╡ be972bad-8c88-4aa5-a5a2-50a829aa18c0
p_val = sum(χs_null .≥ χ_obs) / n_sims

# ╔═╡ a6661bc2-b568-412f-90a0-a608293ab71b
md"since the p-value is quite big, we do not reject the null hypothesis under, say, a reasonable significance level of $\alpha=0.05$.

this doesn't mean the null hypothesis is true. if there really were underlying preference differences among the drinks, a more sizable data set (i.e. more customers) would give us the power to reject the null hypothesis. we show this next.
"

# ╔═╡ 79b2f0f3-309c-4bf9-93a1-b3fbdd94b807
md"## more data

suppose we instead ran a longer experiment over the time span of four days---giving a bigger data set. we construct the data so the distribution is the same shape as before.
"

# ╔═╡ 7ac73aeb-0201-4d32-a8ec-b59fc5279d01
data_big = Dict("🍷" => 14*4, "🍺" => 7*4, "🍸" => 9*4)

# ╔═╡ b5529b25-f9df-4e58-949e-9c4c03fec64f
UnicodePlots.barplot(drinks, [data_big[drink] for drink in drinks])

# ╔═╡ 0bf4b4fd-73c3-49d8-b73b-2d96284517c0
begin
	χ_obs_big = calculate_χ(data_big, p_null)
	χs_null_big = [
		calculate_χ(simulate_null(sum(values(data_big)), p_null), p_null)
		for s = 1:n_sims
	]
	viz_χ_distn(χs_null_big, χ_obs_big)
end

# ╔═╡ 3dc7fdf9-271c-441a-a536-709f31b864ef
md"with bigger data, we have sufficient evidence to reject the null hypothesis under a significance value of $\alpha=0.05$. we have strong evidence that drinks are not equally preferable."

# ╔═╡ 773db6a7-6f75-4a29-a49b-b7426e5b8d91
md"## caveats

* we assumed independent decisions among the customers. this is implicitly part of the null hypothesis. probably the drinks chosen by earlier customers influence the choices of later customers. for example, customers may see me drinking 🍷 then think \"the cool ppl are choosing 🍷, I'll order that too.\" people might order in groups too---and influence each others' decisions. 
* a single customer might order several drinks over the night. maybe the customer prefers variety or sticks to the same drink. the Chumbawamba song _Tubthumping_ suggests variety.
* drink preferences might change based on time of day, weather, season, and day of the week. this challenges the reliability of data collected over just a few days (owing to sampling bias).
* introducing a new drink might shift preferences. look up Luce's choice axiom, regularity, and the decoy effect.
"

# ╔═╡ 6a64fc28-a79d-4d6e-89b3-91b9fb516dbd
md"## Pearson's $\chi^2$ test

this resembles Pearson's $\chi^2$ test but we used absolute differences instead of squared differences in our test statistic.
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
# ╠═ccc5e09d-da8a-47df-a3cd-521676597ab2
# ╟─88ee974b-c023-4d04-ad44-861b307d9d26
# ╠═9ca6382b-033f-4fe1-9539-8be73b781dd9
# ╟─b7e8e6da-a6dd-4a86-b0f0-82628f4063b1
# ╟─638a7b16-d482-4c80-8e35-efce14408ba0
# ╠═6f6f0451-751c-433a-8c5a-5d2048a671cb
# ╠═46266490-6da4-48c6-a722-64b94fb6ea81
# ╠═8cafed20-2c06-4812-b2e7-e7441ee69ab6
# ╟─196a379d-cae2-4aca-a1fb-84a2690aeb84
# ╟─e14e0adc-872e-4304-bb14-452fea9759c4
# ╠═0f79145e-3209-4003-96b0-e6ce6568b45d
# ╠═7c5d8b53-a0aa-41d0-9f97-69cd45af8305
# ╟─6cb646c7-ff52-420d-9da5-32624caef9b3
# ╠═3a6ec4ff-b6d2-41ed-a390-588b3341f9e0
# ╠═100a41a0-4b68-4c47-8be8-661a637867f3
# ╠═1b93b0ce-4ed5-463d-9743-b443a2426bc0
# ╠═d379cf62-27c6-4908-9271-b32bfe80d0c7
# ╟─c97eab42-ce38-425c-824e-0c476d273d11
# ╠═be972bad-8c88-4aa5-a5a2-50a829aa18c0
# ╟─a6661bc2-b568-412f-90a0-a608293ab71b
# ╟─79b2f0f3-309c-4bf9-93a1-b3fbdd94b807
# ╠═7ac73aeb-0201-4d32-a8ec-b59fc5279d01
# ╠═b5529b25-f9df-4e58-949e-9c4c03fec64f
# ╠═0bf4b4fd-73c3-49d8-b73b-2d96284517c0
# ╟─3dc7fdf9-271c-441a-a536-709f31b864ef
# ╟─773db6a7-6f75-4a29-a49b-b7426e5b8d91
# ╟─6a64fc28-a79d-4d6e-89b3-91b9fb516dbd
