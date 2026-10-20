### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ ab2c1074-cc27-11f1-83ad-358a865b7c3d
begin
	import Pkg; Pkg.activate()
	using CairoMakie, Distributions, SpecialFunctions

	set_theme!(theme_light())
	update_theme!(fontsize=15)
end

# ╔═╡ 5d47e1be-f9ac-4942-b3bf-1169a9ed80c7
md"

# importance sampling

by Cory Simon & Qia Ke.

😱 suppose we are concerned with a probability density $\pi(x)$ over a sample space $\Omega$ with two issues:

1. we can't directly draw samples from it
2. we only know the shape of the density, up to a scaling factor. 

i.e.,

$$\pi(x) =  \gamma(x)/Z, \quad Z:=\int_{\Omega} \gamma(x)dx$$

with $\gamma(x)$ a function we can evaluate pointwise while the normalizing factor $Z$ is the value of an integral that (a) lacks a closed-form solution and (b) is high-dimensional, making it difficult to compute with numerical quadrature.

🎯 we want to estimate the expectation of a function:

$$\mathbb{E}_{x\sim \pi(x)}[f(x)]=\int_\Omega \pi(x) f(x) dx$$

as well as the normalizing factor $Z$.

!!! example
    as a toy example, suppose $\pi(x)$ is a Beta distribution with $\alpha=3$ and $\beta=6$ over $\Omega=[0, 1]$ and pretend we don't know the normalizing factor. we only know 
    
    $$\pi(x)\propto x^2(1-x)^5=\gamma(x)$$

    and let's say we want to know both the normalizing factor $Z$ (which is the Beta function evaluated at $\alpha=3$ and $\beta=6$) and the expectation of $f(x)=x$ under the beta distribution which is just its mean ($\alpha/(\alpha+\beta)$).

!!! note
    if we _could_ draw samples from $\pi(x)$ readily, we would just draw $N$ samples $\{x_i\}_{i=1}^N$ then approximate the expectation of $f(x)$ as:

    $$\mathbb{E}_{x\sim \pi(x)}[f(x)]\approx \frac{1}{N}\sum_{i=1}^N f(x_i).$$

    alas...
"

# ╔═╡ 0c563ea7-c248-42b8-801f-d1876038b9a5
begin
	α = 3.0
	β = 6.0
	function γ(x)
		if 0 ≤ x ≤ 1
			return x^(α - 1) * (1 - x)^(β - 1)
		else
			return 0.0
		end
	end
end

# ╔═╡ 9195342c-4bdb-46d0-aa75-4fe2b311f5dd
f(x) = x

# ╔═╡ 5db6f3ef-b573-4637-b948-41a2b5c49069
Z_true = beta(α, β) # hide

# ╔═╡ 2e0eab51-e0ff-4556-929f-742709df277b
𝔼f_true = α / (α + β) # hide

# ╔═╡ d6601702-dc6b-4512-a73d-902d475a8f62
md"🔨 we _do_ have an *importance density* $\eta(x)$ that we can draw samples from. crucially, $\pi(x) > 0 \implies \eta(x) > 0$ i.e. the support of $\eta(x)$ contains the support of $\pi(x)$. ideally, we engineer/craft $\eta(x)$ to be as close to $\pi(x)$ as we can.

let $\{x_i\}_{i=1}^N$ be a set of $N$ samples from the importance distribution, which we call _particles_.

!!! example
	let the importance distribution be a Gaussian distribution with mean 1/2 and standard deviation of 1/2. so 
	$$\eta(x) = \mathcal{N}\!\left(x;\, \tfrac{1}{2},\, \tfrac{1}{2}^2\right) = \frac{2}{\sqrt{2\pi}} \exp\left(-\frac{1}{2}\left(\frac{x - 1/2}{1/2}\right)^2\right)$$. and let's draw $N=1000$ samples from it.
"

# ╔═╡ f19d3174-008b-44de-8f77-89f1700ff5fc
η(x) = 2 / sqrt(2 * π) * exp(-((x - 0.5) / 0.5)^2 / 2)

# ╔═╡ 1b8e1902-5dd3-42ad-8599-74538c7d70b8
begin
	N = 1000
	xᵢs = rand(Normal(0.5, 0.5), N)

	hist(xᵢs, axis=(; xlabel="x", ylabel="# samples", title="samples from η(x)"))
end

# ╔═╡ 8f6a78a7-fed4-445e-943e-ac77c749e223
w(x) = γ(x) / η(x)

# ╔═╡ b850cb61-b9ea-45df-b864-3a7fc27a49ab
ws = [w(xᵢ) for xᵢ in xᵢs]

# ╔═╡ 3bf89820-c75e-4a54-aba5-02d0decf9944
w̃s = ws / sum(ws)

# ╔═╡ c6d01c1f-a7a0-4904-bbd4-88fef7feac1c
# approximation
Z = mean(ws)

# ╔═╡ c36f6626-69fa-4be4-8f7d-7325caf66b53
# approximation
𝔼f = sum(w̃s .* f.(xᵢs))

# ╔═╡ 3d24db8d-3841-47a1-9950-f9b52caec1af
md"

❓ how can we leverage these samples from $\eta(x)$ to compute the $Z$ and $\mathbb{E}_{x\sim \pi(x)}[f(x)]$ quantities we want?

💡 write them as expectations over $\eta(x)$ that we _do_ have samples from.

$$\mathbb{E}_{x\sim \pi(x)}[f(x)]=\int_\Omega \pi(x) f(x) \frac{\eta(x)}{\eta(x)} dx=\frac{1}{Z} \int_\Omega \eta(x) w(x) f(x) dx = \frac{1}{Z} \mathbb{E}_{x\sim \eta(x)}[w(x)f(x)] $$

with weight function

$$w(x)=\frac{\gamma(x) }{\eta(x)}$$

that we can readily evaluate.

likewise,

$$Z=\int_\Omega \gamma(x) \frac{\eta(x)}{\eta(x)}dx=\int_\Omega \eta(x) w(x) dx = \mathbb{E}_{x\sim \eta(x)}[w(x)].$$

by approximating $\eta(x)$ with the samples we have from it as a sum of Dirac delta functions at those samples

$$\eta(x)\approx \frac{1}{N}\sum_{i=1}^N\delta(x-x_i)$$

then stuffing this into the integrals above, we get the normalization factor as the average importance weight:

$$Z \approx \frac{1}{N}\sum_{i=1}^N w(x_i)$$

and the expectation of the function $f(x)$ as

$$\mathbb{E}_{x\sim \pi(x)}[f(x)] \approx \frac{1}{Z} \frac{1}{N}\sum_{i=1}^N w(x_i) f(x_i) \approx \sum_{i=1}^N \tilde{w}_i f(x_i)$$

where the *normalized importance weight* associated with particle $i$ is

$$\tilde{w}_i := w(x_i) / \sum_{j=1}^N w(x_j)$$

making the normalized importance weights sum to one.

!!! example
	let's define the weight function, compute the importance weights, and then use these as weights on the particles to compute the normalization factor and expectation of $x$ under our beta distribution!

	estimated Z = $(round(Z, sigdigits=3)) (true: $(round(Z_true, sigdigits=3))) 
	
	estimated 𝔼[x] = $(round(𝔼f, sigdigits=3)) (true: $(round(𝔼f_true, sigdigits=3)))

	👍

	we also show a weighted histogram where the particles are weighted according to the normalized importance weight. we get an empirical approximation of the desired beta distribution! 🙌
"

# ╔═╡ 43f81a70-0061-487c-b9d8-9d60efc72d96
begin
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="x", ylabel="density")
	hist!(ax, xᵢs, weights=w̃s, normalization=:pdf, bins=range(0, 1, length=30),
		  label="weighted particles", color=Cycled(3))
	xs = range(-0.5, 1.5, length=300)
	lines!(ax, xs, γ.(xs) / Z_true, label="π(x)", linewidth=3)
	lines!(ax, xs, η.(xs), label="η(x)", linewidth=3, linestyle=:dash)
	axislegend(ax)
	fig
end

# ╔═╡ 3e4b200a-999a-4d1a-abb4-543c4649c35c
hist(w̃s, axis=(;xlabel="normalized importance weight, w̃ᵢ", ylabel="# particles"))

# ╔═╡ c5fcad79-31ce-4d36-83c3-681f09d2de89
ess = 1 / sum(w̃s .^ 2)

# ╔═╡ 22e70101-17a6-4558-95a6-24002007e682
md"👓 a measure of the effective sample size (ESS) is

$${\rm ESS} = \frac{1}{\sum_{i=1}^N\tilde{w}_i^2}$$

which is close to $N$ if the weights are roughly the same (good!) and closer to one if one weight dominates (bad!).

if the particles sampled from the importance density $\eta(x)$ do not densely cover the high-density regions of the target density $\pi(x)$, then the ESS will be small, and $\eta(x)$ should be adjusted to be closer to the target density.

!!! example
	let's look at the distribution of weights and compute the ESS for our running example. there are lots of zero weights because many samples from our Gaussian fall outside the support of the beta distribution. the ESS is $(round(Int, ess)).
"

# ╔═╡ 23068adb-1c24-4344-b329-e2116579d036
md"
🤯 we used samples from a Gaussian distribution, and re-weighted them, to compute the normalizing factor of a beta distribution and its mean!

## references

Del Moral, Pierre, Arnaud Doucet, and Ajay Jasra. \"Sequential Monte Carlo Samplers.\" Journal of the Royal Statistical Society Series B: Statistical Methodology 68.3 (2006): 411-436.
"

# ╔═╡ Cell order:
# ╠═ab2c1074-cc27-11f1-83ad-358a865b7c3d
# ╟─5d47e1be-f9ac-4942-b3bf-1169a9ed80c7
# ╠═0c563ea7-c248-42b8-801f-d1876038b9a5
# ╠═9195342c-4bdb-46d0-aa75-4fe2b311f5dd
# ╠═5db6f3ef-b573-4637-b948-41a2b5c49069
# ╠═2e0eab51-e0ff-4556-929f-742709df277b
# ╟─d6601702-dc6b-4512-a73d-902d475a8f62
# ╠═f19d3174-008b-44de-8f77-89f1700ff5fc
# ╠═1b8e1902-5dd3-42ad-8599-74538c7d70b8
# ╟─3d24db8d-3841-47a1-9950-f9b52caec1af
# ╠═8f6a78a7-fed4-445e-943e-ac77c749e223
# ╠═b850cb61-b9ea-45df-b864-3a7fc27a49ab
# ╠═3bf89820-c75e-4a54-aba5-02d0decf9944
# ╠═c6d01c1f-a7a0-4904-bbd4-88fef7feac1c
# ╠═c36f6626-69fa-4be4-8f7d-7325caf66b53
# ╠═43f81a70-0061-487c-b9d8-9d60efc72d96
# ╟─22e70101-17a6-4558-95a6-24002007e682
# ╠═3e4b200a-999a-4d1a-abb4-543c4649c35c
# ╠═c5fcad79-31ce-4d36-83c3-681f09d2de89
# ╟─23068adb-1c24-4344-b329-e2116579d036
