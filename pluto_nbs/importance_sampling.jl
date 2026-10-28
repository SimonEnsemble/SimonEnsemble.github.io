### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ ab2c1074-cc27-11f1-83ad-358a865b7c3d
begin
	import Pkg; Pkg.activate(; io=devnull)
	using CairoMakie, Distributions, SpecialFunctions

	set_theme!(theme_light())
	update_theme!(fontsize=16)
end

# ╔═╡ 5d47e1be-f9ac-4942-b3bf-1169a9ed80c7
md"

# importance sampling

by Cory Simon

😱 suppose we are interested in working with a probability density $\pi(x)$ over a sample space $\Omega$, but we face two challenges. first, we can't directly draw samples from it. second, we only know its shape, i.e. we only know $\pi(x)$ up to a scaling factor:

$$\pi(x) =  \gamma(x)/Z, \quad Z:=\int_{\Omega} \gamma(x)dx$$

where $\gamma(x)$ is a function we can evaluate pointwise. the unknown normalizing factor $Z$ for $\gamma(x)$ is the value of an integral that (a) lacks a closed-form solution and (b) is high-dimensional, making it infeasible to compute with numerical quadrature.

🎯 specifically, suppose we want to estimate the expectation of a function $f(x)$ under the density $\pi(x)$:

$$\mathbb{E}_{x\sim \pi(x)}[f(x)]=\int_\Omega \pi(x) f(x) dx$$

as well as the normalizing factor $Z$.

!!! example
    as a toy example, suppose $\pi(x)$ is a Beta distribution with $\alpha=3$ and $\beta=6$ over $\Omega=[0, 1]$ and pretend we don't know the normalizing factor. we only know 
    
    $$\pi(x)\propto x^2(1-x)^5=\gamma(x)$$

    we want to know both the normalizing factor $Z$ (secret: the Beta function evaluated at $\alpha=3$ and $\beta=6$) and the expectation of $f(x)=x$ under the beta distribution which is just its mean (secret: $\alpha/(\alpha+\beta)$). (yes, we _could_ compute $Z$ accurately and easily with numerical quadrature, but this is a toy example!)

!!! note
    if we _could_ draw samples from $\pi(x)$ readily, we would just draw $N$ samples $\{x_i\}_{i=1}^N$ from it then approximate the expectation of $f(x)$ as:

    $$\mathbb{E}_{x\sim \pi(x)}[f(x)]\approx \frac{1}{N}\sum_{i=1}^N f(x_i).$$

    alas...
"

# ╔═╡ 870e04b0-e3b9-44e7-bd61-c21dd2fd010a
α, β = 3.0, 6.0 # params of Beta distribution

# ╔═╡ 0c563ea7-c248-42b8-801f-d1876038b9a5
# shape of the Beta distribution
function γ(x)
	if 0 ≤ x ≤ 1
		return x^(α - 1) * (1 - x)^(β - 1)
	else
		return 0.0
	end
end

# ╔═╡ 9195342c-4bdb-46d0-aa75-4fe2b311f5dd
# function we want to average over the Beta distribution
f(x) = x

# ╔═╡ 5db6f3ef-b573-4637-b948-41a2b5c49069
Z_true = beta(α, β) # secret normalizing factor. can we recover it?

# ╔═╡ 2e0eab51-e0ff-4556-929f-742709df277b
𝔼f_true = α / (α + β) # secret expectation of f(x). can we recover it?

# ╔═╡ d6601702-dc6b-4512-a73d-902d475a8f62
md"🔨 now, suppose we _do_ have an *importance density* $\eta(x)$ that we can draw samples from. crucially, $\pi(x) > 0 \implies \eta(x) > 0$ i.e. the support of $\eta(x)$ contains the support of $\pi(x)$. ideally, we engineer/craft $\eta(x)$ to be as close to $\pi(x)$ as we can make it.

let $\{x_i\}_{i=1}^N$ be a set of $N$ samples from the importance density $\eta(x)$, which we call _particles_.

!!! example \"Example (continued)\"
	we choose the importance distribution to be a Gaussian distribution
	$$\eta(x) = \mathcal{N}\! \left(x \mid \mu, \sigma^2\right) = \dfrac{1}{\sigma \sqrt{2\pi}} \exp\left(-\frac{1}{2}\left(\dfrac{x - \mu}{\sigma}\right)^2\right)$$
	which we can sample from easily.
	we set its mean $\mu=1/2$ and standard deviation $\sigma=1/3$ so the bulk of the density hopefully covers our Beta distribution and its mean is nearby what we think the mean of the Beta distribution is.
	then, we draw $N=5000$ samples from this importance distribution that we crafted, which we call particles.
"

# ╔═╡ 182e19ef-1049-43c7-8911-bb0a2375a65a
μ, σ = 1/2, 1/3

# ╔═╡ f19d3174-008b-44de-8f77-89f1700ff5fc
# importance density
η(x) = 1.0 / (σ * sqrt(2 * π)) * exp(-((x - μ) / σ)^2 / 2)

# ╔═╡ 1b8e1902-5dd3-42ad-8599-74538c7d70b8
begin
	N = 5000
	# particles
	xᵢs = rand(Normal(μ, σ), N)

	# viz dist'n of particles
	hist(xᵢs, axis=(; xlabel="x", ylabel="# samples", title="samples from η(x)"))
end

# ╔═╡ 8f6a78a7-fed4-445e-943e-ac77c749e223
# weight function
w(x) = γ(x) / η(x)

# ╔═╡ b850cb61-b9ea-45df-b864-3a7fc27a49ab
# weights on particles
ws = [w(xᵢ) for xᵢ in xᵢs]

# ╔═╡ 3bf89820-c75e-4a54-aba5-02d0decf9944
# normalized weights
w̃s = ws / sum(ws)

# ╔═╡ c6d01c1f-a7a0-4904-bbd4-88fef7feac1c
# approximation to the normalizing factor for γ(x) to turn it into the target density
Z = mean(ws)

# ╔═╡ c36f6626-69fa-4be4-8f7d-7325caf66b53
# approximation to the expectation of f(x) under the target density π(x)
𝔼f = sum(w̃s .* f.(xᵢs))

# ╔═╡ 3d24db8d-3841-47a1-9950-f9b52caec1af
md"

❓ how can we leverage these samples from $\eta(x)$---the particles---to compute the $Z$ and $\mathbb{E}_{x\sim \pi(x)}[f(x)]$ quantities we want?

💡 write them as expectations over $\eta(x)$ that we _do_ have samples from.

$$\mathbb{E}_{x\sim \pi(x)}[f(x)]=\int_\Omega \pi(x) f(x) \frac{\eta(x)}{\eta(x)} dx=\frac{1}{Z} \int_\Omega \eta(x) w(x) f(x) dx = \frac{1}{Z} \mathbb{E}_{x\sim \eta(x)}[w(x)f(x)] $$

with weight function

$$w(x)=\frac{\gamma(x) }{\eta(x)}$$

that we can readily evaluate. 

likewise,

$$Z=\int_\Omega \gamma(x) \frac{\eta(x)}{\eta(x)}dx=\int_\Omega \eta(x) w(x) dx = \mathbb{E}_{x\sim \eta(x)}[w(x)].$$


!!! note \"notation recap\"
	*  $\pi(x)$: target density, where we only know its shape
	*  $f(x)$: a function whose average we want to know, over the target density
	*  $\gamma(x)$: function we can evaluate point-wise, giving the shape of the target density
	*  $Z$: unknown normalizing factor for $\gamma(x)$ to turn it into the target density $\pi(x)$.
	*  $\eta(x)$: importance density that we construct and can sample from. ideally, it resembles the target density. necessarily, its support includes the support of the target density.
	*  $w(x)$: a weight function, which is higher (lower) in regions where the importance density underestimates (overestimates) the target density.

next, we approximate $\eta(x)$ using the samples we have from it---the particles $\{x_i\}_{i=1}^N$---as a sum of Dirac delta functions at those samples:

$$\eta(x)\approx \frac{1}{N}\sum_{i=1}^N\delta(x-x_i)$$

stuffing this approximation into the integrals above concerning expectations with respect to the importance density, we get the normalization factor as the average importance weight:

$$Z \approx \frac{1}{N}\sum_{i=1}^N w(x_i)$$

and the expectation of the function $f(x)$ as

$$\mathbb{E}_{x\sim \pi(x)}[f(x)] \approx \frac{1}{Z} \frac{1}{N}\sum_{i=1}^N w(x_i) f(x_i) \approx \sum_{i=1}^N \tilde{w}_i f(x_i)$$

where the *normalized importance weight* associated with particle $i$ is

$$\tilde{w}_i := w(x_i) / \sum_{j=1}^N w(x_j)$$

making the normalized importance weights sum to one.
intuitively, the normalized importance weights correct the samples from the importance density. where the importance density $\eta(x)$ over-estimates the target density $\pi(x)$, the weights will be smaller to down-weight those particles. where the importance density under-estimates the target density, the weights will be larger to up-weight those particles.

the estimate of $Z$ is unbiased while the estimate of $E_{x\sim \pi(x)}[f(x)]$ is biased but consistent.

!!! example \"Example (continued)\"
	below, we define the weight function, compute the importance weights, then use these as weights on the particles to compute the normalization factor and expectation of $x$ under our Beta distribution where we pretended its normalizing factor and mean were unknown!

	estimated Z = $(round(Z, sigdigits=3)) (true: $(round(Z_true, sigdigits=3))) 
	
	estimated 𝔼[x] = $(round(𝔼f, sigdigits=3)) (true: $(round(𝔼f_true, sigdigits=3)))

	👍

	we also show a weighted histogram where the particles are weighted according to the normalized importance weights. we get an empirical approximation of the desired Beta distribution! 🙌
"

# ╔═╡ 43f81a70-0061-487c-b9d8-9d60efc72d96
begin
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="x", ylabel="density")
	hist!(ax, xᵢs, weights=w̃s, normalization=:pdf, bins=range(0, 1, length=30),
		  label="weighted particles", color=Cycled(3))
	xs = range(-0.5, 1.5, length=300)
	lines!(ax, xs, γ.(xs) / Z_true, label="target density, π(x)", linewidth=3)
	lines!(ax, xs, η.(xs), label="importance density, η(x)", linewidth=3, linestyle=:dash)
	axislegend(ax)
	fig
end

# ╔═╡ 3e4b200a-999a-4d1a-abb4-543c4649c35c
hist(N * w̃s, axis=(; xlabel="N w̃ᵢ", ylabel="# particles"))

# ╔═╡ 68d23791-3366-4221-baf8-d70ce400a007
sum(w̃s .≈ 0.0) # many weights are exactly zero

# ╔═╡ c5fcad79-31ce-4d36-83c3-681f09d2de89
ess = 1 / sum(w̃s .^ 2)

# ╔═╡ 22e70101-17a6-4558-95a6-24002007e682
md"👓 a measure of the effective sample size (ESS) of our weighted particles for approximating the target density is:

$${\rm ESS} = \frac{1}{\sum_{i=1}^N\tilde{w}_i^2}$$

the ESS is 
* closer to $N$ (good!) if the weights are roughly the same
* closer to one (bad) if one weight dominates

the closer the ESS is to the number of particles $N$, the more closely our engineered importance density $\eta(x)$ approximates the target density $\pi(x)$. i.e., if the particles sampled from the importance density $\eta(x)$ do not densely cover the high-density regions of the target density $\pi(x)$, then the ESS will be small, and $\eta(x)$ should be adjusted to be closer to the target density.

!!! example \"Example (continued)\"
	we look at the distribution of the normalized weights and compute the ESS for our running example. there are lots of zero and near-zero weights due to both shape-mismatch between the importance and target density. the weight is exactly zero when a sample from the importance distribution falls outside of the support of the target distribution, which happens here. the ESS is $(round(Int, ess)) out of $(N) particles. not bad. 

	if you download the notebook, you can play with $\mu$ and $\sigma$ to get higher or lower ESSs depending on how well the Gaussian matches the Beta target distribution.
"

# ╔═╡ 23068adb-1c24-4344-b329-e2116579d036
md"

🤯 in our example, we used samples from a Gaussian distribution, and re-weighted them, to compute the normalizing factor of a _Beta_ distribution and its mean! pretty cool.

!!! note \"Application to Bayesian inference\"
	importance sampling is applicable to Bayesian inference of parameters $\theta$ given data $\mathcal{D}$ and a data-generating model $p(\mathcal{D}\mid \theta)$. the target density $\pi(\theta):=p(\theta \mid \mathcal{D})$ is the posterior. the shape of the posterior is given by Bayes' theorem as $\gamma(\theta):=p( \mathcal{D}\mid\theta) p(\theta)$---the product of the likelihood and prior $p(\theta)$. the normalizing factor is the evidence $Z:=p(\mathcal{D})$. note, choosing the importance density as the prior makes the weights equal to the likelihood.
	
## references

Del Moral, Pierre, Arnaud Doucet, and Ajay Jasra. \"Sequential Monte Carlo Samplers.\" Journal of the Royal Statistical Society Series B: Statistical Methodology 68.3 (2006): 411-436.

## thank you

to Qia Ke for helpful discussions and feedback on the first draft of this article.
"

# ╔═╡ Cell order:
# ╠═ab2c1074-cc27-11f1-83ad-358a865b7c3d
# ╟─5d47e1be-f9ac-4942-b3bf-1169a9ed80c7
# ╠═870e04b0-e3b9-44e7-bd61-c21dd2fd010a
# ╠═0c563ea7-c248-42b8-801f-d1876038b9a5
# ╠═9195342c-4bdb-46d0-aa75-4fe2b311f5dd
# ╠═5db6f3ef-b573-4637-b948-41a2b5c49069
# ╠═2e0eab51-e0ff-4556-929f-742709df277b
# ╟─d6601702-dc6b-4512-a73d-902d475a8f62
# ╠═182e19ef-1049-43c7-8911-bb0a2375a65a
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
# ╠═68d23791-3366-4221-baf8-d70ce400a007
# ╠═c5fcad79-31ce-4d36-83c3-681f09d2de89
# ╟─23068adb-1c24-4344-b329-e2116579d036
