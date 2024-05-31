---
title: modeling the shape of a chemical plume
date: 2024-05-28
tags: 
    - math 
---
_author_: Cory Simon

[WORK IN PROGRESS]

## modeling a chemical plume

we wish to mathematically model the average, steady-state shape of a chemical plume whose source is a point and continuous.

specifically, we wish to model the function $c(\mathbf{x})$ [g/m$^n$], the average concentration of a chemical in the air of a spatially homogeneous environment $\mathbb{R}^n$ ($n\in\\{2,3\\}$), with $\mathbf{x}$ a point in the environment $\mathbb{R}^n$.

we account for four pieces of physics:
1. the chemical is continuously released into the environment at a constant rate $R$ [g/min] from a point source at location $\mathbf{x}_0\in\mathbb{R}^n$.
2. wind transports the chemical downwind through advection, with $\mathbf{v}\in\mathbb{R}^n$ [m/min] the (constant) mean wind vector.
3. the chemical diffuses, owing to both molecular diffusivity and [dominant] turbulent diffusivity, with diffusion coefficient $D$ [m$^2$/min].
4. the chemical decays, owing to e.g. reaction with humidity or photodegradation (via ultraviolet radiation), with $\tau$ [min] the mean lifespan of the chemical.

## the steady-state diffusion-advection-decay equation
a model for the average [over time] concentration of the chemical in the air, $c(\mathbf{x})$ for $\mathbf{x}\in\mathbb{R}^n$, is then the steady-state diffusion-advection-decay [partial differential] equation with a point-source term:
$$D {\boldsymbol \nabla}_{\mathbf{x}}^2 c(\mathbf{x}) - \mathbf{v} \cdot \{\boldsymbol \nabla}\_{\mathbf{x}}c(\mathbf{x}) - \tau^{-1}c(\mathbf{x}) + R \delta(\mathbf{x}-\mathbf{x}_0) = 0$$
respectively, the terms model isotropic diffusivity, advection by wind, decay, and introduction of the chemical into the environment. (here, $\delta(\cdot)$ is the Dirac delta function.)

{{<figure
    src="/blog/plume/eqn.jpeg"
    caption="the model of the chemical plume at steady state."
>}}

### transformation to the modified Helmholtz equation
we transform the steady-state diffusion-advection-decay equation into a modified Helmholtz problem via the transformation:
$$c(\mathbf{x})=:u(\mathbf{x})e^{\mathbf{v}\cdot\mathbf{x} /(2D)}.$$

by the product rule:
$${\boldsymbol \nabla}\_{\mathbf{x}}c(\mathbf{x})= ({\boldsymbol \nabla}\_{\mathbf{x}}u(\mathbf{x}) )e^{\mathbf{v}\cdot\mathbf{x} /(2D)}+ u(\mathbf{x})\mathbf{v}(2D)^{-1}e^{\mathbf{v}\cdot\mathbf{x} /(2D)}.$$
and
$${\boldsymbol \nabla}_{\mathbf{x}}^2c(\mathbf{x}) = {\boldsymbol \nabla}\_{\mathbf{x}} \cdot {\boldsymbol \nabla}\_{\mathbf{x}}c (\mathbf{x})=
({\boldsymbol \nabla}\_{\mathbf{x}}^2u(\mathbf{x}) )e^{\mathbf{v}\cdot\mathbf{x} /(2D)}+ ({\boldsymbol \nabla}\_{\mathbf{x}} u(\mathbf{x}))\mathbf{v} D^{-1}e^{\mathbf{v}\cdot\mathbf{x} /(2D)} + u(\mathbf{x})(2D)^{-2} \mathbf{v}\cdot \mathbf{v} e^{\mathbf{v}\cdot\mathbf{x} /(2D)}.
$$

stuffing these into the original diffusion-advection-decay eqn gives the modified Helmholtz problem in $u(\mathbf{x})$ for $\mathbf{x}\in\mathbb{R}^n$:
$${\boldsymbol \nabla}_{\mathbf{x}}^2u - \kappa^2 u = -\frac{R}{D} \delta(\mathbf{x}-\mathbf{x}_0) e^{-\mathbf{v}\cdot\mathbf{x}/(2D)} = -\frac{R}{D} \delta(\mathbf{x}-\mathbf{x}_0) e^{-\mathbf{v}\cdot\mathbf{x_0}/(2D)}$$
with
$$\kappa^2:= \frac{\lVert \mathbf{v} \rVert^2+4 \tau^{-1}D}{4D^2}.$$
the latter equality holds because the Dirac delta flattens the exponential term except for at $\mathbf{x}_0$.

we now focus on finding the solution $u(\mathbf{x})$ to this modified Helmholtz problem, from which $c(\mathbf{x})$ follows.

## solution 
### the Fourier transform of the modified Helmholtz equation

we take the Fourier transform of the modified Helmholtz problem, giving an algebraic equation in the Fourier transform of $u(\mathbf{x})$, $\tilde{u}(\boldsymbol \omega)$:
$$-4 \pi^2\lVert \boldsymbol \omega \rVert^2 \tilde{u}(\mathbf{\omega}) - \kappa^2 \tilde{u}(\mathbf{\omega})= -\frac{R}{D}e^{-2\pi i \boldsymbol\omega \cdot \mathbf{x}_0}e^{-\mathbf{v}\cdot\mathbf{x_0}/(2D)}$$

and solve for the solution in the frequency domain:

$$\tilde{u}(\mathbf{\omega}) = \dfrac{\frac{R}{D}e^{-\mathbf{v}\cdot\mathbf{x}_0 / (2D)}}{4\pi^2 \lVert \boldsymbol \omega \rVert^2 + \kappa^2}e^{-2\pi i \boldsymbol\omega \cdot \mathbf{x}_0}.$$

### the inverse Fourier transform

to obtain $u(\mathbf{x})$, we take the inverse Fourier transform of $\tilde{u}(\boldsymbol \omega)$:

$$u(\mathbf{x}) = \int_{\mathbb{R}^n} \dfrac{\frac{R}{D}e^{-\mathbf{v}\cdot\mathbf{x}_0 / (2D)}}{4\pi^2 \lvert \boldsymbol \omega \rvert^2 + \kappa^2}e^{-2\pi i \boldsymbol\omega \cdot \mathbf{x}_0}e^{2\pi i \boldsymbol \omega \cdot \mathbf{x}} d\boldsymbol\omega.$$

we account for the shift:

$$u(\mathbf{x}-\mathbf{x}_0) = \frac{R}{D}e^{-\mathbf{v}\cdot\mathbf{x}_0 / (2D)} 
\int\_{\mathbb{R}^n} 
\dfrac{
e^{2\pi i \boldsymbol \omega \cdot \mathbf{x}}
}{4\pi^2 \lvert \boldsymbol \omega \rvert^2 + \kappa^2}d\boldsymbol\omega.$$

the approach to evaluate the integral 

$$I(\mathbf{x}):=
\int\_{\mathbb{R}^n} 
\dfrac{
e^{2\pi i \boldsymbol \omega \cdot \mathbf{x}}
}{4\pi^2 \lvert \boldsymbol \omega \rvert^2 + \kappa^2}d\boldsymbol\omega$$

is qualitatively different depending on the dimension of the space, $n$. 

#### case $n=3$

we write the integral in $\boldsymbol \omega \in \mathbb{R}^3$ in spherical coordinates $(\rho, \phi, \theta)$ with $\rho=\lVert \boldsymbol\omega\rVert$ the radius, $\phi$ the polar angle, and $\theta$ the azimuthal angle. aligning the $z$-axis with the vector $\mathbf{x}$, we can write:
$$\mathbf{x}\cdot \boldsymbol \omega= \lVert \mathbf{x} \lVert \lVert \boldsymbol \omega \lVert \cos \phi$$
and the integral becomes (recall, the volume element in spherical coordinates is $\rho^2\sin \phi d\rho d\phi d\theta$):

$$
I(\mathbf{x})=
\int_0^{2\pi}
\int_0^{\pi}
\int_0^{\infty}
\dfrac{
e^{2\pi i \rho \lVert \mathbf{x}\rVert \cos \phi }
}{4\pi^2 \rho^2 + \kappa^2}
\rho^2 \sin \phi 
d\rho
d\phi
d\theta
$$

right away, we can knock out the integral in $\theta$ as $2\pi$ and tackle the integral in $\phi$ with a simple substitution $\xi := \cos \phi$ giving:
$$
I(\mathbf{x})=
\int_0^{\infty}
\frac{1}{i \lVert \mathbf{x} \rVert}
\left(e^{2\pi i \rho \lVert \mathbf{x}\rVert } - e^{-2\pi i \rho \lVert \mathbf{x}\rVert }\right)
\dfrac{\rho
}{4\pi^2 \rho^2 + \kappa^2}
d\rho
$$

breaking the integral into a sum of two integrals 
$$
I(\mathbf{x})=
\frac{1}{i \lVert \mathbf{x} \rVert}\left(
\int_0^{\infty}
e^{2\pi i \rho \lVert \mathbf{x}\rVert }
\dfrac{\rho
}{4\pi^2 \rho^2 + \kappa^2}
d\rho-
\int_0^{\infty}
e^{-2\pi i \rho \lVert \mathbf{x}\rVert }
\dfrac{\rho
}{4\pi^2 \rho^2 + \kappa^2}
d\rho
\right)
$$

and doing a $\hat{\rho}:=-\rho$ substitution in the second integral allows us to write this as a single integral:
$$
I(\mathbf{x})=
\frac{1}{i \lVert \mathbf{x} \rVert}
\int_{-\infty}^{\infty}
e^{2\pi i \rho \lVert \mathbf{x}\rVert }
\dfrac{\rho
}{4\pi^2 \rho^2 + \kappa^2}
d\rho.
$$

we tackle this integral by involving it in a contour integration in the complex plane $\mathcal{C}$, with the contour $\gamma(R) \subset \mathcal{C}$ a half-circle that lies on the real line from $-R$ to $R$ followed by the upper semicircle of radius $R$, centered at $0$. we choose $R> \kappa/(2\pi)$ so the contour $\gamma$ encloses the pole (one of two) $z=\kappa/(2\pi)i$ of the integrand 

$$g(z):=e^{2\pi i z \lVert \mathbf{x}\rVert }
\dfrac{z
}{4\pi^2 z^2 + \kappa^2}\, z \in \mathbb{C}.$$

{{<figure
    src="/blog/plume/contour.jpeg"
    caption="the contour $\gamma$ together with the two poles of $g(z)$."
>}}

then, we can write:
$$
\frac{1}{i \lVert \mathbf{x} \rVert}
\oint_{\gamma(R)}g(z)dz = 
\frac{1}{i \lVert \mathbf{x} \rVert}
\int_{-R}^R g(z) dz +
\frac{1}{i \lVert \mathbf{x} \rVert} \int_{\\{Re^{i\theta} :\, \theta \in [0, \pi]\\}} g(z)dz
$$

now, as $R\rightarrow\infty$, the second integral becomes $I(\mathbf{x})$ we're looking for. and,
via Jordan's lemma, $\int_{\\{Re^{i\theta} :\, \theta \in [0, \pi]\\}}g(z) dz \rightarrow 0$ as $R\rightarrow \infty$.
so:
$$
\frac{1}{i \lVert \mathbf{x} \rVert}
\lim_{R\rightarrow\infty} \oint_{\gamma(R)}g(z)dz = 
I(\mathbf{x})
$$

via Cauchy's residue theorem:
$$\oint_{\gamma}g(z)dz=2\pi i \text{Res}(g, \kappa/(2\pi)i)$$
since the contour $\gamma$ encircles the pole $\kappa/(2\pi)i$ of the complex function $g(z)$ in a counter-clockwise fashion.

calculating the residue:
$$\text{Res}(g, \kappa/(2\pi)i)=\lim_{z\rightarrow \kappa/(2\pi)i}(z-\kappa/(2\pi)i)g(z)=\frac{1}{8\pi^2}e^{-\kappa \lVert \mathbf{x} \rVert}$$

finally, then:
$$
I(\mathbf{x})=
\int\_{\mathbb{R}^n} 
\dfrac{
e^{2\pi i \boldsymbol \omega \cdot \mathbf{x}}
}{4\pi^2 \lvert \boldsymbol \omega \rvert^2 + \kappa^2}d\boldsymbol\omega = 
\frac{1}{4\pi \lVert \mathbf{x} \rVert}e^{-\kappa\lVert \mathbf{x}\rVert}.
$$

and we have:

$$\boxed{u(\mathbf{x})=\frac{R}{4\pi D} \frac{1}{\lVert \mathbf{x}-\mathbf{x}_0\rVert} e^{-\kappa \lVert \mathbf{x}-\mathbf{x}_0\rVert}e^{-\mathbf{v} \mathbf{x}\_0/(2D)}}$$

going back to $c(\mathbf{x})$, finally, we have the shape of the chemical plume in 3D!

$$\boxed{c(\mathbf{x})=\frac{R}{4\pi D} \frac{1}{\lVert \mathbf{x}-\mathbf{x}_0\rVert} e^{-\kappa \lVert \mathbf{x}-\mathbf{x}_0\rVert}e^{\mathbf{v} \cdot (\mathbf{x}-\mathbf{x}\_0)/(2D)}}$$

### case $n=2$

TODO

### visualization

TODO

## appendix

### the Fourier transform

#### definition
the Fourier transform of a function $f(\mathbf{x})$ is defined as:
$$\mathcal{F}[f(\mathbf{x})] ({\boldsymbol \omega}):=\int_{\mathbb{R}^n} f(\mathbf{x}) e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x} =: \tilde{f}({\boldsymbol \omega})$$
where ${\boldsymbol \omega}\in\mathbb{R}^n$.

#### inverse Fourier transform
for recovering the function $f(\mathbf{x})$ from its Fourier transform $\tilde{f}(\boldsymbol \omega)$, the inverse Fourier transform is:
$$f(\mathbf{x})=\mathcal{F}^{-1}[f({\boldsymbol \omega})] (\mathbf{x}) :=\int_{\mathbb{R}^n} \tilde{f}({\boldsymbol \omega}) e^{2\pi i {\boldsymbol \omega}}d\mathbf{x}.$$

#### Fourier transform of a Dirac delta function

via the sifting property of the Dirac delta function
$$\mathcal{F}[\delta(\mathbf{x}-\mathbf{x}_0)] ({\boldsymbol \omega}) := \int\_{\mathbb{R}^n} \delta (\mathbf{x} - \mathbf{x}_0) e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x} = e^{-2\pi i \mathbf{x}_0}.$$

#### Fourier transform of derivatives

the Fourier transform of a derivative is:
$$\mathcal{F}\left[\frac{\partial f}{\partial x_i}\right] ({\boldsymbol \omega}) := \int\_{\mathbb{R}^n} \frac{\partial f}{\partial x_i} e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x}.$$
tackling the integral over $x_i \in (-\infty, \infty)$ by parts, with $u:=e^{-2\pi i \boldsymbol \omega}$ and $dv=\frac{\partial f}{\partial x_i} dx_i$ gives, provided $f \rightarrow 0$ as $x_i \rightarrow \pm \infty$:
$$\mathcal{F}\left[\frac{\partial f}{\partial x_i}\right] ({\boldsymbol \omega}) = 2 \pi i \omega_i \tilde{f}(\boldsymbol \omega).$$
hence,
$$\mathcal{F}\left[{\boldsymbol \nabla}_{\mathbf{x}} f \right] ({\boldsymbol \omega}) = 2 \pi i {\boldsymbol \omega} \tilde{f}(\boldsymbol \omega)$$
and
$$\mathcal{F}\left[{\boldsymbol \nabla}\_\mathbf{x}^2 f \right] ({\boldsymbol \omega}) = -4 \pi^2 {\boldsymbol \omega} \cdot {\boldsymbol \omega} \tilde{f}(\boldsymbol \omega) = -4\pi^2 \lVert \boldsymbol \omega \rVert^2 \tilde{f}(\boldsymbol \omega).$$

#### Fourier transform of a translated function

finally, we remark that a translation of the function $f(\mathbf{x})$ by vector $\mathbf{x}_0$ gives:
$$\mathcal{F}[f(\mathbf{x}-\mathbf{x}_0)] ({\boldsymbol \omega}) := \int\_{\mathbb{R}^n} f (\mathbf{x} - \mathbf{x}_0) e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x} = e^{-2\pi i \mathbf{x}_0} \tilde{f}(\boldsymbol \omega).$$
shown after a substitution $\mathbf{x}^\prime:=\mathbf{x}-\mathbf{x}_0$. hence,
$$\mathcal{F}^{-1}[f({\boldsymbol \omega})e^{-2\pi i \mathbf{x}_0}] (\mathbf{x}) = f(\mathbf{x}-\mathbf{x}_0).$$

## references
* notes on the Fourier transform by Brad Osgood [here](https://see.stanford.edu/materials/lsoftaee261/chap8.pdf)
* Vergassola, M., Villermaux, E., & Shraiman, B. I. (2007). ‘Infotaxis’ as a strategy for searching without gradients. Nature, 445 (7126), 406-409.
* Yukawa Potential hw assignment from Dept. of Physics at Montana State University. [here](https://www.physics.montana.edu/avorontsov/teaching/phsx545/documents/problems07s.pdf)
