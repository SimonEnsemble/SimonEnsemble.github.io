---
title: the shape of a chemical plume
date: 2024-05-28
tags: 
    - math 
---
_author_: Cory Simon

[WORK IN PROGRESS]

# setup

we wish to mathematically model the average, steady-state shape of a chemical plume whose source is a point and continuous.
* $R$ [g/min]: the constant rate at which the chemical is continuously released
* $n\in\{2,3\}$: dimension of the space $\mathbb{R}^n$ in which the chemical plume exists
* $\mathbf{x}_0\in\mathbb{R}^n$: [m] the location of the source of the chemical
* $\tau$ [min]: mean life-span of the chemical (owes to e.g., reaction with humidity or degradation under ultra-violet light)
* $\mathbf{v}$ [m/min]: mean wind vector
* $D$ [m$^2$/min]: diffusion coefficient (includes molecular diffusivity and [larger] turbulence diffusivity)

# the steady-state diffusion-advection-decay equation
a model for the average [over time] concentration of the chemical in the air, $c(\mathbf{x})$ for $\mathbf{x}\in\mathbb{R}^n$, is then the steady-state diffusion-advection-decay [partial differential] equation with a point-source term:
$$D {\boldsymbol \nabla}_{\mathbf{x}}^2 c - \mathbf{v} \cdot \{\boldsymbol \nabla}\_{\mathbf{x}}c - \tau^{-1}c + R \delta(\mathbf{x}-\mathbf{x}_0) = 0$$
respectively, the terms model isotropic diffusivity, advection by wind, decay, and introduction of the chemical into the environment. 

{{<figure
    src="/blog/plume/eqn.jpeg"
    caption="the model of the chemical plume."
>}}

# transformation to the modified Helmholtz equation
we transform the steady-state diffusion-advection-decay equation into a modified Helmholtz problem via the transformation:
$$c(\mathbf{x})=:u(\mathbf{x})e^{\mathbf{v}\cdot\mathbf{x} /(2D)}.$$

by the product rule:
$${\boldsymbol \nabla}_{\mathbf{x}}c= ({\boldsymbol \nabla}\_{\mathbf{x}}u(\mathbf{x}) )e^{\mathbf{v}\cdot\mathbf{x} /(2D)}+ u(\mathbf{x})\mathbf{v}(2D)^{-1}e^{\mathbf{v}\cdot\mathbf{x} /(2D)}.$$

and
$${\boldsymbol \nabla}_{\mathbf{x}}^2c = {\boldsymbol \nabla}\_{\mathbf{x}} \cdot {\boldsymbol \nabla}\_{\mathbf{x}}c=
({\boldsymbol \nabla}\_{\mathbf{x}}^2u(\mathbf{x}) )e^{\mathbf{v}\cdot\mathbf{x} /(2D)}+ ({\boldsymbol \nabla}\_{\mathbf{x}} u(\mathbf{x}))\mathbf{v} D^{-1}e^{\mathbf{v}\cdot\mathbf{x} /(2D)} + u(\mathbf{x})(2D)^{-2} \mathbf{v}\cdot \mathbf{v} e^{\mathbf{v}\cdot\mathbf{x} /(2D)}.
$$

stuffing these into the original diffusion-advection-decay eqn gives the modified Helmholtz problem in $u(\mathbf{x})$:

$${\boldsymbol \nabla}_{\mathbf{x}}^2u - \kappa^2 u = -\frac{R}{D} \delta(\mathbf{x}-\mathbf{x}_0) e^{-\mathbf{v}\cdot\mathbf{x}/(2D)} = -\frac{R}{D} \delta(\mathbf{x}-\mathbf{x}_0) e^{-\mathbf{v}\cdot\mathbf{x_0}/(2D)}$$

with

$$\kappa^2:= \frac{\lVert \mathbf{v} \rVert^2+4 \tau^{-1}D}{4D^2}.$$

we now focus on finding the solution to this problem in $u(\mathbf{x})$ from which $c(\mathbf{x})$ follows.

# Fourier transform of the modified Helmholtz equation

we take the Fourier transform of the modified Helmholtz problem:

$$-4 \pi^2 i D \lVert \boldsymbol \omega \rVert \tilde{u}(\mathbf{\omega}) + \kappa^2 \tilde{u}(\mathbf{\omega})= \frac{R}{D}e^{-2\pi i \boldsymbol\omega \cdot \mathbf{x}_0}e^{-\mathbf{v}\cdot\mathbf{x_0}/(2D)}$$

and solve for the solution in the frequency domain:

$$\tilde{u}(\mathbf{\omega}) = \dfrac{\frac{R}{D}e^{-\mathbf{v}\cdot\mathbf{x}_0}}{4\pi^2 \lVert \boldsymbol \omega \rVert^2 + \kappa^2}e^{-2\pi i \boldsymbol\omega \cdot \mathbf{x}_0}.$$

to go back to the time domain, we use the inverse Fourier transform:

$$u(\mathbf{x}) = \int_{\mathbb{R}^n} \dfrac{\frac{R}{D}e^{-\mathbf{v}\cdot\mathbf{x}_0 / (2D)}}{4\pi^2 \lvert \boldsymbol \omega \rvert^2 + \kappa^2}e^{-2\pi i \boldsymbol\omega \cdot \mathbf{x}_0}e^{2\pi i \boldsymbol \omega \cdot \mathbf{x}} d\boldsymbol\omega.$$

we account for the shift:

$$u(\mathbf{x}-\mathbf{x}_0) = \frac{R}{D}e^{-\mathbf{v}\cdot\mathbf{x}_0 / (2D)} \int\_{\mathbb{R}^n} 
\dfrac{
e^{2\pi i \boldsymbol \omega \cdot \mathbf{x}}
}{4\pi^2 \lvert \boldsymbol \omega \rvert^2 + \kappa^2}d\boldsymbol\omega.$$

the approach to evaluate this integral is qualitatively different depending on the dimension of the space, $n$.

## case $n=3$

we write the integral in $\boldsymbol \omega \in \mathbb{R}^3$ in spherical coordinates $(\rho, \phi, \theta)$ with $\rho=\lVert \boldsymbol\omega\rVert$ the radius, $\phi$ the polar angle, and $\theta$ the azimuthal angle. aligning the $z$-axis with the vector $\mathbf{x}$, we can write:
$$\mathbf{x}\cdot \boldsymbol \omega= \lVert \mathbf{x} \lVert \lVert \boldsymbol \omega \lVert \cos \phi$$
and the integral becomes (recall, the volume element in sphereical coordinates is $\rho^2\sin \phi d\rho d\phi d\theta$):

$$
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
\int_0^{\infty}
\frac{1}{i \lVert \mathbf{x} \rVert}
\left(e^{2\pi i \rho \lVert \mathbf{x}\rVert } - e^{-2\pi i \rho \lVert \mathbf{x}\rVert }\right)
\dfrac{\rho
}{4\pi^2 \rho^2 + \kappa^2}
d\rho
$$

breaking the integral into a sum of two integrals and doing a $\hat{\rho}:=-\rho$ substitution in the second allows us to write this as a single integral:
$$
\int_{-\infty}^{\infty}
\frac{1}{i \lVert \mathbf{x} \rVert}
e^{2\pi i \rho \lVert \mathbf{x}\rVert }
\dfrac{\rho
}{4\pi^2 \rho^2 + \kappa^2}
d\rho
$$

we tackle this integral using contour integration in the complex plane, with $\gamma \subset \mathcal{C}$ a half-circle that lies on the real line from $-R$ to $R$ followed by the upper half semicircle of radius $R$, centered at $0$. we choose $R> \kappa/(2\pi)$ so the semi-circle encloses the pole of the integrand.

{{<figure
    src="/blog/plume/contour.jpeg"
    caption="the model of the chemical plume."
>}}

$$
\frac{1}{i \lVert \mathbf{x} \rVert}
\oint_{\gamma}
e^{2\pi i \lVert \mathbf{x}\rVert z}
\dfrac{z
}{4\pi^2z^2  + \kappa^2}
dz
=:
\frac{1}{i \lVert \mathbf{x} \rVert}
\oint\_{\gamma}g(z)dz
$$

now, 
$$
\oint_{\gamma}f(z)dz = \int_{-R}^R f(z) dz + \int_{\\{Re^{i\theta} :\, \theta \in [0, \pi]\\}} f(z)dz
$$

via Jordan's lemma, $\int_{\\{Re^{i\theta} :\, \theta \in [0, \pi]\\}}f(z) dz \rightarrow 0$.

via Cauchy's residue theorem:
$$\oint_{\gamma}f(z)dz=2\pi i \text{Res}(f, \kappa/(2\pi)i)$$
since the contour $\gamma$ encircles the pole $\kappa/(2\pi)i$ of the complex function $f(z)$.

the residue is:
$$\text{Res}(f, \kappa/(2\pi)i)=\lim_{z\rightarrow \kappa/(2\pi)i}(z-\kappa/(2\pi)i)f(z)=\frac{1}{8\pi^2}e^{-\kappa \lVert \mathbf{x} \rVert}$$

finally, then:
$$
\int\_{\mathbb{R}^n} 
\dfrac{
e^{2\pi i \boldsymbol \omega \cdot \mathbf{x}}
}{4\pi^2 \lvert \boldsymbol \omega \rvert^2 + \kappa^2}d\boldsymbol\omega = 
\frac{1}{4\pi \lVert \mathbf{x} \rVert}e^{-\kappa \mathbf{x}}.
$$

and we have:

$$\boxed{u(\mathbf{x})=\frac{R}{4\pi D} \frac{1}{\lVert \mathbf{x}-\mathbf{x}_0\rVert} e^{-\kappa \lVert \mathbf{x}-\mathbf{x}_0\rVert}e^{\mathbf{v} \cdot (\mathbf{x}-\mathbf{x}\_0)/(2D)}}$$

## case $n=2$

TODO

## visualization

TODO

# appendix

## the Fourier transform

### definition
the Fourier transform of a function $f(\mathbf{x})$ is defined as:
$$\mathcal{F}[f(\mathbf{x})] ({\boldsymbol \omega}):=\int_{\mathbb{R}^n} f(\mathbf{x}) e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x} =: \tilde{f}({\boldsymbol \omega})$$
where ${\boldsymbol \omega}\in\mathbb{R}^n$.

### inverse Fourier transform
for recovering the function $f(\mathbf{x})$ from its Fourier transform $\tilde{f}(\mathbf{\omega})$, the inverse Fourier transform is:
$$f(\mathbf{x})=\mathcal{F}^{-1}[f({\boldsymbol \omega})] (\mathbf{x}) :=\int_{\mathbb{R}^n} \tilde{f}({\boldsymbol \omega}) e^{2\pi i {\boldsymbol \omega}}d\mathbf{x}.$$

### Fourier transform of a Dirac delta function

the Fourier transform of a Dirac delta function is:
$$\mathcal{F}[\delta(\mathbf{x}-\mathbf{x}_0)] ({\boldsymbol \omega}) := \int\_{\mathbb{R}^n} \delta (\mathbf{x} - \mathbf{x}_0) e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x} = e^{-2\pi i \mathbf{x}_0}.$$
via the sifting property of the delta function.

### Fourier transform of derivatives

the Fourier transform of a derivative is:
$$\mathcal{F}\left[\frac{\partial f}{\partial x_i}\right] ({\boldsymbol \omega}) := \int\_{\mathbb{R}^n} \frac{\partial f}{\partial x_i} e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x}.$$
tackling the integral over $x_i \in (-\infty, \infty)$ by parts, with $u:=e^{-2\pi i \boldsymbol \omega}$ and $dv=\frac{\partial f}{\partial x_i} dx_i$ gives:
$$\mathcal{F}\left[\frac{\partial f}{\partial x_i}\right] ({\boldsymbol \omega}) = 2 \pi i \omega_i \tilde{f}(\boldsymbol \omega).$$
hence,
$$\mathcal{F}\left[{\boldsymbol \nabla} f \right] ({\boldsymbol \omega}) = 2 \pi i {\boldsymbol \omega} \tilde{f}(\boldsymbol \omega)$$
and
$$\mathcal{F}\left[{\boldsymbol \nabla}^2 f \right] ({\boldsymbol \omega}) = -4 \pi^2 {\boldsymbol \omega} \cdot {\boldsymbol \omega} \tilde{f}(\boldsymbol \omega) = -4\pi^2 \lVert \boldsymbol \omega \rVert^2 \tilde{f}(\boldsymbol \omega).$$

### Fourier transform of a translated function

finally, we remark that a translation of the function $f(\mathbf{x})$ by vector $\mathbf{x}_0$ gives:
$$\mathcal{F}[f(\mathbf{x}-\mathbf{x}_0)] ({\boldsymbol \omega}) := \int\_{\mathbb{R}^n} f (\mathbf{x} - \mathbf{x}_0) e^{-2\pi i {\boldsymbol \omega}}d\mathbf{x} = e^{-2\pi i \mathbf{x}_0} \tilde{f}(\boldsymbol \omega).$$
shown after a substitution $\mathbf{x}^\prime:=\mathbf{x}-\mathbf{x}_0$. hence,
$$f(\mathbf{x})=\mathcal{F}^{-1}[f({\boldsymbol \omega})e^{-2\pi i \mathbf{x}_0}] (\mathbf{x}) = f(\mathbf{x}-\mathbf{x}_0).$$


# references
* [notes on the Fourier transform by Brad Osgood](https://see.stanford.edu/materials/lsoftaee261/chap8.pdf)
* Vergassola, M., Villermaux, E., & Shraiman, B. I. (2007). ‘Infotaxis’ as a strategy for searching without gradients. Nature, 445 (7126), 406-409.
* Yukawa Potential hw assignment from Dept. of Physics at Montana State University. [here](https://www.physics.montana.edu/avorontsov/teaching/phsx545/documents/problems07s.pdf)
