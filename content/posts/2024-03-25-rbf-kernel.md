---
title: the RBF kernel's feature map
tags:
    - math
date: 2024-03-25
---
_author_: Cory Simon

# overview
the radial basis function (RBF) kernel between two vectors $x\in\mathbb{R}^n$ and $x^\prime \in \mathbb{R}^n$ is:
$$k(x, x^\prime) := e^{-\lVert x - x^\prime \rVert / (2\gamma^2)}$$

_explicitly_ evaluating $k(x, x^\prime)$ is _implicitly_ equivalent to (i) mapping both these vectors into an infinite-dimensional Hilbert space with feature map $\phi:\mathbb{R}^n \rightarrow \mathbb{H}$ then (ii) taking the dot product of these resultant vectors. i.e.
$$k(x, x^\prime)=\phi(x) \cdot \phi(x^\prime).$$

here, we derive the feature map 
$$\phi(x)=e^{-\lVert x\rVert^2/(2\gamma^2)} \left(
\frac{1}{\sqrt{k_1!\cdots k_n!}}\frac{1}{\gamma^j} \prod_{i=1}^n x_i^{k_i}
\right)_{j=0,...,\infty;\\, k_1+\cdots +k_n = j ;\\, k_1, ..., k_n \geq 0}.$$

# derivation

first, note $x\cdot x = \lVert x\rVert^2$ and expand:
$$\lVert x - x^\prime \rVert^2 = (x-x^\prime) \cdot (x - x^\prime)=x\cdot x + x^\prime \cdot x^\prime -2 x \cdot x^\prime.$$
use this to rewrite the kernel:
$$k(x, x^\prime)=e^{-\lVert x\rVert^2/(2\gamma^2)}e^{-\lVert x^\prime\rVert^2/(2\gamma^2)}e^{x \cdot x^\prime/\gamma^2}$$

now do a Taylor expansion of the last exponential. 
\begin{align}
e^{x\cdot x^\prime / \gamma^2}&=\sum_{j=0}^\infty \dfrac{(x\cdot x^\prime / \gamma^2)^j}{j!} \\\\
&=\sum_{j=0}^\infty \dfrac{1}{\gamma^{j}} \dfrac{1}{\gamma^j}\dfrac{1}{j!} (x\cdot x^\prime)^j.
\end{align}
apply the multinomial formula to the $(x\cdot x^\prime)^j$ term written out in terms of the entries of the vectors:
$$(x\cdot x^\prime)^j = \left(\sum_{i=1}^n x_i x^\prime_i\right)^j=\sum_{\substack{k_1+\cdots+ k_n=j\\\\ k_1, ..., k_n \geq 0}} \dfrac{j!}{k_1! \cdots k_n!} \prod_{i=1}^n (x_ix^\prime_i)^{k_i}$$ 

going back to the kernel function (the $1/j!$ in the Taylor expansion cancels the $j!$ in the multinomial formula):
$$k(x, x^\prime)=\sum_{j=0}^\infty \sum_{\substack{k_1+\cdots+ k_n=j\\\\ k_1, ..., k_n \geq 0}} 
\left(e^{-\lVert x\rVert^2/(2\gamma^2)}\frac{1}{\gamma^j} \frac{1}{\sqrt{k_1!\cdots k_n!}} \prod_{i=1}^n x_i^{k_i}\right)
\left(e^{-\lVert x^\prime\rVert^2/(2\gamma^2)}\frac{1}{\gamma^j} \frac{1}{\sqrt{k_1!\cdots k_n!}} \prod_{i=1}^n x_i^{\prime k_i}\right)
$$

hey, that's just the dot product $\phi(x^\prime)\cdot \phi(x)$ with the feature map above!

the insight to factor out eg. $\gamma^{-2j}=\gamma^{-j}\gamma^{-j}$ was our aim to write the kernel clearly as a dot product:
$$k(x, x^\prime) = \sum_\ell \phi_\ell(x) \phi_\ell(x^\prime).$$
