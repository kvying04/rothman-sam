---
title: Math of W76 Mapping
author: Kevin Ying
date: 12 March 2026
---
# 0. Motivation
A main issue with dealing with the *w76* mutation currently is that we cannot confirm certain genotype/phenotype pairs when propagating from the *F2* to *F6* generation. The two driving factors underpinning the relationship between the genotype at *F2* and the phenotype at *F6* is the combination of whether the *w76* trait is dominant or recessive, and whether *w76* is a maternal effect mutation. 

\pagebreak

# 1. Math Behind Propagation
Assuming a "self" from a heterozygous animal (`+/-`) at generation $g=0$, the proportion of the total progeny coming from that heterozygous animal at $g=g_{\text{final}}$ which is homozygous one-way (either `+/+` or `-/-`) due to propagation at generation $g=k$ for all $1 \le k \le g_\text{final}$ can be modeled by: 
$$p(k) = \frac{1}{4}(\frac{1}{2})^{k-1} = \frac{1}{2^{k+1}}$$

So, for example, the proportion of total progeny that will be homozygous `+/+` that reached homozygosity at generation $g=1$ (i.e. after one generation of propagation) will be $p(1) = \frac{1}{4}(\frac{1}{2})^0 = \frac{1}{4}$. This means that $25\%$ of the total population coming from the original heterozygous individual at $g=0$ at $g_{\text{final}}$ is homozygous solely due to the propagation at the first generation. 

This means that the total proportion of animals that must be homozygous one-way (either `+/+` or `-/-`) after $g$ cycles of propagation must be the summation of all the products of each generation:
$$h(g)=\sum_{k=1}^g p(k)=\sum_{k=1}^g \frac{1}{2^{k+1}}$$

For example, given 4 cycles of propagation ($g=4$):
$$h(4)=\sum_{k=1}^4 p(k) =\frac{1}{2^2}+\frac{1}{2^3}+\frac{1}{2^4}+\frac{1}{2^5}=0.46875$$

In our application, starting from *F2* (where $50\%$ of the population will be heterozygous `+/-`), if we allow it to propagate to *F6* (which means propagating for $g=4$ cycles), we get that $46.875\%$ of the original proportion of heterozygous animals ($50\%$ in our case) will be homozygous one-way. As such, $46.875\% \times 50\% = 23.4375\%$ of the total number of progeny of *F2* (again, including the progeny from the $50\%$ of non-heterozygous animals) will be homozygous one-way. This means around $47\%$ of animals are homozygous in general, being either `+/+` or `-/-`. 

\pagebreak

# 2. Expected Values
## 2.1 Suppose Maternal Effect
We can use the functions outlined in part 1 (i.e. $p(*)$ and $h(*)$) to see how the results respond when considering *w76* as a maternal effect mutation. 

If we assume that the maternal effect is true, the phenotypical distribution at $g=g_\text{final}$ should be expected to be reflective of the genotypical distribution at $g=0$. Therefore, there are two cases to take into account: if *w76* is dominant or recessive. Doing the math for our case (individuals at *F2* can either be `+/-` or `-/-`):
$$\begin{align*}\text{proportion showing loss}_{\text{rec}}^\text{mat}&=\text{prop}_\text{loss,hom,F2} = 50\% \\ \text{proportion showing loss}_{\text{dom}}^\text{mat}&=\text{prop}_\text{loss,hom,F2} + \text{prop}_\text{loss,het,F2} \\ &= 50\% + 50\% = 100\%\end{align*}$$

<!-- We can use the functions outlined in part 1 (namely $p(*)$ and $h(*)$) to see how the results respond to whether *w76* is a maternal effect mutation. 

Suppose the case where the maternal effect is true. Then, if a worm who has a genotype $n$ has progeny, all progeny will also show the phenotype of $n$. In our case, $n$ stands for the ability for loss of *uaDf5*. This means that we can safely use the function $h(*)$, which models the amount of homozygous one-way progeny that appear, as synonymous to the proportion of progeny that will show a phenotype of loss. 

Note a small but important distinction: the application of the function $h(*)$ is dependent on whether the genotype $n$ is dominant or recessive:
$$\begin{align*}\text{proportion showing loss}_{\text{rec}}^\text{mat}&=\text{prop}_\text{loss,hom} = h(g) \\ \text{proportion showing loss}_{\text{dom}}^\text{mat}&=\text{prop}_\text{loss,hom} + \text{prop}_\text{loss,het} = h(g) + (1-2h(g)) \end{align*}$$

Applying our math:
$$\begin{align*}\text{proportion showing loss}_{\text{rec}}^\text{mat}&=\text{prop}_\text{loss,hom} = h(g) = 0.46875 = 46.875\% \\ \text{proportion showing loss}_{\text{dom}}^\text{mat}&=\text{prop}_\text{loss,hom} + \text{prop}_\text{loss,het} = h(g) + (1-2h(g)) = 0.46875 + (1-2(0.46875))=0.53125 = 53.125\% \end{align*}$$ -->

## 2.2 Do Not Suppose Maternal Effect
Now suppose that the maternal effect is not true. Then, one fact about *w76* "kicks in": the mutation has been observed to behave such that once a loss genotype has been established in generation $g=0$, it takes four generations of "self-ing" (until $g=4$) to be able to see a phenotypic response to that loss genotype. Frustratingly, this means that there may exist mismatch between the (unseen) genotype of the ancestor at $g=0$ and the (measurable) phenotype of the progeny at $g=4$. Under this assumption, we cannot use $h(*)$ and assume association between expected proportion of progeny containing some genotype and the measurable proportion of progeny showing the loss phenotype. 

However, since we expect to observe loss only from lines that have kept the loss genotype (`-/-` for recessive; `+/-` or `-/-` for dominant) for four generations, we can use the function $p(*)$ to give an expected proportion of loss. 

Let's first do the calculations assuming that *w76* is recessive. Then, we know that, for $g_{\text{final}} = 4$, certain loss-genotype progeny will not actually show loss. Specifically, those who only get the genotype at $g=2,3,4$ (corresponding to those who become homozygous-loss for the first time at the *F4*, *F5* or *F6* generations) will not actually show loss. Therefore, only individuals from lines who kept the loss genotype from *F2* to *F3* ($g=1$) would show loss of *uaDf5*. 

Mathematically, $p(1) = \frac{1}{4}$. Therefore, we may expect exactly $25\%$ of progeny from heterozygotic ancestors to show loss:[^1]
$$\begin{align*}\text{proportion showing loss}^\text{mat-}_\text{rec}&=\text{prop}_\text{loss,hom,F6} + \text{prop}_\text{loss,hom,F2}\\&= p(1)\cdot 0.5 + 50\% \\ &= 62.5\%\end{align*}$$

Let's now do the same calculation assuming that *w76* is dominant. Then, the total proportion of progeny from heterozygous ancestors that should show loss should be equal to the proportion of progeny that get the homozygous loss phenotype by $g=1$ plus the proportion of progeny that retain the heterozygous genotype by the end of $g=4$. Therefore:[^4]
$$\begin{align*}\text{proportion showing loss}^\text{mat-}_\text{dom}&=\text{prop}_\text{loss,hom,F2} + \text{prop}_\text{loss,hom,F6} + \text{prop}_\text{loss,het,F6}\\ &= 0.5 + [p(1) + (1-2h(4))] \cdot 0.5 \\&= 0.5 + [0.25 + (1-2(0.46875))]\cdot 0.5 \\ &= 0.65625= 65.625\%\end{align*}$$

## 2.3 Summary of Expected Values
In sum, the following is true, where the percentage correlates to the proportion of worms that should show loss at *F6*:

| | $\textbf{dom}$ | $\textbf{rec}$ |
| - | :-: | :-: |
| $\textbf{mat}^+$ | $100\%$ | $50\%$ |
| $\textbf{mat}^-$ | $65.625\%$ | $62.5\%$ |

\pagebreak

# 3. Addressing Maternal Effect: Variance
## 3.1 Context: Sample Sizes
Theoretically, given a large enough $n$, the true nature of *w76* can be revealed just by performing the above crossing experiment. However, to be able to prove that a proportion is statistically significant (and "different enough" from the other percentages), the following are the $n$ values needed to differentiate each percentage pair at $\alpha = 0.05$ and $\beta = 0.2$ (via a two-sided two-sample z-test for sample proportion):

| | $\textbf{dom,mat}^+$ | $\textbf{rec,mat}^+$ | $\textbf{dom,mat}^-$ | $\textbf{rec,mat}^-$ |
| - | :-: | :-: | :-: | :-: |
| $\textbf{dom,mat}^+$ | - | $11$ | $18$ | $16$ |
| $\textbf{rec,mat}^+$ | $11$ | - | $156$ | $247$ |
| $\textbf{dom,mat}^-$ | $18$ | $156$ | - | $3700$ |
| $\textbf{rec,mat}^-$ | $16$ | $247$ | $3700$ | - |

Evidently, some of these sample size requirements make it unfeasible to seek to entirely identify the characteristics of *w76* purely through the crossing experiment, since it would imply needing to do $n$ trials of the cross: for example, to be able to tell whether a result is from a recessive non-maternal characteristic or a dominant non-maternal characteristic, one would have to perform $n=3700$ trials (e.g. of 20 worms each).

Since the $n$ range is not exceedingly unreasonable, this could be argument to do line-wise lysate scoring rather than the gel experiment. This way, we would be able to get percentages, not boolean values, for *uaDf5* loss. This way, instead of each $n_i$ being a trial of $x$ worms, one could raise and score $n$ lines from *F2* hermaphrodites, theoretically cutting down the work substantially. 

## 3.2 Propagation Variance: Assuming Maternal Effect
Here, we will attempt to analyze the variance of the distribution of loss-phenotype progeny in each of the combinations listed above (e.g. $\text{dom,mat}^+$). Let's first assume that there does exist the maternal effect in *w76*. 

$\textbf{dom,mat}^+$: Since the maternal effect dictates that progenitive genotype has no impact on the shown phenotype, we know that the only variance that can exist in this case must be from the *F1* to *F2* cycle, referred to from now as $\sigma^2_{\text{mat}}$. In this case however, since we know that heterozygous (`+/-`) and homozygous-loss (`-/-`) are the only two genotypes, and since both convey a loss phenotype, we know that $\sigma^2_{\text{mat,dom}}=0$. 

$\textbf{rec,mat}^+$: The recessive maternal case is similar to the dominant maternal case, except for that, instead of having $100\%$ of its genotypic distribution at *F2* be conveying loss, only $50\%$ should in this case. Theoretically at *F2* however, there still should be a guaranteed even split between heterozygous and homozygous-loss genotype individuals. 

## 3.3 Propagation Variance: Multitype Galton-Watson Process
The difficulty of considering non-maternal cases is that the simple solution of taking the terminal expected genotypic distribution (e.g. $(0.46875,0.46875,0.0625)$ for *F6*) and calculating variance under a multinomial distribution does not work. This is because the self-ing process is skewed: a multinomial distribution expects that all progeny have the same "access" to the expected genotypic distribution, but this is not the case. Individuals that come from homozygous-loss *F5* mother will have a different distribution to "choose from" when compared to an individual from a heterozygous *F5* mother, for example.

Instead, we must approach this using a stochastic process known as a Galton-Watson branching process. Specifically, it will be a multitype Galton-Watson process (MGW), which allows us to model the development of a multitype population (e.g. types like `+/+`, `+/-`, `-/-`) when reproduction follows probabilistic law (e.g. heterozygous self-ing leads to a $(0.25, 0.5, 0.25)$ distribution). 

First, we must establish a progenital count distribution (i.e. how many "children" does each individual produce) for the MGW. Let's set $\pi$ to be the random variable representing the count of progeny from a worm, which we will model as a normal distribution:
$$\pi \sim N(\mu = 250, \sigma^2 = 35)$$
Again, this can be adjusted depending on observed brood behavior. 

Then, because we are tracking multiple types (namely homozygous-keep, heterozygous, and homozygous-loss), we can represent the expected distribution of an individual's progeny as a Markov matrix:
$$M = \begin{pmatrix}1 & 0 & 0 \\ 0.25 & 0.5 & 0.25 \\ 0 & 0 & 1\end{pmatrix}$$


We will use this simplified matrix first before considering *w76*'s four-generation property to prove that the MGW algorithm works as expected. 

The propagation matrix when adjusted for expected progeny growth ($\mathbb E [\pi] = \mu_\pi = 250$) will be the following:
$$M_{\text{prop}} = M\cdot \mu_\pi = \begin{pmatrix}250 & 0 & 0 \\ 62.5 & 125 & 62.5 \\ 0 & 0 & 250\end{pmatrix}$$

For example, given that an individual worm is heterozygous (second row), then its progeny will be expected to be distributed such that there are 62.5 homozygous-keep progeny, 62.5 homozygous-loss progeny, and 125 heterozygous progeny. 

Then, let's set an initial condition (i.e. how many individuals we want to start out with at generation $g=0$): arbitrarily, we can set the population at $g=0$ to be
$$N^{(0)}=20$$

The above is the total population $N$. We can also define the distribution of the individuals at generation $g=0$, given its initial proportion distribution $\vec{\textbf{p}}^{(0)}$, to be: 
$$\vec{\textbf{v}}^{(0)} = \vec{\textbf{p}}^{(0)} \cdot N^{(0)} = [0.5,0.5,0] \cdot 20= [10, 10, 0]$$

This states that there are 10 homozygous-loss, 10 heterozygous, and 0 homozygous-keep worms at generation 0 (in our case, *F2*). Notice that $N$ will always be equal to the sum of elements in $\vec{\textbf{v}}$. This means, conversely, that $\vec{\textbf{v}} \div N$ will give the proportional distribution, not the population distribution, of the individuals at some generartion. 

The expected number of individuals at any generation $g$ can be modeled using the propagation matrix $M_\text{prop}$ and the initial distribution vector:
$$\vec{\textbf{v}}^{(g)} = \vec{\textbf{v}}^{(0)}M_\text{prop}^g$$

Then, from what we have established, 
$$\vec{\textbf{p}}^{(g)} = \vec{\textbf{v}}^{(g)} \div N^{(g)} = \vec{\textbf{v}}^{(0)}M^g_{\text{prop}} \div N^{(g)}$$

Let's evaluate $\vec{\textbf{p}}^{(g)}$ at $g=4$ (our proportion distribution at *F6*) with initial distribution $[0,20,0]$ to confirm that our derivation was correct:
$$\begin{align*}\vec{\textbf{p}}^{(4)} = \begin{pmatrix}0 \\ 10 \\ 0\end{pmatrix} \begin{pmatrix}250 & 0 & 0 \\ 62.5 & 125 & 62.5 \\ 0 & 0 & 250\end{pmatrix}^4 \cdot \frac{1}{N^{(g)}} \\ = [0.46875, 0.0625, 0.46875] \end{align*}$$

This distribution looks familiar!

We can now find the variance of this system. Let's consider, for the sake of simplicity, that the variance of *F2* (i.e. is *F2* truly always $50\%$ homozygous?) is 0. That is, it is a guaranteed $50\%$ split always.[^2]

Under this, the only variance is the variance stemming from propagation for type $i$, or $\sigma^2_{\text{prop}, i}$, coming from the MGW:
$$\text{Var}(i) = \frac{p^2_i}{N^{(0)}}(\frac{\sigma^2_\pi}{\mu^2_\pi-\mu_\pi})$$

This is why this was so important to establish a distribution for brood size! Further, note that all $p_i$ are just elements in $\vec{\textbf{p}}^{(g)}$.

Now we can look at statistical properties that differentiate each phenotype category (e.g. $\textbf{dom,mat}^-$). Let's assume that *F2* has $50\%$ homozygous-loss and $50\%$ heterozygous individuals. Again, this comes with the caveat that the four-generation requirement for propagation has not been integrated into any of the calculations.

$\textbf{dom,mat}^-$: The expected proportion showing *uaDf5* loss assuming that *w76* is dominant will be 
$$\mathbb E[p_\text{dom}] = p_{\text{hom,loss}} + p_{\text{het}} = 0.734 + 0.031 = 76.5\%$$
The variance will be the process of recalculating $\text{Var}(i)$ for $i = p_{\text{dom}} = p_{\text{het}} + p_{\text{hom,loss}}$ instead of any individual phenotype; therefore, $\text{Var}(i_\text{dom}) \approx 0.00001516 = 1.516 \times 10^{-5}$. 

$\textbf{rec,mat}^-$: The expected proportion will be: 
$$\mathbb E[p_\text{dom}] = p_{\text{hom,loss}}= 0.734 = 73.4\%$$
Lastly, this distribution's variance will be $\text{Var}(i_\text{dom}) \approx 0.000001545 = 1.545 \times 10^{-6}$. 

Note that the two expected proportions calculated match the naive distribution of *F6* when assuming that the four-generation requirement does not exist.[^3]

## 3.4 Considering *w76*'s Four-Generation Property
The following examines Markov transition matrices for the non-maternal dominant and recessive cases under the *w76* four-generation property. 

We must first establish the concept of how we should track whether an individual's ancestry has "completed" the four-generation requirement to show loss. As such, let's set the four-generation requirement as $r=4$. Next, we need to analyze the cases for which each individual's ancestry will indicate loss:
- Under the dominant assumption (that *w76* is a dominant trait), an individual's ancestry needs to have (consecutively) four generations of being either heterozygous (`+/-`) or homozygous-loss (`-/-`) to be able to show loss. This immediately deviates from the case in $\textbf{3.3}$, where individuals are assumed to immediately show loss upon reaching homozygosity. 
- Under the recessive assumption, only those who have been homozygous-loss for at least four generations will show loss. Therefore, in this matrix, we need not track the generation of heterozygosity. 

Here are the matrices for $r=4$, with explanations following:
[WIP]

<!-- $$\begin{align*}M_{\text{dom,mat}^-}&=\begin{pmatrix} 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\ 0 & 0.25 & 0 & 0 & 0.5 & 0 & 0.25 & 0 \\ 0 & 0 & 0.25 & 0 & 0 & 0.5 & 0.25 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0.25 & 0.75 \\ 0 & 0 & 0 & 0 & 0 & 0 & 1 &0  \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \end{pmatrix} \\ M_{\text{rec,mat}^-}&=\begin{pmatrix} 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 \\ 0.25 & 0 & 0 & 0.5 & 0.25 & 0 \\ 0 & 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 \end{pmatrix}\end{align*}$$ -->

The rows/columns in these cases, in order, are:
$$M_{\text{dom,mat}^-}=\begin{pmatrix}\text{hom,0} \\ \text{hom,1} \\ \text{hom,2} \\ \text{hom,3} \\ \text{het,0} \\ \text{het,1} \\ \text{het,2} \\ \text{het,3} \\ \text{keep} \\ \text{loss}\end{pmatrix}, M_{\text{rec,mat}^-} = \begin{pmatrix}\text{hom,0} \\ \text{hom,1} \\ \text{hom,2} \\ \text{hom,3} \\ \text{het} \\ \text{keep} \\ \text{loss}\end{pmatrix}$$

This way, we can track (a) what stage of the four-generation requirement any individual is in, and (b) whether any individual is, at that stage, heterozygous or homozygous. 

We should adapt the input vector $\vec{\textbf{v}}^{(0)}$ as well to match the dimensions of the new $M$'s: 
$$\begin{align*}\vec{\textbf{v}}^{(0)}_{\text{dom,mat}^-} &= [0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0] \times N^{(0)}\\ \vec{\textbf{v}}^{(0)}_{\text{rec,mat}^-} &= [0.5, 0, 0, 0, 0, 0.5, 0, 0] \times N^{(0)}\end{align*}$$

The results from $M_{\text{dom,mat}^-}$ are: [WIP]


The results from $M_{\text{rec,mat}^-}$ are: [WIP]

## 3.5 Experimental Variance
The variance across trials of $n=20$ worms (referred to from now as $\sigma^2_{exp}$) was tested using a simulation in Python 3.14. It was seen that this variance was rather negligible; this means that, across experiments, there should not be much variance as to the *F6* distribution of genotypes. This intuitively makes sense: as self-ing makes the progeny trend towards homozygosity, the possible permutations for all progenital individuals decreases (associated with the quantity of heterozygous individuals). Numerically, $\sigma^2_{exp}$ was always on the magnitude of $5\times 10^{-6}$ (or $\approx 0.000005$). Note that this algorithm makes the assumption that each worm at each stage will have a number of progeny modeled by $\pi\sim N(250,35)$.


[^1]: One large concern of mine is that, if we are propagating from *F1*, some progeny should have the time to be able to keep the loss genotype for at least 4 generations before switching back to a keep genotype. A worthwhile thing to keep in mind, though I'm not sure how it affects any experiments. 

[^2]: As a note, at *F2*, the variance described can be represented as $\sigma^2 = \frac{p(1-p)}{n}$ where $p$ is the probability of being some genotype (e.g. heterozygous). For our scenario, $p=0.5$ always. If we assume that there is variance, $\sigma^2_\text{mat,rec}$ will be equal to this value. 

[^3]: The proof is left as an exercise to the reader :^)

[^4]: This is wrong; I can't compute it without establishing the Markov transition matrix first. 