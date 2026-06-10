# mle_fit 📈

**High-Performance Modern Fortran Library for Empirical Power-Law Fitting**

This project implements robust and optimized Maximum Likelihood Estimation (MLE) methods for fitting heavy-tailed and power-law distributions to empirical data, originally based on the theoretical framework by [Clauset, Shalizi, and Newman (2009)](http://arxiv.org/abs/0706.1062). 

The project was initially focused on Python-Fortran integration, developed as part of the thesis [*Complexidade e Emergência em Linguagens Naturais*](https://drive.google.com/file/d/1VF6BgdaVsb-DG5UdC0VnKhsgc-dNbC1J/view?usp=sharing) at the Federal University of Viçosa (UFV), Brazil. Now, it has been completely rewritten in **Modern Fortran (OOP)** to leverage SIMD vectorization and advanced algorithmic optimizations.

## ✨ Key Features

* **Extreme Performance:** Implements an $O(1)$ dynamic cumulative $\alpha$ calculation and an AVX/SIMD-vectorized Kolmogorov-Smirnov (KS) statistic evaluation, fitting datasets with hundreds of thousands of elements in a fraction of a second.
* **OOP API:** Clean, intuitive Object-Oriented design using the `empirical_pl` class.
* **Standard MLE (Clauset et al.):** Classic KS-distance minimization to find the optimal lower bound $x_{\min}$ and exponent $\alpha$.
* **Tail-Weighted MLE:** A regularized KS-statistic approach that prevents the algorithm from overestimating $x_{\min}$ in noisy empirical datasets (see theoretical details below).
* **FPM Ready:** Fully integrated with the Fortran Package Manager (`fpm`).

## 📦 Installation

To use `mle_fit` in your Fortran project, simply add it as a dependency in your `fpm.toml` file:

```toml
[dependencies]
mle_fit = { git = "https://github.com/timotheosf/mle_fit.git" }
```
Then, build your project using the release profile for maximum optimization:
```bash
fpm build --profile release
```

## ⚙️ Core Functionalities

`mle_fit` utilizes Fortran classes to provide a clean and user-friendly API. The main class is `empirical_pl`, which handles data processing and fitting seamlessly. Current capabilities include:

* **Polymorphic Data Input:** The class can ingest raw arrays of almost any numerical type (e.g., `integer(i4)`, `integer(i8)`, `real(sp)`, `real(dp)`). It automatically detects if the data is discrete or continuous, converting it to double-precision internally for optimized handling.
* `fast_fit`: Executes the standard Clauset et al. search algorithm to find the best $x_{\min}$ and $\alpha$ using pure KS-distance minimization.
* `wfast_fit`: Executes the *Weighted* fast fit, applying the decay loss function (regularization) to prevent tail shrinkage in the presence of statistical noise.

## 🚀 Quick Start

Here is a basic example of how to fit a power-law to an array of empirical data:

```fortran
    ! 0. Import the modules
    use kinds_mod, only: dp
    use power_law_mle_fit_mod, only: empirical_pl
    implicit none

    type(empirical_pl) :: pl
    real(dp) :: x_min, alpha
    
    ! 1. Load or generate your empirical data
    ! data = ... -> data array can be of any numerical type
    
    ! 2. Initialize the empirical_pl object (automatically sorts the data)
    call pl%init( data )
    
    ! 3. Standard Clauset et al. Fit
    call pl%fast_fit( xmin=x_min, alpha=alpha )
    print *, "Standard MLE -> x_min:", x_min, " | alpha:", alpha
    
    ! 4. Regularized Fit (Recommended for noisy empirical data)
    call pl%wfast_fit( xmin=x_min, alpha=alpha )
    print *, "Weighted MLE -> x_min:", x_min, " | alpha:", alpha
```

## 🔬 Theoretical Highlights: The Tail-Weighted MLE

The standard Clauset et al. method minimizes the KS distance between the empirical CDF and the theoretical power-law CDF. However, in highly noisy empirical datasets, the standard KS statistic often overestimates the $x_{\min}$, because smaller tails are statistically easier to fit.

To solve this, we introduce a **Loss Function** defined as:
$$ \mathcal{L}[F, x_{\min}] = \mathcal{K}[F] - \lambda \left( \frac{n_{\text{tail}}(x_{\min})}{N} \right)^2,$$
where:
* $\mathcal{K}[F]$ is the standard KS statistic.
* $n_{\text{tail}} / N$ is the fraction of data kept in the power-law tail.
* $\lambda$ is the regularization hyperparameter (default $\approx 0.05$).

The quadratic reward heavily penalizes the algorithm for discarding large portions of the dataset due to small local fluctuations, balancing the statistical measure, $\mathcal{K}[F]$, and the tail length, $n_{\text{tail}}$. As the tail gets smaller, the quadratic term decays to zero, naturally restoring the standard behavior: $\mathcal{L}[F, x_{\min}]\to\mathcal{K}[F]$.

## 🗺️ Roadmap & Future Work

While the Fortran numerical core is complete and highly optimized, the following features are planned for upcoming releases:
- [ ] **Goodness-of-fit:** Re-implementation of the Monte Carlo p-value, log-likelihood, and Vuong's closeness test.
- [ ] **Other distributions:** Implementation of exponential, log-normal, normal distributions, and power-law with exponential cut-off.
- [ ] **Python Bindings:** Re-introducing Python wrappers (via `f2py` or `ctypes`).
- [ ] **Plotting Utilities:** Python helper scripts for standard Log-Log PDF/CCDF visualizations using `matplotlib`.

## 📜 License & Citation

*Note: This project is in active development.*

*(Add your license type here, e.g., MIT License)*