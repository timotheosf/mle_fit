# mle_fit 📈

**High-Performance Modern Fortran Library for Empirical Power-Law Fitting**

This project implements robust and optimized Maximum Likelihood Estimation (MLE) methods for fitting heavy-tailed and power-law distributions to empirical data, originally based on the theoretical framework by [Clauset, Shalizi, and Newman (2009)](http://arxiv.org/abs/0706.1062). 

The project was initially focused on Python-Fortran integration, developed as part of the thesis [*Complexidade e Emergência em Linguagens Naturais*](https://drive.google.com/file/d/1VF6BgdaVsb-DG5UdC0VnKhsgc-dNbC1J/view?usp=sharing) at the Federal University of Viçosa (UFV), Brazil. Now, it has been completely rewritten in **Modern Fortran (OOP)** to leverage SIMD vectorization, multi-threading, and advanced algorithmic optimizations.

## ✨ Key Features

* **Extreme Performance:** Implements an $O(1)$ dynamic cumulative $\alpha$ calculation and a vectorized Kolmogorov-Smirnov (KS) statistic evaluation, fitting datasets with hundreds of thousands of elements in a fraction of a second.
* **OpenMP Parallelization:** The computationally expensive Monte Carlo goodness-of-fit test (p-value) is fully parallelized utilizing OpenMP.
* **OOP API:** Clean, intuitive Object-Oriented design using the `empirical_pl` class.
* **Unified Fitting Interface:** A single, powerful `fast_fit` method that handles both standard Clauset et al. approaches and advanced weighted adjustments.
* **Tail-Weighted MLE:** A regularized KS-statistic approach that prevents the algorithm from overestimating $x_{\min}$ in noisy empirical datasets (see theoretical details below).
* **FPM Ready:** Fully integrated with the Fortran Package Manager (`fpm`).

## 📦 Installation

To use `mle_fit` in your Fortran project, simply add it as a dependency in your `fpm.toml` file:

```toml
[dependencies]
mle_fit = { git = "https://github.com/timotheosf/mle_fit.git" }
```

> **Note:** To leverage the OpenMP parallel processing for p-value generation, ensure your compiler supports it and add the appropriate flag (e.g., `-fopenmp` for GCC) in your build command or `fpm.toml` compiler options.

Then, build your project using the release profile for maximum optimization:
```bash
fpm build --profile release --compiler gfortran --flag "-O3 -fopenmp"
```

## ⚙️ Core Functionalities

`mle_fit` utilizes Fortran classes to provide a clean and user-friendly API. The main class is `empirical_pl`, which handles data processing and fitting seamlessly. Current capabilities include:

* **Polymorphic Data Input:** Ingest raw arrays of almost any numerical type (`integer(i1..i8)`, `real(sp)`, `real(dp)`). It automatically detects if the data is discrete or continuous, converting it to double-precision internally for optimized handling.
* `init`: Receives and automatically sorts the empirical data using state-of-the-art radix/ordinal sort.
* `fast_fit`: Executes the robust search algorithm. It supports optional arguments like `lambda_in` (for penalization) and `use_weight` (for Anderson-Darling style KS weighting).
* `p_value`: Generates synthetic synthetic datasets in parallel to compute the statistical goodness-of-fit.
* `report`: Instantly prints a clean summary of the estimated parameters.

## 🚀 Quick Start

Here is a basic example of how to fit a power-law to an array of empirical data:

```fortran
program main
    use kinds_mod, only: dp
    use power_law_mle_fit_mod, only: empirical_pl
    implicit none

    type(empirical_pl) :: pl
    real(dp), allocatable :: data(:)
    real(dp) :: p_val
    
    ! 1. Load or generate your empirical data (can be any numerical type)
    ! allocate(data(100000))
    ! data = ... 
    
    ! 2. Initialize the empirical_pl object (data is internally sorted)
    call pl%init( data )
    
    ! 3. Standard Clauset et al. Fit (pure KS-distance)
    call pl%fast_fit( use_weight=.false. )
    
    ! 4. Regularized Fit (Recommended for noisy data - uses AD weight and lambda penalty)
    call pl%fast_fit( use_weight=.true., lambda_in=0.05_dp )
    
    ! 5. Calculate Goodness-of-Fit (P-value) using OpenMP (e.g., 1000 iterations)
    p_val = pl%p_value( N_samples=1000 )
    
    ! 6. Print the results
    call pl%report()

end program main
```

## 🔬 Theoretical Highlights: The Tail-Weighted MLE

The standard Clauset et al. method minimizes the KS distance between the empirical CDF and the theoretical power-law CDF. However, in highly noisy empirical datasets, the standard KS statistic often overestimates the $x_{\min}$, because smaller tails are statistically easier to fit.

To solve this, we introduce a **Loss Function** defined as:
$$ \mathcal{L}[F, x_{\min}] = \mathcal{K}_w[F] - \lambda \left( \frac{n_{\text{tail}}(x_{\min})}{N} \right)^2,$$
where:
* $\mathcal{K}_w[F]$ is the KS statistic (optionally weighted by the Anderson-Darling variance approach to emphasize tail behavior).
* $n_{\text{tail}} / N$ is the fraction of data kept in the power-law tail.
* $\lambda$ is the regularization hyperparameter (default $\approx 0.05$).

The quadratic reward heavily penalizes the algorithm for discarding large portions of the dataset due to small local fluctuations, balancing the statistical measure, $\mathcal{K}[F]$, and the tail length, $n_{\text{tail}}$. As the tail gets smaller, the quadratic term decays to zero, naturally restoring the standard behavior: $\mathcal{L}[F, x_{\min}]\to\mathcal{K}_w[F]$.

## 🗺️ Roadmap & Future Work

While the Fortran numerical core is complete and highly optimized, the following features are planned for upcoming releases:
- [ ] **Advanced Statistics:** Implementation of the Log-Likelihood ratio test and Vuong's closeness test to compare fits.
- [ ] **Other Distributions:** Addition of exponential, log-normal, normal distributions, and power-law with exponential cut-off.
- [ ] **Python Bindings:** Re-introducing Python wrappers using `f2py` or `ctypes` to bring Fortran performance to Python data pipelines.
- [ ] **Plotting Utilities:** Helper scripts for standard Log-Log PDF/CCDF visualizations.

## 📜 License & Citation

*Note: This project is in active development.*

