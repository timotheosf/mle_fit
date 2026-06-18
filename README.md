# mle_fit 📈

**High-Performance Modern Fortran Library for Empirical Distributions Fitting**

This project implements robust and optimized Maximum Likelihood Estimation (MLE) methods for fitting heavy-tailed and probability distributions to empirical data, originally based on the theoretical framework for power-laws by [Clauset, Shalizi, and Newman (2009)](http://arxiv.org/abs/0706.1062). 

Originally focused on Python-Fortran integration as part of the thesis [*Complexidade e Emergência em Linguagens Naturais*](https://drive.google.com/file/d/1VF6BgdaVsb-DG5UdC0VnKhsgc-dNbC1J/view?usp=sharing) at the Federal University of Viçosa (UFV), Brazil. Now, it has been completely rewritten in **Modern Fortran (OOP)** to leverage SIMD vectorization, multi-threading, and advanced algorithmic optimizations. 

*"It's bigger! It's badder! Ladies and gentlemen, it's too much for Mr. Incredible!"*

## ✨ Key Features

* **Polymorphic Architecture:** Built on a Strategy Pattern using Fortran Abstract Interfaces. The central MLE engine is completely agnostic, treating any probability distribution (Power Law, Lognormal, Exponential) as a plug-and-play module.
* **Extreme Performance:** Implements an $O(1)$ dynamic cumulative parameter calculation and a vectorized Kolmogorov-Smirnov (KS) statistic evaluation, fitting datasets with hundreds of thousands of elements in a fraction of a second.
* **OpenMP Parallelization:** The computationally expensive Monte Carlo goodness-of-fit test (p-value) is fully parallelized utilizing OpenMP with thread-safe RNGs.
* **Data Science Ready:** Natively exports fit results and parameters to Pandas-compatible `.csv` files or direct Numpy/Python variable files (`.py`), eliminating the need for complex C/Fortran bindings.
* **Tail-Weighted & Fixed MLE:** Includes regularized KS-statistic approaches for noisy data, and allows fitting distributions from an arbitrary, user-defined $x_{\min}$ cut-off.

## 🏗️ Architecture & Logical Tree

The library is designed with strict separation of concerns. The mathematical definitions of the distributions are isolated from the heavy-lifting optimization loops:

```text
src/
 ├── mle_kinds.f90        (Foundation: precision types, dynamic data arrays, timers)
 │
 ├── dist_interface.f90   (The Mold: 'empirical_distribution' abstract interface)
 │
 ├── pl_mod.f90           (Distribution: extends abstract class with Power-Law O(1) math)
 ├── exp_mod.f90          (Distribution: extends abstract class with Exponential math...)
 │
 └── mle_fit.f90          (The Engine: handles data, KS minimization, and OpenMP P-values)
```

## 📦 Installation

To use `mle_fit` in your Fortran project, simply add it as a dependency in your `fpm.toml` file:

```toml
[dependencies]
mle_fit = { git = "https://github.com/timotheosf/mle_fit.git" }
```

> **Note:** To leverage the OpenMP parallel processing for p-value generation, ensure your compiler supports it and add the appropriate flag (e.g., `-fopenmp` for GCC) in your build command or `fpm.toml` compiler options.

Build your project using the release profile for maximum optimization:
```bash
fpm build --profile release --compiler gfortran --flag "-O3 -fopenmp"
```

## 🚀 Quick Start

Here is how you use the new Object-Oriented API to fit a Power-Law to your data and export it to Python:

```fortran
program main
    use mle_kinds_mod
    use pl_mod          ! Import the Distribution
    use mle_fit_mod     ! Import the Engine
    implicit none

    type(mle) :: engine
    type(power_law) :: my_pl
    integer(i4) :: file_unit
    
    ! 1. Initialize the specific distribution
    call my_pl%start_pl()
    
    ! 2. Load your empirical data (auto-detects integer/real and sorts it)
    open(newunit=file_unit, file='my_data.dat', action='read', status='old')
    call engine%data%from_file( file_unit )
    close(file_unit)
    
    ! 3. Fit the data by dynamically binding the distribution to the engine
    call engine%fast_fit( distribuiton=my_pl )
    
    ! Optional: Force the fit starting from an arbitrary x_min cut-off
    ! call engine%fixed_fit( distribuiton=my_pl, fixed_xmin=100.0_dp )
    
    ! 4. Calculate Goodness-of-Fit (P-value) using OpenMP (e.g., 2500 iterations)
    call engine%p_value( N_samples=2500 )
    
    ! 5. Print a beautiful report to the terminal
    call engine%report()

    ! 6. Export results seamlessly to Python!
    call engine%write_pandas("results_dataframe.csv")  ! Load with pd.read_csv()
    call engine%write_numpy("results_vars.py")         ! Load with 'import results_vars'

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

