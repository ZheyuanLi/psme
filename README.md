# Penalized Splines Mixed-Effects Models (psme)


You can install this **R** package using the following **R** code:

```
## you may need to first install package "remotes" from CRAN
remotes::install_github("ZheyuanLi/psme")
```

# Caveat

This package has a dependency chain like: **mgcv** -> **psme** <- **lme4** <- **Matrix**.

It is important to keep your installed versions of **lme4** and **Matrix** up to date, otherwise **psme** may not work. For example, if you use a **Matrix** version released in 2025 but an **lme4** version released in 2024, you will be in trouble.
