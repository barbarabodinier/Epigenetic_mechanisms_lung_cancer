untar("sharp_1.0.1.tar.gz")
devtools::install("sharp", upgrade = "always")
devtools::install_version("glmnet", version = "4.1") # required to reproduce results
