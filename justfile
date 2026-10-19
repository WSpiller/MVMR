docs:
    Rscript -e "devtools::document()"
check: docs
    Rscript -e "devtools::check()"
install: docs
    Rscript -e "pkg <- pkgbuild::build(); install.packages(pkg, repos = NULL, type = 'source'); unlink(pkg)"
dev:
    Rscript -e "pak::local_install_dev_deps()"
