
## Creamos un repositorio para el proyecto
usethis::create_project("~/DE_project")

## Generamos un script de R con las especificaciones del setup
usethis::use_r("01-Setup")

## Directorios
dir_plots <- here::here("figuras")
dir_rdata <- here::here("processed-data")

## Creamos los directorios para las figuras y archivos
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

## Conectamos con GitHub
usethis::create_github_token()

## Inicializamos el repositorio de Git
usethis::use_git()

## Conectamos el repositorio local de Git con los servidores de GitHub
usethis::use_github()
