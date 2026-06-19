build_and_deploy <- function() {
  
  deploy_dir <- "../ShinyScab"
  dir.create(deploy_dir, showWarnings = FALSE)
  
  file.copy("app.R", file.path(deploy_dir, "app.R"), overwrite = TRUE)
  file.copy(c("run_scenario.R", "plot_scenario.R"), deploy_dir, overwrite = TRUE)
  file.copy("DESCRIPTION", file.path(deploy_dir, "DESCRIPTION"), overwrite = TRUE)
  
  rsconnect::deployApp(appDir = deploy_dir)
}

unlink("../ShinyScab", recursive = TRUE)
build_and_deploy()