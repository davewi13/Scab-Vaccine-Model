build_and_deploy <- function() {
  
  deploy_dir <- "../ShinyScab"
  
  dir.create(deploy_dir, showWarnings = FALSE)
  
  file.copy("app.R", file.path(deploy_dir, "app.R"), overwrite = TRUE)
  file.copy(c("run_scenario.R", "plot_scenario.R"),
            deploy_dir,
            overwrite = TRUE)
  
  
  file.copy(
    from = "../SimInf",
    to = deploy_dir,
    recursive = TRUE,
    overwrite = TRUE
  )
  
  
  setwd(deploy_dir)
  rsconnect::deployApp()
}

build_and_deploy()
