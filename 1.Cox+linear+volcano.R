## =========================================================
## migraine proteomics pipeline
## - Cox + LM + Volcano (FDR threshold)
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(data.table)
  library(survival)
  library(ggplot2)
  library(ggrepel)
})

fdr_cutoff <- 0.05

cox_covars <- c(
  "age","sex","ethn","Qualification","bmi","Socioeconomic",
  "Smoking_status","Alcohol_consumption",
  "screen time (TV)","screen time (computer)","sleep duration",
  "diabetes_status","CVD_status"
)

## Add backticks for non-syntactic names
bt <- function(x) ifelse(make.names(x) != x, paste0("`", x, "`"), x)

## ---- 1) Covariate harmonization ----
define_covariates <- function(df){
  df %>%
    mutate(
      age = as.numeric(age),
      bmi = as.numeric(bmi),
      Socioeconomic = as.numeric(Socioeconomic),
      Qualification = as.numeric(Qualification),
      
      sex  = factor(sex),
      ethn = factor(ethn),
      Smoking_status      = factor(Smoking_status),
      Alcohol_consumption = factor(Alcohol_consumption),
      
      `screen time (TV)`       = as.numeric(`screen time (TV)`),
      `screen time (computer)` = as.numeric(`screen time (computer)`),
      `sleep duration`         = as.numeric(`sleep duration`),
      
      diabetes_status = factor(diabetes_status),
      CVD_status      = factor(CVD_status),
      
      baseline_date = as.Date(baseline_date)
    ) %>%
    mutate(across(where(is.character), ~na_if(., "")))
}

## ---- 2) Data preparation (re-define migraine_status each time) ----
prep_data <- function(df,
                      mode = c("cox_incident", "lm_prevalent"),
                      time_var = "followup_years",
                      incident_var = "incident_migraine",
                      prevalent_var = "prevalent_migraine"){
  mode <- match.arg(mode)
  df <- define_covariates(df)
  
  if (mode == "cox_incident"){
    df %>%
      filter(.data[[prevalent_var]] == 0) %>%                 # baseline-free
      mutate(
        migraine_status = as.integer(.data[[incident_var]] == 1),
        migraine_years  = as.numeric(.data[[time_var]])
      ) %>%
      filter(!is.na(migraine_years) & migraine_years > 0,
             !is.na(migraine_status))
  } else {
    df %>%
      mutate(
        migraine_status = as.integer(.data[[prevalent_var]] == 1) # cross-sectional
      ) %>%
      filter(!is.na(migraine_status))
  }
}

## ---- 3) Batch Cox ----
run_cox_batch <- function(dat0, protein_cols, covars, min_n = 50, min_event = 5){
  cov_terms <- bt(covars)
  
  map_dfr(protein_cols, function(p){
    df <- dat0 %>%
      transmute(
        migraine_years, migraine_status,
        Protein = as.numeric(.data[[p]]),
        across(all_of(covars))
      ) %>%
      mutate(Protein_z = as.numeric(scale(Protein))) %>%
      filter(complete.cases(.))
    
    if (nrow(df) < min_n || sum(df$migraine_status == 1) < min_event){
      return(data.frame(Protein = p, n = nrow(df), nevent = sum(df$migraine_status == 1),
                        HR = NA_real_, p = NA_real_))
    }
    
    form_txt <- paste0(
      "Surv(migraine_years, migraine_status) ~ Protein_z + ",
      paste(cov_terms, collapse = " + ")
    )
    
    fit <- coxph(as.formula(form_txt), data = df)
    s <- summary(fit)
    
    data.frame(
      Protein = p,
      n = s$n, nevent = s$nevent,
      HR = s$coef["Protein_z","exp(coef)"],
      p  = s$coef["Protein_z","Pr(>|z|)"]
    )
  })
}

## ---- 4) Batch linear regression ----
run_lm_batch <- function(dat0, protein_cols, covars, min_n = 50){
  cov_terms <- bt(covars)
  
  map_dfr(protein_cols, function(p){
    df <- dat0 %>%
      transmute(
        Protein = as.numeric(.data[[p]]),
        migraine_status,
        across(all_of(covars))
      ) %>%
      filter(complete.cases(.))
    
    if (nrow(df) < min_n){
      return(data.frame(Protein = p, n = nrow(df),
                        beta = NA_real_, se = NA_real_, pvalue = NA_real_))
    }
    
    form_txt <- paste0(
      "Protein ~ migraine_status + ",
      paste(cov_terms, collapse = " + ")
    )
    
    fit <- lm(as.formula(form_txt), data = df)
    s <- summary(fit)$coefficients
    
    data.frame(
      Protein = p,
      n = nrow(df),
      beta   = s["migraine_status","Estimate"],
      se     = s["migraine_status","Std. Error"],
      pvalue = s["migraine_status","Pr(>|t|)"]
    )
  })
}

## ---- 5) Postprocess: FDR + mapping + export ----
postprocess_map_export <- function(df, p_col = c("p","pvalue"), reverse_map, out_prefix){
  p_col <- match.arg(p_col)
  
  df2 <- df %>%
    mutate(
      p_fdr = p.adjust(.data[[p_col]], method = "BH"),
      Protein_raw  = trimws(Protein),
      Protein_name = unname(reverse_map[Protein_raw])
    )
  
  fwrite(df2, paste0(out_prefix, "_FDR_mapped.csv"))
  return(df2)
}

## =========================================================
## ====================== RUN (Examples) ====================
## =========================================================

## Auto-detect protein columns
protein_cols <- grep("^Prot_", names(dat), value = TRUE)

## A) Cox overall (incident)
dat_cox <- prep_data(dat, mode="cox_incident",
                     time_var="followup_years",
                     incident_var="incident_migraine",
                     prevalent_var="prevalent_migraine")

cox_df <- run_cox_batch(dat_cox, protein_cols, cox_covars)

cox_out <- postprocess_map_export(cox_df, p_col="p", reverse_map=reverse_map, out_prefix="cox_overall")

## B) LM overall (prevalent)
dat_lm <- prep_data(dat, mode="lm_prevalent", prevalent_var="prevalent_migraine")
lm_df  <- run_lm_batch(dat_lm, protein_cols, cox_covars)

lm_out <- postprocess_map_export(lm_df, p_col="pvalue", reverse_map=reverse_map, out_prefix="lm_overall")
