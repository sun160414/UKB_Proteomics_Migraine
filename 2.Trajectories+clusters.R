suppressPackageStartupMessages({
  library(MatchIt)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

## ---------------------------------------------------------
## Helpers
## ---------------------------------------------------------
bt <- function(x) ifelse(make.names(x) != x, paste0("`", x, "`"), x)

safe_z <- function(x, status){
  # Z = (case - mean(ctrl)) / sd(ctrl); return NA if sd=0 or missing
  case_val  <- x[status == 1]
  ctrl_vals <- x[status == 0]
  if (length(case_val) != 1) return(NA_real_)
  sd_ctrl <- sd(ctrl_vals, na.rm = TRUE)
  if (is.na(sd_ctrl) || sd_ctrl == 0) return(NA_real_)
  (case_val - mean(ctrl_vals, na.rm = TRUE)) / sd_ctrl
}

## ---------------------------------------------------------
## 1) Prepare matched dataset (always redefine migraine_status)
## ---------------------------------------------------------
prep_for_match <- function(df,
                           covars,
                           incident_var = "incident_migraine",
                           time_var     = "migraine_years",
                           followup_var = "followup_years"){
  
  df %>%
    mutate(
      migraine_status = as.integer(.data[[incident_var]] == 1),
      # If you already have migraine_years, keep it; otherwise use followup_years
      migraine_years  = ifelse(!is.null(df[[time_var]]),
                               as.numeric(.data[[time_var]]),
                               as.numeric(.data[[followup_var]]))
    ) %>%
    mutate(
      # Light type harmonization (edit here if your coding differs)
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
      CVD_status      = factor(CVD_status)
    ) %>%
    filter(!is.na(migraine_status), !is.na(migraine_years), migraine_years > 0) %>%
    # keep only needed columns + proteins later (no hard-coded indices)
    as.data.frame()
}

## ---------------------------------------------------------
## 2) Run matching (configurable)
## ---------------------------------------------------------
run_matching <- function(df_match,
                         covars,
                         exact_var = "sex",
                         ratio = 10,
                         method = "nearest",
                         distance = "mahalanobis"){
  
  fml <- as.formula(
    paste0("migraine_status ~ ", paste(bt(setdiff(covars, exact_var)), collapse = " + "))
  )
  
  MatchIt::matchit(
    formula  = fml,
    data     = df_match,
    method   = method,
    distance = distance,
    ratio    = ratio,
    exact    = as.formula(paste0("~ ", exact_var))
  )
}

## ---------------------------------------------------------
## 3) Add time_scale: assign each subclass the case event time (robust)
## ---------------------------------------------------------
add_time_scale <- function(mdata,
                           subclass_col = "subclass",
                           status_col   = "migraine_status",
                           years_col    = "migraine_years"){
  
  mdata %>%
    group_by(.data[[subclass_col]]) %>%
    mutate(
      case_time  = .data[[years_col]][.data[[status_col]] == 1][1],
      time_scale = -abs(case_time)  # "Years before diagnosis" (negative)
    ) %>%
    ungroup() %>%
    as.data.frame()
}

## ---------------------------------------------------------
## 4) Residualize proteins by covariates (returns eid + *_res)
## ---------------------------------------------------------
residualize_proteins <- function(mdata,
                                 covars,
                                 id_col = "eid",
                                 protein_pattern = "^Prot_",
                                 protein_cols = NULL){
  
  if (is.null(protein_cols)){
    protein_cols <- grep(protein_pattern, names(mdata), value = TRUE)
  }
  
  rhs <- paste(bt(covars), collapse = " + ")
  out <- data.frame(eid = mdata[[id_col]])
  
  for (y in protein_cols){
    keep <- c(id_col, covars, y)
    df <- mdata[, keep, drop = FALSE]
    df <- df[complete.cases(df), , drop = FALSE]
    if (nrow(df) < 20) next
    
    fml <- as.formula(paste0(bt(y), " ~ ", rhs))
    fit <- lm(fml, data = df)
    
    tmp <- data.frame(
      eid = df[[id_col]],
      res = residuals(fit)
    )
    names(tmp)[2] <- paste0(y, "_res")
    out <- dplyr::left_join(out, tmp, by = "eid")
  }
  out
}

## ---------------------------------------------------------
## 5) Case-vs-control Z per subclass (wide matrix + Year)
## ---------------------------------------------------------
subclass_z_matrix <- function(df_res,
                              subclass_col = "subclass",
                              status_col   = "migraine_status",
                              year_col     = "time_scale",
                              res_suffix   = "_res"){
  
  res_cols <- grep(paste0(res_suffix, "$"), names(df_res), value = TRUE)
  if (length(res_cols) == 0) stop("No residual columns found (e.g., *_res).")
  
  df_g <- df_res %>% group_by(.data[[subclass_col]])
  
  # Z for each residual column
  z_wide <- map_dfc(res_cols, function(rc){
    df_g %>% summarise(z = safe_z(.data[[rc]], .data[[status_col]]), .groups="drop") %>%
      select(z)
  })
  
  # Add subclass + Year (one per group)
  meta <- df_g %>%
    summarise(
      subclass = first(.data[[subclass_col]]),
      Year     = first(.data[[year_col]]),
      .groups = "drop"
    )
  
  z_wide <- bind_cols(meta, z_wide)
  
  # Clean column names: Prot_xxx_res -> Prot_xxx
  names(z_wide) <- sub(paste0(res_suffix, "$"), "", names(z_wide))
  
  z_wide
}

## ---------------------------------------------------------
## 6) Loess smoothing (long + wide export), clustering, plots
## ---------------------------------------------------------
loess_smooth_long <- function(z_wide,
                              year_col = "Year",
                              id_col   = "subclass",
                              step = 0.2,
                              span = 0.75){
  
  long <- z_wide %>%
    pivot_longer(cols = -c(all_of(id_col), all_of(year_col)),
                 names_to = "Protein", values_to = "Estimate") %>%
    filter(!is.na(Estimate))
  
  yr_min <- suppressWarnings(min(long[[year_col]], na.rm = TRUE))
  yr_max <- suppressWarnings(max(long[[year_col]], na.rm = TRUE))
  grid_year <- seq(yr_min, yr_max, by = step)
  
  sm <- map_dfr(unique(long$Protein), function(p){
    dfp <- long %>% filter(Protein == p)
    if (nrow(dfp) < 20) return(NULL)
    fit <- loess(Estimate ~ .data[[year_col]], data = dfp, span = span)
    pred <- predict(fit, newdata = data.frame(Year = grid_year))
    data.frame(Year = grid_year, Protein = p, Estimate_loess = as.numeric(pred))
  })
  
  sm
}

cluster_from_loess <- function(df_loess_sig, k = 3){
  mat <- df_loess_sig %>%
    pivot_wider(names_from = Year, values_from = Estimate_loess) %>%
    column_to_rownames("Protein") %>%
    as.matrix()
  
  mat[is.na(mat)] <- 0  # simple fallback; replace with imputation if you prefer
  
  hc <- hclust(dist(mat), method = "ward.D")
  data.frame(Protein = rownames(mat),
             cluster = factor(cutree(hc, k = k), levels = 1:k))
}

plot_heatmap <- function(df_loess_sig, cluster_df,
                         file_pdf = "Fig2A_heatmap.pdf",
                         width = 6, height = 14){
  
  heat_mat <- df_loess_sig %>%
    left_join(cluster_df, by = "Protein") %>%
    arrange(cluster, Protein) %>%
    select(-cluster) %>%
    pivot_wider(names_from = Year, values_from = Estimate_loess) %>%
    column_to_rownames("Protein") %>%
    as.matrix()
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("#2b6cb0", "white", "#c53030"))
  
  ht <- Heatmap(
    heat_mat,
    name = "Z score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 7),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 7),
    column_title = "Years before diagnosis",
    heatmap_legend_param = list(direction = "horizontal"),
    rect_gp = grid::gpar(col = "black", lwd = 0.15)
  )
  
  pdf(file_pdf, width = width, height = height, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "top")
  dev.off()
  
  invisible(ht)
}

plot_cluster_trajectories <- function(df_loess_sig, cluster_df,
                                      file_pdf = "Fig2B_trajectories.pdf",
                                      width = 6, height = 18){
  
  plot_df <- df_loess_sig %>%
    left_join(cluster_df, by = "Protein") %>%
    filter(!is.na(cluster))
  
  cluster_colors <- c("1"="#C0392B","2"="#E67E22","3"="#2980B9")
  
  plot_one <- function(cl){
    dfc <- plot_df %>% filter(cluster == cl)
    n_pro <- n_distinct(dfc$Protein)
    col <- cluster_colors[as.character(cl)]
    
    ggplot(dfc, aes(x = Year, y = Estimate_loess, group = Protein)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_line(color = col, alpha = 0.30, linewidth = 0.3) +
      stat_summary(aes(group = 1), fun = mean, geom = "line",
                   linewidth = 1.2, color = col) +
      annotate("text",
               x = min(dfc$Year, na.rm = TRUE),
               y = max(dfc$Estimate_loess, na.rm = TRUE),
               label = paste0("n = ", n_pro),
               hjust = 0, vjust = 1, size = 4) +
      labs(title = paste0("Cluster ", cl),
           x = "Years before diagnosis", y = "Z score") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank())
  }
  
  p <- plot_one("1") / plot_one("2") / plot_one("3")
  
  pdf(file_pdf, width = width, height = height, useDingbats = FALSE)
  print(p)
  dev.off()
  
  invisible(p)
}

export_cluster_lists <- function(cluster_df, prefix = "cluster_protein"){
  # Long
  long <- cluster_df %>%
    arrange(cluster, Protein) %>%
    group_by(cluster) %>%
    mutate(rank_in_cluster = row_number()) %>%
    ungroup()
  fwrite(long, paste0(prefix, "_list_long.csv"))
  
  # Collapsed
  collapsed <- cluster_df %>%
    arrange(cluster, Protein) %>%
    group_by(cluster) %>%
    summarise(n_proteins = n(),
              proteins = paste(Protein, collapse = ", "),
              .groups = "drop")
  fwrite(collapsed, paste0(prefix, "_collapsed.csv"))
  
  # Wide
  wide <- cluster_df %>%
    arrange(cluster, Protein) %>%
    group_by(cluster) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    select(cluster, idx, Protein) %>%
    pivot_wider(names_from = cluster, values_from = Protein) %>%
    arrange(idx)
  fwrite(wide, paste0(prefix, "_list_wide.csv"))
}

## =========================================================
## ===================== PIPELINE RUN =======================
## =========================================================

## ---- User inputs ----
covars <- c(
  "age","sex","ethn","Qualification","bmi","Socioeconomic",
  "Smoking_status","Alcohol_consumption",
  "screen time (TV)","screen time (computer)","sleep duration",
  "diabetes_status","CVD_status"
)

## 0) Prepare data for matching
df_match <- prep_for_match(mydata_base, covars = covars,
                           incident_var = "incident_migraine",
                           followup_var = "followup_years")

## 1) Matching
m.out <- run_matching(df_match, covars = covars, exact_var = "sex", ratio = 10)
m.data <- match.data(m.out)

## Quick checks (optional)
print(table(m.data$migraine_status))
print(table(table(m.data$subclass)))

## 2) Add time_scale (robust, no rep(each=11))
m.data <- add_time_scale(m.data)

fwrite(m.data, "match_incident_migraine.csv")

## 3) Residualize proteins (auto-detect Prot_)
res_mat <- residualize_proteins(m.data, covars = covars,
                                id_col = "eid", protein_pattern = "^Prot_")

## Merge back key meta columns needed downstream (minimal)
meta_keep <- c("eid","subclass","migraine_status","time_scale")
df_res <- left_join(res_mat, m.data[, meta_keep, drop = FALSE], by = "eid")

fwrite(df_res, "resid_ukb_protein.csv")

## 4) Subclass-level Z matrix (case vs controls)
z_wide <- subclass_z_matrix(df_res)
# z_wide: subclass + Year + Protein columns

## 5) LOESS smoothing (long)
df_loess <- loess_smooth_long(z_wide, step = 0.2, span = 0.75)

## 6) Protein name mapping (optional; keep both raw + mapped)
df_loess2 <- df_loess %>%
  mutate(Protein_raw = str_trim(Protein),
         Protein_name = unname(reverse_map[Protein_raw]))

fwrite(df_loess,  "df_loess_raw.csv")
fwrite(df_loess2, "df_loess_mapped.csv")

## 7) Filter to significant proteins (sig must have column 'Protein')
df_loess_sig <- df_loess %>% semi_join(sig, by = "Protein")

## 8) Clustering (k=3) + plots + exports
cluster_df <- cluster_from_loess(df_loess_sig, k = 3)
fwrite(cluster_df, "cluster_df.csv")

plot_heatmap(df_loess_sig, cluster_df, file_pdf = "Fig2A_heatmap.pdf", width = 6, height = 14)
plot_cluster_trajectories(df_loess_sig, cluster_df, file_pdf = "Fig2B_trajectories.pdf", width = 6, height = 18)

export_cluster_lists(cluster_df, prefix = "cluster_protein")

save.image("Data3.Rdata")