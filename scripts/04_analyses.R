library(readr)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(vcfR)
library(adegenet)
library(maps)
library(sf)
library(phangorn)
library(StAMPP)
library(ggtree)
library(ggrepel)
library(glue)
library(patchwork)
library(clusterProfiler)
library(org.Gg.eg.db)
library(AnnotationDbi)

# read in samples
samples <- read_tsv("/home/k14m234/erythrura/scripts/data/samples.tsv")
pts <- st_as_sf(samples, coords = c("long", "lat"), crs = 4326)

# plot by species and country
map <- map_data("world")
trich_range <- st_read("/home/k14m234/erythrura/scripts/data/shp/Erythrura_trichroa_22719712.shp", quiet = TRUE) |> st_make_valid()
pap_range <- st_read("/home/k14m234/erythrura/scripts/data/shp/Erythrura_papuana_22719721.shp", quiet = TRUE) |> st_make_valid()
ranges <- bind_rows(trich_range, pap_range)
world <- st_as_sf(maps:::map("world", plot = FALSE, fill = TRUE))
world <- st_make_valid(world) |> st_transform(4326)
bb <- st_bbox(pts)
pad <- 2  # degrees
xlim <- c(bb["xmin"] - pad, bb["xmax"] + pad)
ylim <- c(bb["ymin"] - pad, bb["ymax"] + pad)

ranges2 <- ranges %>%
  mutate(species = case_when(
    str_detect(tolower(SCINAME), "papuana")  ~ "papuana",
    str_detect(tolower(SCINAME), "trichroa") ~ "trichroa",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(species))

# ensure pts species matches too
pts2 <- pts %>%
  mutate(species = tolower(species))

xlim <- unname(as.numeric(xlim))
ylim <- unname(as.numeric(ylim))

# bounding box to deal with polygon issue
bbox <- sf::st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                      ymin = ylim[1], ymax = ylim[2]),
                    crs = sf::st_crs(4326))

# crop world basemap
world_crop <- world %>%
  sf::st_make_valid() %>%
  sf::st_crop(bbox)

# crop shpfiles
ranges_crop <- ranges2 %>%
  sf::st_make_valid() %>%
  sf::st_crop(bbox)

# fix palette
pal_named <- c(papuana = "yellowgreen",
               trichroa = "dodgerblue1")

# jittering
pts2   <- sf::st_set_crs(pts2, 4326)
pts_p  <- sf::st_transform(pts2, 3857)
pts_jit <- pts_p %>%
  mutate(geometry = sf::st_jitter(geometry, amount = 20000))

# plot
p1 <- ggplot() +
  geom_sf(data = world_crop, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = ranges_crop, aes(fill = species), color = NA, alpha = 0.25) +
  geom_sf(data = pts_jit, aes(fill = species), shape = 21, color = "black",
          size = 2.5, alpha = 0.9, stroke = 0.2) +
  scale_fill_manual(values = pal_named, name = "Species") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "grey92", linewidth = 0.2))

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/map.pdf",
  plot = p1,
  device = cairo_pdf,
  width = 6,
  height = 4,
  units = "in"
)

# initial pca
vcf <- read.vcfR("/home/k14m234/erythrura_assembly/results/GCF_005870125.1/QC/erythrura.pruned.vcf.gz")
dna <- vcfR2DNAbin(vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
erythrura <- DNAbin2genind(dna)
erythrura_scaled <- scaleGen(erythrura,NA.method="mean",scale=F)
pca <- prcomp(erythrura_scaled, center=F,scale=F)
write.csv(pca$x, "/home/k14m234/erythrura/scripts/data/pca.csv")
screeplot(pca)
pcs <- data.frame(pca$x[,1:3])
pcs$prep_ID <- rownames(erythrura@tab)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100

# % variance explained
var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100

# scores + metadata
df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(samples, by = "sample_id")

# lock factor order so colors don’t flip
df$species <- factor(df$species, levels = names(pal_named))

# write pcs
write_tsv(df, "/home/k14m234/erythrura/scripts/data/pca_sub.tsv")

# use same axis limits for both plots
xlim <- range(df$PC1, na.rm = TRUE)
ylim <- range(df$PC2, na.rm = TRUE)

# set shp values
shape_vals <- c(21, 22, 23, 24, 25)
names(shape_vals) <- unique(df$country)

# plot with all samples
p2 <- ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(fill = species, shape = country),
             size = 2.5, alpha = 0.9, color = "black", stroke = 0.2) +
  scale_fill_manual(values = pal_named, name = "Species", drop = FALSE) +
  scale_shape_manual(values = shape_vals, name = "Country") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, colour = "black")),
    shape = guide_legend(override.aes = list(fill = "white"))
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]))

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/pca_all.pdf",
  plot = p2,
  device = cairo_pdf,
  width = 5,
  height = 4,
  units = "in"
)

# subset to drop far away countries
drop_countries <- c("australia", "palau")
ids <- adegenet::indNames(erythrura)
samples2 <- samples %>%
  mutate(country_lc = tolower(country))
keep_ids <- samples2 %>%
  filter(!country_lc %in% drop_countries) %>%
  pull(sample_id)
keep_ids <- intersect(keep_ids, ids)
erythrura_sub <- erythrura[keep_ids, ]

# subset genind BEFORE scaling + PCA
erythrura_sub <- erythrura[keep_ids, ]

# scale + PCA
erythrura_sub_scaled <- scaleGen(erythrura_sub, NA.method = "mean", scale = FALSE)
pca_sub <- prcomp(erythrura_sub_scaled, center = FALSE, scale. = FALSE)

# % variance explained for THIS PCA
var_exp_sub <- (pca_sub$sdev^2) / sum(pca_sub$sdev^2) * 100

# export as df
df_sub <- as.data.frame(pca_sub$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(samples, by = "sample_id")

# write scores if you want
write_tsv(df_sub, "/home/k14m234/erythrura/scripts/data/pca_sub.tsv")

# optional scree
screeplot(pca)

p3 <- ggplot(df_sub, aes(PC1, PC2)) +
  geom_point(aes(fill = species, shape = country),
             size = 2.5, alpha = 0.9, color = "black", stroke = 0.2) +
  scale_fill_manual(values = pal_named, name = "Species", drop = FALSE) +
  scale_shape_manual(values = shape_vals, name = "Country") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, colour = "black")),
    shape = guide_legend(override.aes = list(fill = "white"))
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp_sub[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp_sub[2]))

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/pca_sub.pdf",
  plot = p3,
  device = cairo_pdf,
  width = 5,
  height = 4,
  units = "in"
)

# admixture plots!
fam_path  <- "/home/k14m234/erythrura/results/admixture/erythrura.admix.fam"
runs_root <- "/home/k14m234/erythrura/results/admixture/runs"
read_fam <- function(path) {
  read_table(path, col_names = c("FID","IID","PID","MID","SEX","PHENO"),
              col_types = cols(.default = col_character()))
}

parse_species <- function(iid) {
  # based on your IDs: ..._papuana or ..._trichroa
  case_when(
    str_detect(iid, "_papuana$")  ~ "papuana",
    str_detect(iid, "_trichroa$") ~ "trichroa",
    TRUE ~ NA_character_
  )
}
parse_cv <- function(log_text) {
  # e.g. "CV error (K=2): 0.72541"
  m <- str_match(log_text, "CV error \\(K=\\d+\\):\\s*([0-9.]+)")
  as.numeric(m[,2])
}

read_q <- function(q_path) {
  # read without colnames, let readr determine number of cols
  q <- readr::read_table2(q_path, col_names = FALSE, col_types = readr::cols())
  # rename to Q1..QK based on actual K
  names(q) <- paste0("Q", seq_len(ncol(q)))
  q
}

fam <- read_fam(fam_path) %>%
  mutate(species = parse_species(IID))

q_files <- list.files(runs_root, pattern = "\\.Q$", recursive = TRUE, full.names = TRUE)

runs <- tibble(q_path = q_files) %>%
  mutate(
    # parse K and rep from path
    K   = as.integer(str_match(q_path, "/K(\\d+)/")[,2]),
    rep = as.integer(str_match(q_path, "/rep(\\d+)/")[,2]),
    dir = dirname(q_path),
    log_path = file.path(dir, "run.log")
  ) %>%
  filter(!is.na(K), !is.na(rep)) %>%
  distinct(K, rep, .keep_all = TRUE) %>%
  arrange(K, rep)

runs <- runs %>% 
  mutate(log_txt = map_chr(log_path, ~ if (file.exists(.x)) read_file(.x) else ""), 
    cv = map_dbl(log_txt, parse_cv)) %>% select(-log_txt)

best_runs <- runs %>%
  group_by(K) %>%
  arrange(cv, rep, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

check <- best_runs %>%
  mutate(ncol_q = purrr::map_int(q_path, ~ ncol(readr::read_table2(.x, col_names = FALSE)))) %>%
  mutate(ok = (ncol_q == K))
print(check)

write_tsv(runs, file.path("/home/k14m234/erythrura/results/admixture", "cv_all_reps.tsv"))
write_tsv(best_runs, file.path("/home/k14m234/erythrura/results/admixture", "cv_best_by_K.tsv"))

best_q_long <- best_runs %>%
  mutate(q_df = purrr::map(q_path, read_q)) %>%
  mutate(
    q_long = purrr::map2(q_df, K, ~{
      # .x is wide Q table with exactly K columns for this run
      stopifnot(ncol(.x) == .y)
      names(.x) <- paste0("Q", seq_len(.y))
      dplyr::bind_cols(fam %>% dplyr::select(IID, species), .x) %>%
        tidyr::pivot_longer(dplyr::starts_with("Q"),
                            names_to = "cluster", values_to = "q") %>%
        dplyr::mutate(cluster = as.integer(stringr::str_remove(cluster, "Q")))
    })
  ) %>%
  dplyr::select(K, rep, q_long) %>%
  tidyr::unnest(q_long)
write_tsv(best_q_long, file.path("/home/k14m234/erythrura/results/admixture", "bestQ_long.tsv"))

best_q_wide <- best_q_long %>%
  select(K, rep, IID, species, cluster, q) %>%
  pivot_wider(names_from = cluster, values_from = q, names_prefix = "Q")

write_tsv(best_q_wide, file.path("/home/k14m234/erythrura/results/admixture", "bestQ_wide.tsv"))

# plot scatter
p4 <- runs %>%
  ggplot(aes(x = K, y = cv)) +
  geom_point() +
  geom_line(data = best_runs, linewidth = 0.6) +
  theme_bw() +
  labs(x = "K", y = "CV error", title = "ADMIXTURE cross-validation across replicates")

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/admixture_cv.pdf",
  plot = p4,
  device = cairo_pdf,
  width = 5,
  height = 4,
  units = "in"
)

# order individuals (same as before)
order_tbl <- best_q_long %>%
  group_by(K, IID, species) %>%
  summarise(maxQ = max(q), .groups="drop") %>%
  arrange(K, species, desc(maxQ), IID) %>%
  group_by(K) %>%
  mutate(order = row_number()) %>%
  ungroup()

plot_df <- best_q_long %>%
  left_join(order_tbl, by=c("K","IID","species")) %>%
  filter(K %in% c(2,3,4,5)) %>% 
  mutate(IID = factor(IID, levels = unique(IID[order(order)])))

write_tsv(plot_df, "/home/k14m234/erythrura/scripts/data/admixture.tsv")

p5 <- plot_df %>%
  ggplot(aes(x = IID, y = q, fill = factor(cluster))) +
  geom_col(width = 1, color="white") +
  facet_wrap(~ K, ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.15, "in"),
    strip.background = element_blank(),
  ) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Cluster",
       title = "ADMIXTURE best replicate per K")

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/admixture_k_2_5.pdf",
  plot = p5,
  device = cairo_pdf,
  width = 5,
  height = 6,
  units = "in"
)

# subset for k2
plot_k2 <- plot_df %>% filter(K==2)

# plot k2
p6 <- plot_k2 %>%
  ggplot(aes(x = IID, y = q, fill = factor(cluster))) +
  geom_col(width = 1, color="white") +
  scale_fill_manual(values=c("yellowgreen", "dodgerblue1")) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.15, "in"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
  ) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Cluster")

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/admixture_k_2.pdf",
  plot = p6,
  device = cairo_pdf,
  width = 6,
  height = 3,
  units = "in"
)

# manhattan plot of all samples 
fst <- read_delim("erythrura/results/vcftools/species_fst_w50k_s10k.windowed.weir.fst", delim = "\t")
chr_map <- read_tsv("/home/k14m234/erythrura/scripts/data/nc_to_chr.tsv",col_names = c("CHROM", "chr"))
fst_chr <- fst %>% left_join(chr_map, by = "CHROM")
chr_levels <- c("1","1A",as.character(2:29),"Z")
fst_chr$chr <- factor(fst_chr$chr, levels = chr_levels)
fst_chr <- fst_chr %>%
  group_by(chr) %>%
  mutate(chr_len = max(BIN_END)) %>%
  ungroup()
fst_chr$chr <- factor(fst_chr$chr, levels = chr_levels)
fst_chr <- fst_chr %>%
  arrange(chr, BIN_START) %>%
  mutate(chr_offset = cumsum(chr_len) - chr_len,
         bp_cum = BIN_START + chr_offset)
fst_chr <- fst_chr %>%
  filter(chr != "MT")
axis_df <- fst_chr %>%
  group_by(chr) %>%
  summarize(
    mid = mean(bp_cum),
    chr_size = max(BIN_END)
  ) %>%
  filter(chr_size >= 2.5e7)

pal <- rep(c("yellowgreen", "dodgerblue1"), length.out = length(chr_levels))
thr <- quantile(fst_chr$MEAN_FST, 0.9999, na.rm = TRUE)

p7 <- ggplot(fst_chr, aes(x = bp_cum, y = MEAN_FST, color = chr)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = pal, drop = FALSE) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  scale_x_continuous(breaks = axis_df$mid, labels = axis_df$chr) +
  geom_hline(yintercept = thr,
             linetype = "dashed",
             color = "red",
             linewidth = 0.6,
             alpha=0.6) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  labs(x = "Chromosome", y = "Mean FST")

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/manhattan.pdf",
  plot = p7,
  device = cairo_pdf,
  width = 9,
  height = 4,
  units = "in"
)

# manhattan plot of ingroup
fst_in <- read_delim("erythrura/results/vcftools/ingroup_fst_w50k_s10k.windowed.weir.fst", delim = "\t")
fst_in_chr <- fst_in %>% left_join(chr_map, by = "CHROM")
fst_in_chr$chr <- factor(fst_in_chr$chr, levels = chr_levels)
fst_in_chr <- fst_in_chr %>%
  group_by(chr) %>%
  mutate(chr_len = max(BIN_END)) %>%
  ungroup()
fst_in_chr$chr <- factor(fst_in_chr$chr, levels = chr_levels)
fst_in_chr <- fst_in_chr %>%
  arrange(chr, BIN_START) %>%
  mutate(chr_offset = cumsum(chr_len) - chr_len,
         bp_cum = BIN_START + chr_offset)
fst_in_chr <- fst_in_chr %>%
  filter(chr != "MT")
axis_df_in <- fst_in_chr %>%
  group_by(chr) %>%
  summarize(
    mid = mean(bp_cum),
    chr_size = max(BIN_END)
  ) %>%
  filter(chr_size >= 2.5e7)

p8 <- ggplot(fst_in_chr, aes(x = bp_cum, y = MEAN_FST, color = chr)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = pal, drop = FALSE) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  scale_x_continuous(breaks = axis_df$mid, labels = axis_df$chr) +
  geom_hline(yintercept = thr,
             linetype = "dashed",
             color = "red",
             linewidth = 0.6,
             alpha=0.6) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  labs(x = "Chromosome", y = "Mean FST")

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/manhattan_in.pdf",
  plot = p8,
  device = cairo_pdf,
  width = 9,
  height = 4,
  units = "in"
)

# facetplot
fst_chr$panel     <- "All Samples"
fst_in_chr$panel  <- "Excluding Australia + Palau"
fst_both <- bind_rows(fst_chr, fst_in_chr)
p9 <- ggplot(fst_both, aes(x = bp_cum, y = MEAN_FST, color = chr)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = pal, drop = FALSE) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  scale_x_continuous(breaks = axis_df$mid, labels = axis_df$chr) +
  geom_hline(yintercept = thr,
             linetype = "dashed",
             color = "red",
             linewidth = 0.6,
             alpha=0.6) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  labs(x = "Chromosome", y = "Mean FST") +
  facet_wrap(~panel, ncol = 1)

ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/manhattan_facet.pdf",
  plot = p9,
  device = cairo_pdf,
  width = 9,
  height = 6,
  units = "in"
)

# outliers: all windows
fst_out <- read_tsv("/home/k14m234/erythrura/results/fst/all_windows_with_genes.tsv.gz",
                    show_col_types = FALSE) %>%
  mutate(
    CHROM = as.character(CHROM),
    BIN_START = as.integer(BIN_START),
    BIN_END   = as.integer(BIN_END),
    mid = (BIN_START + BIN_END) / 2
)

# map accessions -> chr labels (your existing mapping)
fst_out_chr <- fst_out %>%
  left_join(chr_map, by = "CHROM") %>%
  filter(!is.na(chr)) %>%                    # drop unmapped contigs if any
  mutate(chr = factor(chr, levels = chr_levels)) %>%
  filter(chr != "MT")                        # optional

# per-chromosome table with lengths + offsets
chr_tbl <- fst_out_chr %>%
  group_by(chr) %>%
  summarise(chr_len = max(BIN_END, na.rm = TRUE), .groups="drop") %>%
  arrange(chr) %>%
  mutate(chr_offset = lag(cumsum(chr_len), default = 0))

# reattach and compute bp_cum safely
fst_out_chr <- fst_out_chr %>%
  left_join(chr_tbl, by = "chr") %>%
  arrange(chr, BIN_START) %>%
  mutate(bp_cum = BIN_START + chr_offset)

axis_df_out <- fst_out_chr %>%
  group_by(chr) %>%
  summarise(
    mid = (min(bp_cum, na.rm = TRUE) + max(bp_cum, na.rm = TRUE)) / 2,
    chr_size = max(BIN_END, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(chr_size >= 2.5e7, !is.na(mid))

# threshold (top 0.1% by default)
thr <- quantile(fst_in_chr$MEAN_FST, 0.9999, na.rm = TRUE)

# candidate list
cands <- read_lines("erythrura/config/candidates_beak.txt") %>%
  str_trim() %>%
  discard(~ .x == "" | str_starts(.x, "#")) %>%
  unique()

# find best hit per candidate gene 
cand_hits <- fst_out_chr %>%
  filter(!is.na(genes), genes != ".") %>%
  separate_rows(genes, sep = ",") %>%
  mutate(genes = str_trim(genes)) %>%
  filter(genes %in% cands) %>%
  group_by(genes) %>%
  slice_max(order_by = MEAN_FST, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(gene = genes) %>%
  mutate(fst_pct = percent_rank(MEAN_FST))

# find outliers
outlier_hits <- fst_out_chr %>%
  filter(!is.na(genes), genes != ".") %>%
  separate_rows(genes, sep = ",") %>%
  mutate(genes = str_trim(genes)) %>%
  filter(MEAN_FST > thr) %>%
  group_by(genes) %>%
  slice_max(order_by = MEAN_FST, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(gene = genes) %>%
  mutate(fst_pct = percent_rank(MEAN_FST))

# plot same aesthetic
p10 <- ggplot(fst_out_chr, aes(x = bp_cum, y = MEAN_FST, color = chr)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = pal, drop = FALSE) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  scale_x_continuous(breaks = axis_df_out$mid, labels = axis_df_out$chr) +
  geom_hline(yintercept = thr,
             linetype = "dashed",
             color = "red",
             linewidth = 0.6,
             alpha = 0.6) +
  geom_point(
    data = cand_hits,
    aes(x = bp_cum, y = MEAN_FST),
    inherit.aes = FALSE,
    shape = 21, color = "black", stroke = 0.5, fill = "orange", 
    size = 3, alpha = 0.9
  ) +
  geom_point(
    data = outlier_hits,
    aes(x = bp_cum, y = MEAN_FST),
    inherit.aes = FALSE,
    shape = 21, color = "black", stroke = 0.5, fill = "purple", 
    size = 3, alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = cand_hits,
    aes(x = bp_cum, y = MEAN_FST, label = gene),
    inherit.aes = FALSE,
    size = 3, max.overlaps = 50
  ) +
  ggrepel::geom_text_repel(
    data = outlier_hits,
    aes(x = bp_cum, y = MEAN_FST, label = gene),
    inherit.aes = FALSE,
    size = 3, max.overlaps = 50
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  labs(x = "Chromosome", y = "Mean FST")

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/manhattan_candidates.pdf",
  plot = p10,
  device = cairo_pdf,
  width = 9,
  height = 4,
  units = "in"
)

# outlier go dotplot 
candidate_genes <- fst_out_chr %>%
  filter(MEAN_FST >= thr, genes != ".", !is.na(genes)) %>%
  separate_rows(genes, sep = ",") %>%
  mutate(genes = str_trim(genes)) %>%
  distinct(genes) %>%
  pull(genes)

# check mapping success
mapped <- bitr(candidate_genes,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Gg.eg.db)

length(unique(mapped$SYMBOL)) / length(candidate_genes)

# filter genes to db
all_genes <- fst_out_chr %>%
  filter(genes != ".", !is.na(genes)) %>%
  separate_rows(genes, sep = ",") %>%
  mutate(genes = str_trim(genes)) %>%
  distinct(genes) %>%
  pull(genes)

# enrich!
ego <- enrichGO(
  gene = candidate_genes,
  universe = all_genes,
  OrgDb = org.Gg.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)
ego

# grep for families
special_genes <- grep("BMP|WNT|FGF|SHH|HOX|DLX|MSX|RUNX", candidate_genes, value = TRUE)
special_hits <- fst_out_chr %>%
  filter(!is.na(genes), genes != ".") %>%
  separate_rows(genes, sep=",") %>%
  mutate(genes = str_trim(genes)) %>%
  filter(genes %in% special_genes)

# check chr clustering
special_hits %>%
  dplyr::select(gene=genes, CHROM, chr, BIN_START, BIN_END, MEAN_FST) %>%
  arrange(CHROM, BIN_START)

# check peak around BMP10
bmp10_region <- fst_out_chr %>%
  filter(CHROM == "NC_042587.1",
         BIN_START >= 3600000,x
         BIN_START <= 5600000)

# check ID of special genes
clean_genes <- candidate_genes[!grepl("^LOC|^TRNA", candidate_genes)]
desc <- AnnotationDbi::select(
  org.Gg.eg.db,
  keys = clean_genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL","GENENAME")
)

# plot around BMP 10
p11 <- ggplot(bmp10_region, aes(BIN_START, MEAN_FST)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab("position") +
  ylab("mean FST")

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/bmp10_window.pdf",
  plot = p11,
  device = cairo_pdf,
  width = 5,
  height = 4,
  units = "in"
)

# check size
bmp10_region %>%
  filter(MEAN_FST > 0.5) %>%
  summarise(
    min_bp = min(BIN_START),
    max_bp = max(BIN_END),
    span = max_bp - min_bp
  )

# check other genes
fst_out_chr %>%
  filter(CHROM == "NC_042587.1",
         BIN_START >= 4300000,
         BIN_START <= 5000000,
         genes != ".") %>%
  dplyr::select(BIN_START, MEAN_FST, genes) %>%
  arrange(desc(MEAN_FST))

# network
gl <- vcfR2genlight(vcf)

# Ensure order matches individuals in VCF/genlight:
samples <- samples[match(indNames(gl), samples$sample_id), ]
pop(gl) <- samples$species
st <- stamppConvert(gl, type = "genlight")  
D <- stamppNeisD(st, pop = FALSE)  # matrix
D <- as.matrix(D)

# write distance matrice
d <- D %>% as.data.frame()
write_tsv(d, "/home/k14m234/erythrura/scripts/data/nj_distance.tsv")

# create, root, and ladderize
nj_tree <- nj(D)
nj_tree <- midpoint(nj_tree)
nj_tree <- ladderize(nj_tree)
write.tree(nj_tree, "/home/k14m234/erythrura/scripts/data/nj_distance.tre")

# shapevals
shape_vals <- c(21, 22, 23, 24, 25)
names(shape_vals) <- unique(samples$country)

# ggtree plot
p11 <- ggtree(nj_tree, layout = "daylight", size = 0.5) %<+% samples +
  geom_tippoint(aes(fill = species, shape = country),
                size = 2.5,
                color = "black",
                stroke = 0.25) +
  #geom_tiplab(size = 2.5, align = TRUE, offset = 0.002) +
  scale_fill_manual(values = pal_named, name = "Species") +
  scale_shape_manual(values = shape_vals, name = "Country") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "white"))
  ) +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/nj_tree.pdf",
  plot = p7,
  device = cairo_pdf,
  width = 7,
  height = 5,
  units = "in"
)

# compare model fit
im <- read_tsv("erythrura/scripts/dadi/im_model.tsv", col_names = FALSE)[,-1]
sc <- read_tsv("erythrura/scripts/dadi/sc_model.tsv", col_names = FALSE)[,-1]
si <- read_tsv("erythrura/scripts/dadi/si_model.tsv", col_names = FALSE)[,-1]
colnames(im) <- c("ll_model", "nTri","nPap","tSplit","m12","m21","theta")
colnames(sc) <- c("ll_model", "nTri","nPap","tSplit","tContact","m12","m21","theta")
colnames(si) <- c("ll_model", "nTri","nPap","tSplit", "theta")
im <- im %>% mutate(model = "IM", AIC = -2*ll_model + 2*6) %>% select(model, AIC, ll_model)
sc <- sc %>% mutate(model = "SC", AIC = -2*ll_model + 2*7) %>% select(model, AIC, ll_model)
si <- si %>% mutate(model = "SI", AIC = -2*ll_model + 2*4) %>% select(model, AIC, ll_model)
mod.df <- bind_rows(im, sc, si)

# get best supported model
best_aic <- mod.df[which.min(mod.df$AIC),]
best_ll <- mod.df[which.max(mod.df$ll_model),]

# select top 20 models
df_top <- mod.df %>%
  group_by(model) %>%
  slice_max(ll_model, n = 20)   # keep top 20 per model

# plot model comparison
p8 <- ggplot(df_top, aes(x=model,y=ll_model)) +
  geom_boxplot() +
  geom_jitter(pch=21, size=2) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  ylab("Log-likelihood") +
  xlab("Model")
p8

# export
ggsave(
  filename = "/home/k14m234/erythrura/scripts/figures/model_fits.pdf",
  plot = p8,
  device = cairo_pdf,
  width = 5,
  height = 4,
  units = "in"
)


