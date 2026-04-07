#####################################
#Gene expression analysis in Fagus crenata
#Full file
#####################################

#creator : Valentin Journ\'e; Kyushu University
#contact : journe.valentin@gmail.com
#file created on 2024-09-26
#updated on 2026 - 02

#####################################
#loading package and functions
#####################################
#https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html
#for the install some packages require other dependencies
#see eg. https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html
source('functionsGeneExp.R')
using(
  'here',
  'pheatmap',
  'tidyverse',
  'readxl',
  "ggdendro",
  'factoextra',
  'NbClust',
  'RColorBrewer',
  'UpSetR',
  'cowplot',
  'patchwork',
  "Boruta",
  "data.table",
  "lubridate",
  "betareg",
  "EnhancedVolcano",
  "limma",
  "ComplexHeatmap",
  "ggimage",
  "ggfortify"
)

library(DHARMa)
library(glmmTMB)
library(emmeans)
library(gghalves)
library(ggsignif)
library(lme4) # for mixed models
library(lmerTest) # adds p-values to lme4 output
#set theme for ggplot
theme_set(theme_cowplot())
#general color option
#colors <- colorRampPalette(brewer.pal(9, "BrBG"))(100) #RdBu PuOr
colors <- colorRampPalette(scico::scico(30, palette = "roma", direction = -1))(
  100
)

colors.roma.short <- colorRampPalette(scico::scico(
  30,
  palette = "roma",
  direction = -1
))(
  10
)

base_cols_treeid <- ggsci::pal_iterm("Rose Pine")(6)
TreeID.color.pal = colorRampPalette(base_cols_treeid)(8)


genes.color <- ggsci::pal_aaas()(10)
genes.color.pal = colorRampPalette(genes.color)(13)


RColorBrewer::brewer.pal(n = 12, name = "BuPu")
colors.pvalues <- colorRampPalette(c(
  "#C6DBEF",
  "#9ECAE1",
  "#6BAED6",
  "#4292C6",
  "#2171B5",
  "#08519C",
  "#08306B"
))(100)
#####################################
#Fagus Gene Expression files loading and cleaning
#####################################
#blast file with tair id and fagus crenata genome
blastparabidopsis = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/blastp_Athaliana/blastp_db_Araport11_201606_q_F.crenata_pep_tophit_summary.xlsx'
) %>%
  mutate(Tair.up = str_extract(`Tair ID`, ".*(?=\\.)"))

#plaza data used later for GO data analysis
#tair id with function
plaza.go.data = read.delim(
  "~/Dropbox/F_crenata/Plaza/go.ath.csv",
  header = FALSE,
  comment.char = "#"
) %>%
  dplyr::select(1, 3, 8) %>%
  dplyr::rename(TAIR.ID = 1, GO_ID = 2, FUNCTION = 3)

#isoforms gene table and made a clean version
#each isoform associated with name
#but do not necessearly match a fagus gene
isoform = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/ISO-seq/final/Facr_r1.0_s1_collapse_isoforms_ann.xlsx',
  skip = 1
)

isoform.cleand = isoform %>%
  mutate(name = qseqid) %>%
  dplyr::select(name, ProteinName, sseqid, Subject...58) %>%
  mutate(
    gene = case_when(
      str_detect(Subject...58, "No_hits|#") ~ NA_character_, # convert to NA
      TRUE ~ str_remove(Subject...58, " unnamed protein product$") # remove suffix if present
    )
  ) %>%
  dplyr::select(-sseqid) %>%
  rename(sseqid = gene)


#the monthly year gene expression for fagus individual
#use data table, faster and less memory usage
genes_rpm = data.table::fread(
  "/Users/valentinjourne/Dropbox/F_crenata/RNAseqdata/final/rpm.gn.csv"
)
setnames(genes_rpm, gsub("//", "..", names(genes_rpm), fixed = TRUE))
setnames(genes_rpm, "V1", "qseqid")

#dplyr syntax to format table
#here change index to get information about tree , samples, etc (because it includes date)
data.genes.all <- genes_rpm %>%
  as_tibble() %>%
  pivot_longer(
    -qseqid,
    names_to = 'IndexNo',
    values_to = 'level.exp'
  ) %>%
  mutate(
    date = str_extract(IndexNo, "\\d{6}"),
    TreeID = str_extract(IndexNo, "(?<=Fc)\\d+"),
    Tissue = str_extract(IndexNo, "(L|B|R|R1|R2|R3|R4)$"), # ⬅️ extract last letter as tissue type with R1 and other roots
    dateformat = as.Date(date, format = "%y%m%d"),
    month = month(dateformat),
    year = year(dateformat)
  )

length(unique(data.genes.all$qseqid))

#determine sum gene expression
#do this to remove low gene expression
#i will do this process for all, even if I will focus only on leaves
sumexp = data.genes.all %>%
  group_by(qseqid, Tissue) %>%
  summarise(sum = sum(level.exp), mean = mean(level.exp))

hist(log(sumexp$mean))
summary(sumexp)
summary(sumexp$sum)
sumexp %>%
  ggplot(aes(x = log2(mean))) +
  geom_histogram() +
  facet_grid(. ~ Tissue) +
  ggpubr::theme_cleveland()

breakpoints <- c(0, 1, 10, 100, 1000, Inf)

sumexp %>%
  mutate(
    categories_gene_exp = cut(
      mean,
      breaks = breakpoints,
      labels = c('0-1', '1-10', '10-100', '100-1000', 'Above 1000'),
      right = FALSE
    )
  ) %>%
  group_by(categories_gene_exp, Tissue) %>%
  summarise(count = n(), percentage = (count / nrow(sumexp)) * 100) %>%
  ggplot(aes(
    x = categories_gene_exp,
    y = percentage,
    fill = categories_gene_exp
  )) +
  geom_bar(stat = 'identity') +
  geom_text(
    aes(label = paste0(count, " (", round(percentage, 1), "%)")),
    vjust = -0.5
  ) +
  facet_grid(. ~ Tissue) +
  labs(x = "Category", y = "Percentage") +
  scale_fill_brewer(palette = "Set1") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 90))

#see how many genes I would keep below thrshould
thresholds <- c(seq(0.1, 1, by = 0.1), seq(1, 10, by = 1))

sumexp.subset = sumexp %>% filter(Tissue == "L")
#number total genes reported for leaves
dim(sumexp.subset)

removal_summary <- thresholds %>%
  purrr::map_dfr(
    ~ {
      tibble(
        threshold = .x,
        removed = sum(sumexp.subset$mean < .x),
        kept = sum(sumexp.subset$mean >= .x),
        total = nrow(sumexp.subset)
      )
    }
  ) %>%
  mutate(percentage_removed = 100 * removed / total)

removal_summary

ggplot(sumexp.subset, aes(x = mean)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  scale_x_continuous(trans = "log10") +
  geom_vline(xintercept = thresholds, color = "red", linetype = "dashed") +
  geom_text(
    data = data.frame(x = thresholds),
    aes(x = x, y = 4000, label = paste0("<", x)),
    color = "gray30",
    angle = 90,
    vjust = 0,
    hjust = 0,
    size = 3
  )

#now data to remove
to.remove.below.threshold = sumexp %>%
  filter(Tissue == "L") %>%
  filter(!mean < 1) %>%
  dplyr::select(-sum)

#####################################
#Fagus Gene Expression - focus on leaves
#####################################
#now get cumsum for genes
#focus on leaves, and remove october because those are missing samples for half years

data.genes.leaves.monthly <- data.genes.all %>%
  filter(Tissue == "L") %>%
  right_join(
    to.remove.below.threshold %>% filter(Tissue == "L"),
    by = c("Tissue", "qseqid")
  ) %>%
  filter(month != 10) %>%
  mutate(
    level.exp.log10 = log10(0.25 + level.exp),
    level.exp.log2 = log2(1 + level.exp)
  ) %>%
  arrange(qseqid, TreeID, year, month) %>%
  group_by(TreeID, year, qseqid) %>%
  mutate(
    cumsum_log2 = cumsum(level.exp.log2),
    cumsum_log10 = cumsum(level.exp.log10)
  ) %>%
  ungroup()

table(data.genes.leaves.monthly$qseqid)
head(data.genes.leaves.monthly)
#count number of genes before and after filtering
length(unique(data.genes.all$qseqid))
length(unique(data.genes.leaves.monthly$qseqid))

#that will be the dataset
#write_csv(
#  data.genes.leaves.monthly,
#  here("data_clean", 'data.genes.leaves.monthly.csv')
#)

#####################################
#Flowering dataset, loading, cleaning and plot
#####################################
flo.intensity = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/FlowerIntensity.xlsx',
  skip = 4,
  n_max = 93
)
# Function to replace the string in a column
string_to_replace = '<5'
new_value = 4

flo.intensity = flo.intensity %>%
  mutate_at(vars(-`Tree ID`), list(replace_string)) %>%
  mutate_at(vars(-`Tree ID`), as.numeric) %>%
  mutate(`Tree ID` = as.character(`Tree ID`)) %>%
  rename(TreeID = `Tree ID`) %>%
  dplyr::filter(
    TreeID %in% c('60', '202', '206', '224', '263', '264', '209', '210')
  ) %>%
  pivot_longer(
    `2013`:`2023`,
    names_to = 'year',
    values_to = 'flowering.percentage'
  ) %>%
  mutate(year = as.numeric(as.character(year)))

#here I am shifting to one year, because it is the year of flowering, but I know that initiation
#occured the previous year
plotflointensity = flo.intensity %>%
  filter(year > 2013) %>%
  mutate(year = year, year = as_factor(year)) %>%
  ggplot(aes(
    x = year,
    y = flowering.percentage,
    group = TreeID,
    col = TreeID,
    fill = TreeID
  )) +
  geom_line(alpha = .5) +
  geom_point(shape = 21) +
  ylab('Flowering intensity (%)') +
  xlab('') +
  theme(legend.position = c(.8, .8)) +
  scale_color_manual(values = TreeID.color.pal, "Tree") +
  scale_fill_manual(values = TreeID.color.pal, "Tree") #+
#stat_summary(
#  aes(group = 1),
#  fun = mean,
#  geom = "line",
#  linewidth = 1,
#  color = "black"
#)

plotflointensity

cowplot::save_plot(
  'figuresR/seedproduction.pdf',
  plotflointensity,
  nrow = 1,
  ncol = 1
)

plot.inititation = flo.intensity %>%
  filter(year > 2014) %>%
  mutate(year = year - 1, year = as_factor(year)) %>%
  ggplot(aes(
    x = year,
    y = flowering.percentage,
    group = TreeID,
    col = TreeID,
    fill = TreeID
  )) +
  geom_line(alpha = .5) +
  geom_point(shape = 21) +
  ylab('Initiation intensity (%)') +
  xlab('') +
  theme(legend.position = c(.8, .8)) +
  scale_color_manual(values = TreeID.color.pal, "Tree") +
  scale_fill_manual(values = TreeID.color.pal, "Tree") +
  theme(legend.position = "none")

plot.inititation

cowplot::save_plot(
  'figuresR/plot.inititation.pdf',
  plot.inititation,
  nrow = .8,
  ncol = 1.5 / 2
)


#convert to binary mast vs non mast years
#and shift to one year because flowering is the year of flowering but initiation is the previous year
#adn this file will be use later
#I did presence absence, but using a threshold of 5 percent will not change much the results, because most values are either 0 or above 5
average.flo.intensity.individual = flo.intensity %>%
  group_by(year, TreeID) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(year = year - 1) %>%
  ungroup() %>%
  mutate(
    mastONOFF = ifelse(flowering.percentage > 0, 1, 0),
    fac.mastONOFF = as_factor(mastONOFF),
    flowering.percentage.trans = y.transf.betareg(flowering.percentage) / 100
  ) %>%
  ungroup()

#####################################
#WEATHER AND CLIMATE dataset, loading, cleaning and plot
#####################################
#create temperature data summary
temperature.station = read_csv(
  '/Users/valentinjourne/Dropbox/F_crenata/climate/NaebaClimate2011-2023.csv',
  skip = 5
) %>%
  format_climate_temperature()

data.summary.temperature.month = temperature.station %>%
  filter(month %in% c(6, 7, 8, 9, 10)) %>%
  group_by(year, month) %>%
  summarise(
    mean.temp = mean(meanTemp, na.rm = T),
    sd.temp = sd(meanTemp, na.rm = T),
    .groups = "drop"
  ) %>%
  filter(year > 2013 & year < 2023)


# create a continuous time variable
data.summary.temperature.month <- data.summary.temperature.month %>%
  mutate(
    ym = interaction(year, month, sep = "-"),
    ym = factor(ym, levels = unique(ym))
  )

daily <- temperature.station %>%
  filter(year > 2013, year < 2023, month %in% 6:10)

daily_full = daily %>%
  mutate(Date = as.Date(Date)) %>%
  complete(
    Date = seq(min(Date), max(Date), by = "day")
  ) %>%
  mutate(
    year = year(Date),
    month = month(Date)
  ) %>%
  filter(year > 2013, year < 2023, month %in% 6:10)

daily_roll <- daily %>%
  arrange(date.bis) %>%
  mutate(
    temp_roll15days = slider::slide_dbl(
      meanTemp,
      mean,
      .before = 7,
      .after = 7,
      .complete = TRUE,
      na.rm = TRUE
    )
  ) %>%
  complete(
    date.bis = seq(min(date.bis), max(date.bis), by = "day")
  ) %>%
  mutate(
    year = year(date.bis),
    month = month(date.bis)
  ) %>%
  filter(year > 2013, year < 2023, month %in% 6:10) %>%
  mutate(
    year = year(date.bis),
    season_day = row_number(), # continuous index after filtering
    year_f = factor(year)
  )

breaks_df <- daily_roll %>%
  group_by(year) %>%
  summarise(season_day0 = min(season_day), .groups = "drop")


daily.temp = ggplot(daily_roll) +
  geom_line(
    aes(season_day, meanTemp, group = year_f),
    alpha = .9,
    color = "darkblue",
    linewidth = .3
  ) +
  geom_line(
    aes(season_day, temp_roll15days, group = year_f),
    linewidth = 1.2,
    color = "deepskyblue"
  ) +
  scale_x_continuous(
    breaks = breaks_df$season_day0,
    labels = breaks_df$year,
    expand = c(0, 0)
  ) +
  xlab("") +
  ylab("Daily T (°C)\n")
daily.temp

cowplot::save_plot(
  'figuresR/daily_temp.pdf',
  daily.temp,
  nrow = .5,
  ncol = 2
)

#now make anomalies plot

temp_anom_year <- temperature.station %>%
  filter(month %in% c(6, 7, 8, 9, 10)) %>%
  group_by(year, month) %>%
  summarise(
    mean.temp = mean(meanTemp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  summarise(
    summer_mean = mean(mean.temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    clim_mean = mean(summer_mean, na.rm = TRUE), # mean over all years
    anomaly = summer_mean - clim_mean
  )

ano.year.temp = ggplot(
  temp_anom_year %>% filter(year > 2013 & year < 2023),
  aes(x = factor(year), y = anomaly, fill = anomaly > 0)
) +
  geom_col(width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("steelblue", "firebrick"),
    labels = c("Cooler", "Warmer"),
    name = ""
  ) +
  labs(
    x = "",
    y = "Growing season \nT anom. (°C)"
  ) +
  theme(
    axis.text.x = element_text(vjust = .5),
    legend.position = c(.7, .4)
  ) +
  ylim(-1, .5)
ano.year.temp

cowplot::save_plot(
  'figuresR/anomalie_temp.pdf',
  ano.year.temp,
  nrow = .5,
  ncol = 2.2
)

#just small betaregregresison
small.model.flo.anomalies = average.flo.intensity.individual %>%
  left_join(temp_anom_year) %>%
  mutate(
    flo.mean.y = y.transf.betareg(flowering.percentage / 100),
    factor.TreeID = as.factor(TreeID)
  )
summary(glmmTMB::glmmTMB(
  flo.mean.y ~ anomaly + (1 | factor.TreeID),
  data = small.model.flo.anomalies,
  family = glmmTMB::beta_family()
))


#####################################
#Data combination and analysis for seasonal variation - gene by gene using clustering and heatmaps
#####################################
#here I am creating a new one because I want to aggregate for year-month to make heatmaps
#and pca analysis
summary_all_genes_monthly <- data.genes.leaves.monthly %>%
  arrange(qseqid, TreeID, year, month) %>%
  group_by(year, month, qseqid, dateformat) %>%
  summarise(mean.log2.month.year = mean(level.exp.log2)) %>%
  ungroup()

#summarise per isoform
isoform_summary <- summary_all_genes_monthly %>%
  group_by(qseqid) %>%
  summarise(
    sum_log2 = sum(mean.log2.month.year, na.rm = TRUE),
    .groups = "drop"
  ) %>% #here add this to keep the higher sum2 log
  arrange(qseqid, desc(sum_log2)) %>%
  group_by(qseqid) %>%
  mutate(
    rank = row_number(),
    diff_to_top = sum_log2 - first(sum_log2)
  )

isoform_summary

#caclulate timpoint and remove october from heatmaps
wider.monthly.expression.allgenes = summary_all_genes_monthly %>%
  filter(month != 10) %>%
  group_by(qseqid) %>%
  mutate(sum0 = sum(mean.log2.month.year)) %>%
  filter(sum0 > 0) %>%
  dplyr::select(-sum0) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  ) %>%
  dplyr::select(timepoint, qseqid, mean.log2.month.year) %>% #sseqid
  pivot_wider(
    names_from = timepoint,
    values_from = c(mean.log2.month.year)
  ) %>%
  as.data.frame()

row.names(
  wider.monthly.expression.allgenes
) <- wider.monthly.expression.allgenes$qseqid
wider.monthly.expression.allgenes$qseqid <- NULL

data_subset_norm.correlation_monthly_gene <- t(apply(
  as.matrix(wider.monthly.expression.allgenes),
  1,
  cal_z_score
))

#determine optimal k for cluster
###############################
#because it is crashing just randomly select a subset (30 percent)
set.seed(123)
sub_idx <- sample(
  1:nrow(data_subset_norm.correlation_monthly_gene),
  round(dim(data_subset_norm.correlation_monthly_gene)[1] * 0.30)
)
data_small <- data_subset_norm.correlation_monthly_gene[sub_idx, ]

d <- dist(t(data_small))
sil_width <- sapply(2:10, function(k) {
  pam_fit <- cluster::pam(d, diss = TRUE, k = k) #transpose or not the matrix (for time clus or gene)
  pam_fit$silinfo$avg.width
})
plot(2:10, sil_width, type = "b", xlab = "k", ylab = "Average silhouette width")

display_labels_pheatmap.all <- rep(c("6", "7", "8", "9"), times = 9)
#try complex map
###############################
#here I tried to make subset before using the full function
#mat = data_subset_norm.correlation_monthly_gene[
#  sample(nrow(data_subset_norm.correlation_monthly_gene), 1000),
#]
mat <- data_subset_norm.correlation_monthly_gene
years <- factor(rep(2014:2022, each = 4))
months <- rep(6:9, times = 9)
ht_opt$message = FALSE

ht <- Heatmap(
  mat,
  name = "z values", # legend title = "z"
  col = colors,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_km = 4, # like cutree_rows = 2, i initially put two, but I checked and seems that 4 better
  column_split = years, # 4-column blocks by year
  column_labels = as.character(months), # show 6,7,8,9 under each column
  column_gap = unit(.5, "mm"),
  show_row_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_dend_reorder = FALSE,
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "z values",
    title_position = "topcenter",
    legend_direction = "horizontal" # Makes the color bar horizontal
  )
)
#draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
pdf("figuresR/pheatmap.all.pdf", width = 9, height = 3)
draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()


#data for Princp Comp Analysis (PCA)
####################################
#explore gene seasonality for individual and seasons
pca.data.genes = data.genes.all %>%
  filter(Tissue == "L" & month != 10) %>% #
  right_join(
    to.remove.below.threshold %>% filter(Tissue == "L"),
    by = c("Tissue", "qseqid")
  ) %>%
  mutate(
    level.exp.log2 = log2(1 + level.exp)
  ) %>%
  arrange(qseqid, TreeID, year, month) %>%
  drop_na(qseqid) %>%
  right_join(
    isoform_summary %>%
      filter(diff_to_top == 0) %>%
      dplyr::select(qseqid)
  ) %>%
  left_join(
    #now add temperature station matching the date of sampling
    temperature.station %>%
      dplyr::select(
        minTemp,
        maxTemp,
        rainfall,
        sunshine,
        month,
        year,
        date.bis
      ) %>%
      dplyr::rename(date = date.bis) %>%
      mutate(date = as.Date(date)) %>%
      mutate(date = anytime::anydate(date)),
    by = c("dateformat" = "date", "month", "year")
  ) %>%
  ungroup() %>%
  dplyr::select(
    qseqid,
    TreeID,
    dateformat,
    year,
    month,
    level.exp.log2,
    minTemp,
    rainfall,
    sunshine
  ) %>%
  pivot_wider(names_from = qseqid, values_from = level.exp.log2) %>%
  mutate(fac.cat = paste0(TreeID, "_", dateformat))

#do PCA to see cluster genes seasons
expr_matrix_pca <-
  pca.data.genes %>%
  dplyr::select(Facr_v2.5_s1cl000011:tail(names(.), 1)) %>%
  drop_na() %>%
  dplyr::select(-fac.cat)

#PCA for seasonality
pca <- prcomp(expr_matrix_pca, scale. = TRUE)
season.pca.1 = autoplot(
  #based on ggforityfy
  pca,
  data = pca.data.genes %>%
    left_join(
      average.flo.intensity.individual %>%
        dplyr::select(TreeID, year, flowering.percentage)
    ) %>%
    mutate(month = as_factor(month)) %>%
    mutate(
      year = as_factor(year),
      mastyear = if_else(year %in% c(2014, 2017, 2021), "Mast", "Non-mast")
    ),
  colour = 'month',
  size = "flowering.percentage",
  x = 1,
  y = 2
) +
  labs(shape = NULL, size = "Flo. intensity (%)") +
  scale_color_brewer("Months", palette = "BuPu") +
  scale_fill_brewer("Months", palette = "BuPu") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  scale_size_continuous(
    breaks = c(5, 50, 100),
    limits = c(0, 100) # optional but recommended
  ) +
  theme(
    legend.position = c(.6, .8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
season.pca.1

cowplot::save_plot(
  'figuresR/pca.seasons.mast.pdf',
  season.pca.1,
  nrow = 1.2,
  ncol = .8
)

#####################################
#Variability of gene expression
#####################################
#calculate individual CV variablity
data.genes.leaves.monthly
cv_gene = functioncalc_gene_variability_cv(data.genes.leaves.monthly)
#simple plot check comp CV vs kCV
ggplot(cv_gene, aes(x = CV, y = kCV)) + geom_point()
#now calc sync
synchrony_gene = functioncalc_gene_synchrony(
  data.genes.leaves.monthly,
  method = "pearson"
)

gene.cv.syncro = cv_gene %>%
  left_join(synchrony_gene, by = "qseqid") %>%
  mutate(combined.metric = scale(kCV) + scale(abs(mean_synchrony))) #%>%
#left_join(
#  subset.best.logistic.reg %>%
#    dplyr::select(-nobs) %>%
#    filter(term == "cumsum_expr")
#)

#another plots
gene.cv.syncro %>%
  ggplot(aes(x = kCV)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = "Coefficient of Variation (kCV)",
    y = "Number of Genes"
  ) +
  theme_minimal()

gene.cv.syncro %>%
  ggplot(aes(x = mean_synchrony)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = "Synchorny (pearson mean)",
    y = "Number of Genes"
  ) +
  theme_minimal()

#determine band of auqntile low and high 10 percent
low_cut <- quantile(gene.cv.syncro$kCV, 0.1, na.rm = TRUE)
high_cut <- quantile(gene.cv.syncro$kCV, 0.9, na.rm = TRUE)
#low_cut_sync <- quantile(gene.cv.syncro$mean_synchrony, 0.025, na.rm = TRUE)
#high_cut_sync <- quantile(gene.cv.syncro$mean_synchrony, 0.975, na.rm = TRUE)

kCVdis = ggplot(gene.cv.syncro, aes(x = kCV)) +
  geom_histogram(bins = 50, fill = "gray80", color = "white") +
  geom_vline(
    xintercept = c(low_cut, high_cut),
    color = "red",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = low_cut,
    y = 1500,
    label = "Low kCV",
    vjust = 1.5,
    color = "red",
    angle = 90
  ) +
  annotate(
    "text",
    x = high_cut,
    y = 1500,
    label = "High kCV",
    vjust = 1.5,
    color = "red",
    angle = 90
  ) +
  labs(
    x = "kCV",
    y = "Number of Genes"
  ) +
  theme_minimal()

syncdis = ggplot(gene.cv.syncro, aes(x = mean_synchrony)) +
  geom_histogram(bins = 50, fill = "gray80", color = "white") +
  geom_vline(
    xintercept = .3, #c(low_cut_sync, high_cut_sync),
    color = "red",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.3,
    y = 1500,
    label = "Threshold sync",
    vjust = 1.5,
    color = "red",
    angle = 90
  ) +
  labs(
    x = "Sync",
    y = "Number of Genes"
  ) +
  theme_minimal()

kCVdis + syncdis


#####################################
#meaning BVOC / Flowering genes
#####################################
meaningFloGenLFTC = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/RNAseqdata/F.crenata_Flowering_BVOC_genes_summary_geneID.xlsx',
  sheet = 2
)
meaningFloGenLFTC.function = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/RNAseqdata/F.crenata_Flowering_BVOC_genes_summary_geneID.xlsx',
  sheet = 1
) %>%
  dplyr::select(Gene_name, At_gene_id, `F. crenata`, schemes) %>%
  drop_na(Gene_name) %>%
  dplyr::rename(
    Gene_name_Original = Gene_name,
    F.crenata = `F. crenata`,
    AT_gene_id = At_gene_id
  )


flowering.genes.tair.qseqid = blastparabidopsis %>%
  dplyr::select(gene, Tair.up) %>%
  right_join(
    meaningFloGenLFTC.function %>%
      mutate(F.crenata = as.numeric(F.crenata)) %>%
      filter(F.crenata > 0),
    by = c("Tair.up" = "AT_gene_id")
  ) %>% #limit to flo genes
  dplyr::rename(sseqid = gene) %>%
  left_join(
    isoform.cleand %>%
      dplyr::rename(qseqid = 1) %>%
      dplyr::select(qseqid, sseqid) %>%
      distinct()
  ) %>%
  group_by(Gene_name_Original) %>%
  mutate(n = n()) %>%
  mutate(
    Gene_name_Original_new = ifelse(
      n > 1,
      paste0(Gene_name_Original, "_", row_number()),
      Gene_name_Original
    )
  )


#calculate the metric for flowering time series now
cv_flowering <- average.flo.intensity.individual %>%
  group_by(TreeID) %>%
  summarise(
    CVi = my_cv_fun(flowering.percentage),
    .groups = "drop"
  )
cv_flowering_summary <- cv_flowering %>%
  summarise(
    CV = mean(CVi, na.rm = TRUE),
    CV_median = median(CVi, na.rm = TRUE),
    nobs = sum(!is.na(CVi))
  ) %>%
  mutate(
    kCV = sqrt(CV^2 / (1 + CV^2))
  )
flowering_mat <- average.flo.intensity.individual %>%
  dplyr::select(TreeID, year, flowering.percentage) %>%
  pivot_wider(names_from = TreeID, values_from = flowering.percentage)
mat.flo <- cor(
  flowering_mat[, -1],
  use = "pairwise.complete.obs",
  method = "pearson"
)
flowering_synchrony <- tibble(
  mean_synchrony = mean(mat.flo[lower.tri(mat)], na.rm = TRUE),
  percentage.na = sum(is.na(mat.flo)) / prod(dim(mat.flo))
)

#now I got both metrics for CV and sync of flower production
flowering_synchrony
cv_flowering_summary

cv_sync_toghether <- gene.cv.syncro %>%
  mutate(
    Category_kCV = case_when(
      kCV <= low_cut ~ "Low_kCV",
      kCV >= high_cut ~ "High_kCV",
      TRUE ~ "Middle"
    ),
    Category_sync = case_when(
      abs(mean_synchrony) <= 0.3 ~ "Low_sync",
      abs(mean_synchrony) >= 0.3 ~ "High_sync",
      TRUE ~ "Middle"
    )
  ) %>%
  left_join(flowering.genes.tair.qseqid)

#write_csv(
#  cv_sync_toghether %>% dplyr::select(-combined.metric),
#  here::here("data_clean", "cv_sync_together.csv")
#)

Figure.variability.sync.gene.exp = ggplot(
  cv_sync_toghether,
  aes(y = kCV, x = mean_synchrony)
) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  geom_point(
    data = cv_sync_toghether %>%
      #filter(str_detect(schemes, "Flower|flower")) %>%
      drop_na(Gene_name_Original_new),
    aes(col = Gene_name_Original_new)
  ) +
  geom_text_repel(
    data = cv_sync_toghether %>%
      #filter(str_detect(schemes, "Flower|flower")) %>%
      drop_na(Gene_name_Original_new),
    aes(label = Gene_name_Original_new, col = Gene_name_Original_new),
    size = 3,
    max.overlaps = 50,
    fontface = "bold"
  ) +
  guides(color = "none") +
  xlab("Average synchrony") +
  ylab("Coefficient of variation (kCV)") +
  theme(legend.position = "top") +
  geom_point(
    data = tibble(
      kCV = cv_flowering_summary$kCV,
      mean_synchrony = flowering_synchrony$mean_synchrony
    ),
    col = "red",
    shape = 8,
    size = 4
  )

Figure.variability.sync.gene.exp

cowplot::save_plot(
  'figuresR/Figure.variability.sync.gene.exp.pdf',
  Figure.variability.sync.gene.exp,
  nrow = 1.6,
  ncol = 1
)

#####################################
#Data combination and analysis for seasonal variation in relation to flowering - gene by gene using regression
#####################################
#now make some space
rm(genes_rpm)
gc()
# order data
setorder(data.genes.leaves.monthly, TreeID, year, qseqid, month)
unique(data.genes.leaves.monthly$qseqid)

# Create summary of cumsum or mean for each window
#take 1h min at least, MORE!
# Split gene list into batches
#CEHCK TO BE SURE TO HAVE NO GROUPING FACTOR VAR
gene_batches <- data.genes.leaves.monthly %>%
  distinct(qseqid) %>%
  pull(qseqid) %>%
  split(., ceiling(seq_along(.) / 10))
gc()
#try to make it faster by using furrr and parallelization

#ALREDY DONE
batchgene = F
#if you want, but it would take almost a small week to get the results
if (batchgene == T) {
  library(furrr)
  plan(multisession, workers = parallel::detectCores() - 4)
  future::nbrOfWorkers()
  options(future.globals.maxSize = 5 * 1024^3) # Allow 4 GB

  #should not take more than a day with 5000 batch - so please check your gene_batches length
  future_walk2(
    gene_batches,
    seq_along(gene_batches),
    ~ {
      result <- compute_window_summary(.x, data.genes.leaves.monthly)
      qs::qsave(
        result,
        here::here(paste0(
          "format_gene_month_expression_combination/window_summary_batch_",
          .y,
          ".qs"
        ))
      )
    },
    .progress = TRUE
  )
  future::plan(future::sequential)
  future::nbrOfWorkers()
  gc()

  #just load to check some error
  window_summary_all <- list.files(
    path = here("format_gene_month_expression_combination/"),
    pattern = ".qs",
    full.names = TRUE
  ) %>%
    map_dfr(qs::qread)

  #took already a small day
  #betareg is hell
  #https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/11-0029.1
  #the problem I got with betareg is that most values are either 0 (or closed) or 1
  #so many models failed to converge ...
  #so I fitted both here betareg and log reg, but I am more confident with logistic
  #because of the data structure
  #small note, I used the cumsum log2 expression as predictor
  #but mean would provide the same model performance, just slope would change
  ################################################################################
  library(furrr)
  #plan(multisession, workers = parallel::detectCores() - 2)
  # future_walk(
  #   1:length(list.files("format_gene_month_expression_combination")),
  #   ~ run_regression_batch(
  #     batch_index = .x,
  #     input_dir = here::here("format_gene_month_expression_combination/"),
  #     output_dir = here::here("summary_betaregstats_month_geneexp/"),
  #     average.flo.intensity.individual = average.flo.intensity.individual,
  #     modeltype = "betareg",
  #     response_formula = as.formula("flowering.percentage.trans ~ cumsum_expr")
  #   ),
  #   .progress = TRUE
  # )
  # future::plan(future::sequential)
  # future::nbrOfWorkers()
  # gc()

  #test now with logistic regression
  #take 14 hours with 8 cores
  plan(multisession, workers = parallel::detectCores() - 2)
  future_walk(
    1:length(list.files("format_gene_month_expression_combination")),
    ~ run_regression_batch(
      batch_index = .x,
      input_dir = here::here("format_gene_month_expression_combination/"),
      output_dir = here::here("summary_logregstats_month_geneexp/"),
      average.flo.intensity.individual = average.flo.intensity.individual,
      modeltype = "logistic", #if not betareg
      response_formula = as.formula("fac.mastONOFF ~ cumsum_expr")
    ),
    .progress = TRUE
  )
  future::plan(future::sequential)
  future::nbrOfWorkers()
  gc()
}

#I fitted all models at that time so here
#I will focus only on consecutive windows
##########################################
correct_window = c(
  '6',
  '7',
  '8',
  '9',
  '6-7',
  '7-8',
  '8-9',
  '6-7-8',
  '7-8-9',
  '6-7-8-9'
)

#load beta reg outputs and logistic regression outputs
######################################################
#a small fraction wehre not able to get all windows
#that came from the fact that most of the data is either 0 or 1 so betaregression remained unstable and did not converge,
#so I will focus on logistic regression for now
# betareg.fruit.genes.all <- list.files(
#   path = here("summary_betaregstats_month_geneexp/"),
#   pattern = ".qs",
#   full.names = TRUE
# ) %>%
#   map_dfr(qs::qread)
#here I want to keep the top 10 AIC from each sequnce of genes
#then if difference between min and max is below 2, just get the one with highest pseudo R2
#and combine it again to the slope and other stuff
# subset.best.beta.reg = betareg.fruit.genes.all %>%
#   right_join(to.remove.below.threshold %>% dplyr::select(qseqid)) %>%
#   filter(window_months %in% correct_window) %>%
#   group_by(qseqid) %>%
#   arrange(AIC) %>%
#   slice_head(n = 10) %>% #keep the first top 10
#   summarise(
#     best.windows = if (max(AIC) - min(AIC) <= 2) {
#       window_months[which.max(pseudo.r.squared)]
#     } else {
#       window_months[which.min(AIC)]
#     },
#     .groups = "drop"
#   ) %>%
#   left_join(
#     betareg.fruit.genes.all,
#     by = c("qseqid", "best.windows" = "window_months")
#   )

# subset.best.beta.reg %>%
#   dplyr::select(qseqid, best.windows, pseudo.r.squared, AIC) %>%
#   distinct() %>%
#   ggplot(aes(x = pseudo.r.squared, y = AIC, col = best.windows)) +
#   geom_point()
#
# subset.best.beta.reg %>%
#   filter(term == "cumsum_expr") %>%
#   ggplot(aes(x = pseudo.r.squared, y = log10(p.value), col = best.windows)) +
#   geom_point()

#once good month have been identified togetner after the filtering
#split month name to combine it aso keep how many month have veen used
# subset.best.month.reg <- subset.best.beta.reg %>%
#   mutate(
#     months_list = str_split(best.windows, "-", simplify = FALSE) %>%
#       map(as.integer),
#     n_months = map_int(months_list, length)
#   )

#filtering based on logistic regression
######################################################
logreg.fruit.genes.all <- list.files(
  path = here("summary_logregstats_month_geneexp/"),
  pattern = ".qs",
  full.names = TRUE
) %>%
  map_dfr(qs::qread)

subset.best.logistic.reg = logreg.fruit.genes.all %>%
  right_join(to.remove.below.threshold %>% dplyr::select(qseqid)) %>% #add genes with higher expression level
  filter(window_months %in% correct_window) %>%
  group_by(qseqid) %>%
  arrange(AIC) %>%
  slice_head(n = 10) %>% #keep the first top 10
  summarise(
    best.windows = if (max(AIC) - min(AIC) <= 2) {
      window_months[which.max(AUC)]
    } else {
      window_months[which.min(AIC)]
    },
    .groups = "drop"
  ) %>%
  left_join(
    logreg.fruit.genes.all,
    by = c("qseqid", "best.windows" = "window_months")
  )

# comparison.win.log.beta.reg. = subset.best.logistic.reg %>%
#   dplyr::select(qseqid, best.windows, AUC, AIC) %>%
#   distinct() %>%
#   left_join(
#     subset.best.beta.reg %>%
#       dplyr::select(qseqid, best.windows, pseudo.r.squared, AIC) %>%
#       distinct() %>%
#       dplyr::rename(
#         AIC_beta = AIC,
#         pseudo_r2 = pseudo.r.squared,
#         best_windows_beta = best.windows
#       )
#   ) %>%
#   mutate(agreement = ifelse(best.windows == best_windows_beta, 1, 0))
#
# comparison.win.log.beta.reg.$agreement %>%
#   table()

count.best.seson.model = subset.best.logistic.reg %>%
  filter(term == "cumsum_expr") %>%
  group_by(best.windows) %>%
  tally() %>%
  ggplot(aes(n, x = best.windows, col = best.windows, fill = best.windows)) +
  geom_col() +
  coord_flip() +
  scale_color_viridis_d(option = "mako") +
  scale_fill_viridis_d(option = "mako") +
  theme(legend.position = "none") +
  xlab("Best windows") +
  ylab("Count")

subset.best.logistic.reg %>%
  filter(term == "cumsum_expr") %>%
  group_by(best.windows) %>%
  summarise(
    median.AUC = median(AUC),
    median.AIC = median(AIC),
    n = n()
  ) %>%
  arrange(median.AUC)

model.fit.aucaiclog = ggplot(
  subset.best.logistic.reg,
  aes(x = AIC, y = AUC, col = best.windows, fill = best.windows)
) +
  geom_point(alpha = .3, shape = 21, stroke = 1, size = 1) +
  #geom_hline(yintercept = 0.5, col = "black", linetype = "dashed") +
  geom_hline(yintercept = 0.8, col = "black", linetype = "dashed") +
  labs(color = "Best windows", fill = "Best windows") +
  ylab("Area under the ROC curve (AUC)") +
  xlab("Akaike information criterion (AIC)") +
  scale_color_viridis_d(option = "mako") +
  scale_fill_viridis_d(option = "mako") +
  theme(legend.position = "bottom")

FigureSup.model.fit = model.fit.aucaiclog +
  count.best.seson.model
FigureSup.model.fit

cowplot::save_plot(
  'figuresR/FigureSup.model.fit.pdf',
  FigureSup.model.fit,
  nrow = 1.8,
  ncol = 1.8
)
#I Save also as png for easier viewing on inkskape
cowplot::save_plot(
  'figuresR/FigureSup.model.fit.png',
  FigureSup.model.fit + theme(legend.position = "none"),
  nrow = 1.8,
  ncol = 1.8
)

#now same as beta reg , add new columns
subset.best.month.reg <- subset.best.logistic.reg %>%
  mutate(
    months_list = str_split(best.windows, "-", simplify = FALSE) %>%
      map(as.integer),
    n_months = map_int(months_list, length)
  )

########################################
#Machine learning boruta and other check
########################################
#filtering month data before sumution of gene expression
#filtering based on logistic regression
#here I am going back to initial monthly data to sum the expression correctly
#####################################################################
filtered_data.genes.leaves.monthly <- data.genes.leaves.monthly %>%
  inner_join(
    subset.best.month.reg %>%
      select(qseqid, months_list, best.windows, n_months) %>%
      distinct(),
    by = "qseqid"
  ) %>%
  rowwise() %>%
  filter(month %in% months_list) %>%
  ungroup()

#here I will sum the gene expression each genes to the best window
#and combine it to flowering intensity
#####################################################################
filtered_data.genes.leaves.yearly.mastfruiting <- filtered_data.genes.leaves.monthly %>%
  group_by(
    TreeID,
    year,
    qseqid,
    best.windows,
    n_months
  ) %>%
  summarise(
    mean_expr = mean(level.exp.log2, na.rm = TRUE),
    cumsum_expr = sum(level.exp.log2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(average.flo.intensity.individual) #flowering dataset

#Boruta model with all data
###########################
#variable I want to use for boruta predictors
#and then create matrix for boruta
var.of.interest = "cumsum_expr"
gene.masting.wide <- filtered_data.genes.leaves.yearly.mastfruiting %>%
  select(
    TreeID,
    year,
    qseqid,
    {{ var.of.interest }},
    flowering.percentage #
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = qseqid,
    values_from = {{ var.of.interest }}
  )

#predictors and response variable
#for response here, I used directly flowering percentage as it is better to use direct obs
predictors.boruta <- gene.masting.wide %>%
  dplyr::select(-TreeID, -year, -flowering.percentage) # flowering.percentage
response.boruta <- gene.masting.wide$flowering.percentage

#take few minutes
#based on tutorial, and used "general" setup default arguments
set.seed(123)
boruta_result <- Boruta(
  x = predictors.boruta,
  y = response.boruta,
  doTrace = 1,
  pValue = 0.01,
  maxRuns = 2000
)
boruta_result
summary(boruta_result)
plotImpHistory(boruta_result)
#see all first run are shit which is normal when such
#super high number of predictors
#to increase eventually replace F by True for tentative
confirmed_genes.fin <- getSelectedAttributes(
  boruta_result,
  withTentative = T
)
table(boruta_result$finalDecision)

importance.features <- attStats(boruta_result)
importance.features
isoform

#I got the criticism fair that I have way too
#much predictors for small nb of obser (more X than y obs)
##########################################################
run.smaller.ml.boruta = F
if (run.smaller.ml.boruta == T) {
  filtered_data.genes.leaves.monthly.smaller.predictors <- data.genes.leaves.monthly %>%
    inner_join(
      subset.best.month.reg %>%
        filter(AUC > .8) %>%
        dplyr::select(qseqid, months_list, best.windows, n_months) %>%
        distinct(),
      by = "qseqid"
    ) %>%
    rowwise() %>%
    filter(month %in% months_list) %>%
    ungroup()

  unique(filtered_data.genes.leaves.monthly.smaller.predictors$qseqid)
  filtered_data.genes.leaves.yearly.mastfruiting.smaller.predictors <- filtered_data.genes.leaves.monthly.smaller.predictors %>%
    group_by(
      TreeID,
      year,
      qseqid,
      best.windows,
      #categories_gene_exp,
      n_months
    ) %>%
    summarise(
      mean_expr = mean(level.exp.log2, na.rm = TRUE),
      cumsum_expr = sum(level.exp.log2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(average.flo.intensity.individual)

  gene.masting.wide.smaller.predictors <- filtered_data.genes.leaves.yearly.mastfruiting.smaller.predictors %>%
    select(
      TreeID,
      year,
      #best.windows,
      #n_months,
      qseqid,
      {{ var.of.interest }},
      flowering.percentage
    ) %>%
    distinct() %>%
    pivot_wider(
      names_from = qseqid,
      values_from = {{ var.of.interest }}
    )

  predictors.boruta.smaller.predictors <- gene.masting.wide.smaller.predictors %>%
    select(-TreeID, -year, -flowering.percentage)

  response.boruta.smaller.predictors <- gene.masting.wide.smaller.predictors$flowering.percentage
  set.seed(123)
  boruta_result_smaller.predictors <- Boruta(
    x = predictors.boruta.smaller.predictors,
    y = response.boruta.smaller.predictors,
    doTrace = 1,
    pValue = 0.01,
    maxRuns = 2000 #default used 1000
  )

  #plot(boruta_result) #too many features
  plotImpHistory(boruta_result_smaller.predictors) #see all first run are shit which is normal when such
  confirmed_genes.fin.smaller.predictors <- getSelectedAttributes(
    boruta_result_smaller.predictors,
    withTentative = T
  )
  clean_Bor.hvo.smaller.predictors <- reshape_the_Boruta_data(
    boruta_result_smaller.predictors,
    filter.output = T,
    select.features = confirmed_genes.fin.smaller.predictors
  )

  fin.ml.selection.smaller.predictors = clean_Bor.hvo.smaller.predictors %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(value.imp.mean = mean(value)) %>%
    dplyr::select(name, value.imp.mean) %>%
    left_join(
      isoform.cleand %>%
        dplyr::select(name, ProteinName, sseqid)
    ) %>%
    left_join(blastparabidopsis %>% dplyr::rename(sseqid = gene)) %>%
    arrange((value.imp.mean))

  write_csv(
    fin.ml.selection.smaller.predictors,
    "/Users/valentinjourne/Dropbox/F_crenata/temporary_results2025/outputs.genes.importance.smaller.rf.12.2025.csv"
  )

  clean_Bor.hvo.smaller.predictors %>%
    pivot_longer(everything()) %>%
    rename(qseqid = name) %>%
    ggplot(aes(x = fct_reorder(qseqid, value), y = value)) +
    geom_boxplot() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)
    )

  table(confirmed_genes.fin.smaller.predictors %in% confirmed_genes.fin)
  table(confirmed_genes.fin %in% confirmed_genes.fin.smaller.predictors)
  #6 are not present and 33 are similar (two more in the smaller subset predictors)
}
#################################################################################

#Ok after the parenthesis using smaller subset of genes
clean_Bor.hvo <- reshape_the_Boruta_data(
  boruta_result,
  filter.output = T,
  select.features = confirmed_genes.fin
)


fin.ml.selection.last = clean_Bor.hvo %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(value.imp.mean = mean(value)) %>%
  dplyr::select(name, value.imp.mean) %>%
  left_join(
    isoform.cleand %>%
      dplyr::select(name, ProteinName, sseqid)
  ) %>%
  left_join(blastparabidopsis %>% dplyr::rename(sseqid = gene)) %>%
  arrange((value.imp.mean))

#write_csv(
#  fin.ml.selection.last,
#  "/Users/valentinjourne/Dropbox/F_crenata/temporary_results2025/outputs.genes.importance.December2025.csv"
#)
#load the clean one - based on initial test
#fin.ml.selection.meaning = read_csv(
#  "/Users/valentinjourne/Dropbox/F_crenata/temporary_results2025/outputs.genes.importance.December2025.csv"
#)

fin.ml.selection.meaning = fin.ml.selection.last %>%
  left_join(read_table("Tair_gene_name")) %>% #tair id gene name would include specific gene names
  mutate(
    Gene_name = if_else(sseqid == "evm.model.Scaf_2.1813", "APR1", Gene_name)
  ) %>% #add manually APR1 because it is a key gene and I got the annotation from blast)
  group_by(Gene_name) %>%
  mutate(
    n = n(),
    Gene_name = if_else(n == 1, Gene_name, paste0(Gene_name, "_", row_number()))
  ) %>%
  dplyr::select(-n) %>%
  ungroup()

#####################################
#Heatmaps and boruta visualisaiton plot
#####################################
#we removed genes without hits
#because they have no annotation
#we retrieved APR1 using gene sequences

summary_all_genes <- filtered_data.genes.leaves.yearly.mastfruiting %>%
  right_join(
    fin.ml.selection.meaning %>%
      drop_na(sseqid) %>%
      dplyr::select(value.imp.mean, name, Tair.up, Gene_name) %>%
      dplyr::rename(qseqid = name)
  ) %>%
  drop_na(Gene_name) %>%
  group_by(year, Gene_name) %>% #"qseqid"
  summarise(mean.cumsum = mean(cumsum_expr)) %>%
  pivot_wider(names_from = c("year"), values_from = c(mean.cumsum)) %>%
  as.data.frame()

row.names(
  summary_all_genes
) <- summary_all_genes$Gene_name
summary_all_genes$Gene_name <- NULL


data_subset_norm.correlation <- t(apply(
  as.matrix(summary_all_genes),
  1,
  cal_z_score
))

my_hclust_gene.correlation <- hclust(
  dist(data_subset_norm.correlation, method = 'euclidean'),
  method = "ward.D2"
)
#work if smaller subset
fviz_nbclust(
  data_subset_norm.correlation,
  kmeans,
  method = "silhouette",
  k.max = 10
)
nclustertp = 2 #best nb cluster

my_gene_col <- cutree((my_hclust_gene.correlation), k = nclustertp)
my_gene_col = as.data.frame(my_gene_col) %>%
  mutate(cluster = paste0('cluster ', my_gene_col)) %>%
  dplyr::select(cluster)


mycolors = c("black", "#B3B3B3") #colorRampPalette(c(brewer.pal(6, "Set2")))(nclustertp)
cluster = paste0("cluster ", rep(1:nclustertp, each = 1))
my_colourtp = list(
  cluster = setNames(mycolors, cluster)
)

pyearheat.correlation = ComplexHeatmap::pheatmap(
  data_subset_norm.correlation,
  col = colors,
  cluster_rows = my_hclust_gene.correlation,
  cutree_rows = nclustertp,
  cluster_cols = F,
  annotation_colors = my_colourtp,
  annotation_row = my_gene_col,
  row_names_centered = TRUE,
  show_rownames = TRUE,
  #main = "Expression level",
  fontsize_row = 12,
  fontsize_col = 12,
  heatmap_legend_param = list(
    title = "z values",
    title_position = "topcenter",
    legend_direction = "horizontal" # Makes the color bar horizontal
  )
)
#pyearheat.correlation=ComplexHeatmap::draw(pyearheat.correlation, heatmap_legend_side = "bottom")
gg.pyearheat.correlation <- as.ggplot(
  ~ draw(pyearheat.correlation, heatmap_legend_side = "bottom")
)

gg.pyearheat.correlation

cowplot::save_plot(
  'figuresR/heatmaps.genes.final.pdf',
  gg.pyearheat.correlation,
  nrow = 2,
  ncol = 1
)

#plot for boruta
################
ro_1 <- row_order(pyearheat.correlation)
ordered_genes = rownames(data_subset_norm.correlation)[unlist(ro_1)]


boruta.outputs = clean_Bor.hvo %>%
  pivot_longer(everything()) %>%
  dplyr::rename(qseqid = 1) %>%
  right_join(
    fin.ml.selection.meaning %>%
      drop_na(sseqid) %>%
      dplyr::select(value.imp.mean, name, Tair.up, Gene_name) %>%
      dplyr::rename(qseqid = name),
    by = "qseqid" #should have multiple iterations so it is OK
  ) %>%
  mutate(Gene_name = as_factor(Gene_name)) %>%
  mutate(
    Gene_name = factor(Gene_name, levels = rev(ordered_genes))
  ) %>%
  #filter(value != 0) %>% #i removed the 0 because it is noise too many predictors with few observation
  ggplot(aes(x = Gene_name, y = value)) +
  geom_boxplot(
    #outlier.alpha = 0.05,
    outlier.shape = NA,
    fill = colors.roma.short[4],
    col = colors.roma.short[2],
    linewidth = .45,
    alpha = .7
  ) +
  coord_flip() +
  ylab("Variable importance") +
  xlab("") +
  scale_x_discrete(expand = c(0, 0)) +
  theme(
    axis.text.y = element_blank(), #element_text(size = 12),
    plot.margin = margin(5.5, 0, 5.5, 5.5)
  )
boruta.outputs

cowplot::save_plot(
  'figuresR/boruta.outputs.pdf',
  boruta.outputs,
  nrow = 2,
  ncol = .5
)

#boruta.outputs+gg.pyearheat.correlation+ plot_layout(widths = c(1, 2))

#####################################
#Gene ontology GO
#####################################
#now get list of genes for GO analysis
list.genes.cv = cv_sync_toghether %>%
  left_join(
    subset.best.logistic.reg %>%
      filter(term == "cumsum_expr") %>%
      dplyr::select(qseqid, best.windows, AUC, RMSE, AICc)
  ) %>%
  left_join(
    isoform.cleand %>% dplyr::rename(qseqid = name)
  ) %>%
  left_join(blastparabidopsis %>% dplyr::rename(sseqid = gene))

#here it was based on the original data shared
GO.long.format = isoform %>%
  dplyr::select(qseqid, GOs, Subject...38) %>%
  separate_rows(GOs, sep = ",") %>%
  dplyr::rename(GO_ID = GOs, TAIR.ID = Subject...38) %>%
  mutate(
    GO_ID = if_else(GO_ID %in% c("-", "#"), NA, GO_ID),
    TAIR.ID = if_else(TAIR.ID %in% c("No_hits"), NA, TAIR.ID),
    TAIR.ID = gsub("\\..*", "", TAIR.ID)
  )

#ok so I have more TAIR ID than GO ID
#table(is.na(GO.long.format$GO_ID))
#table(is.na(GO.long.format$TAIR.ID))

data.for.plaza = GO.long.format %>%
  dplyr::select(qseqid, TAIR.ID) %>%
  distinct()

GO.data.TAIR = data.for.plaza %>%
  drop_na(TAIR.ID) %>%
  left_join(plaza.go.data)

GO.gene.list = data.genes.leaves.monthly %>%
  dplyr::select(qseqid) %>%
  distinct() %>%
  left_join(GO.data.TAIR)


#now combine my low and high cv to GO genes
combineGO.CVgenes = list.genes.cv %>%
  left_join(GO.data.TAIR)

#number of genes max and mean term
combineGO.CVgenes %>%
  group_by(Category_kCV) %>%
  count(GO_ID) %>%
  summarise(mean_terms = mean(n), max_terms = max(n))


summary_table = GO_analysis(
  list.genes = list.genes.cv %>% filter(AUC > .8), #here is my full list of genes
  GO.gene.list,
  ontology = "BP",
  pvalueCutoff = 0.05
)

data.final.top.gene.GO.function = summary_table %>%
  mutate(adjusted.pvalues = as.numeric(adjusted.pvalues)) %>%
  arrange(raw.pvalues, adjusted.pvalues) %>%
  dplyr::rename(GO_ID = go_id) %>%
  left_join(GO.data.TAIR %>% dplyr::select(GO_ID, FUNCTION) %>% distinct()) %>%
  drop_na(FUNCTION) %>%
  slice_min(order_by = adjusted.pvalues, n = 5, with_ties = F)

colors.pvalues <- colorRampPalette(scico::scico(
  30,
  palette = "vik",
  direction = -1
))(
  100
)

RColorBrewer::brewer.pal(n = 12, name = "PuBu")
#[1] "#FFF7FB" "#ECE7F2" "#D0D1E6" "#A6BDDB" "#74A9CF" "#3690C0" "#0570B0" "#045A8D" "#023858"

colors.pvalues <- colorRampPalette(c(
  "#C6DBEF",
  "#9ECAE1",
  "#6BAED6",
  "#4292C6",
  "#2171B5",
  "#08519C",
  "#08306B"
))(100)

top.5.GO.analysis = data.final.top.gene.GO.function %>%
  mutate(function.go = paste0(FUNCTION, "\n", GO_ID)) %>%
  dplyr::rename(`adj p-val` = adjusted.pvalues) %>%
  ggplot(aes(y = function.go, x = ratio)) +
  #geom_col(fill = "black", width = 0.001) +
  geom_point(aes(col = `adj p-val`, size = ratio)) +
  geom_text(aes(label = universeCounts), vjust = -0.9, color = "black") +
  ylab("") +
  #scale_x_continuous(labels = scales::percent) +
  xlab("Ratio") +
  scale_colour_continuous(palette = colors.pvalues) +
  scale_fill_continuous(palette = colors.pvalues) +
  theme(legend.position = c(.75, .75)) +
  scale_size_continuous(range = c(3, 10)) +
  guides(
    color = guide_colorbar(title = "Adj. \np-value"),
    size = "none"
  ) +
  theme(
    #legend.position = "right",
    legend.box = "vertical",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )

top.5.GO.analysis
cowplot::save_plot(
  'figuresR/Figure.GOanalysis.pdf',
  top.5.GO.analysis,
  nrow = 1,
  ncol = 1.35
)

###################################
#Supplementary figures for the relationship between CV and synchrony
###################################
cv_sync_extremes <- list.genes.cv %>%
  dplyr::select(-sseqid) %>% #to avoid trouble
  filter(AUC > .8) %>%
  left_join(
    fin.ml.selection.meaning %>%
      dplyr::select(name, Gene_name, sseqid) %>%
      dplyr::rename(qseqid = name)
  )

plot1.auc.sync = ggplot(
  cv_sync_extremes,
  aes(x = AUC, y = mean_synchrony)
) +
  geom_point(aes(fill = kCV, col = kCV), shape = 21, stroke = 0.2, size = 3) +
  geom_text_repel(
    data = cv_sync_extremes %>% drop_na(Gene_name, sseqid),
    aes(label = Gene_name),
    size = 3,
    fontface = "bold",
    max.overlaps = 50
  ) +
  ylab("Synchrony") +
  xlab("Area Under the ROC Curve (AUC)") +
  scale_fill_viridis_c(name = "kCV", option = "magma") +
  scale_color_viridis_c(name = "kCV", option = "magma") +
  theme(legend.position = "right")

plot2.auc.sync = ggplot(
  cv_sync_extremes,
  aes(x = AUC, y = kCV)
) +
  ylab("kCV") +
  xlab("Area Under the ROC Curve (AUC)") +
  geom_point(
    aes(fill = mean_synchrony, col = mean_synchrony),
    shape = 21,
    stroke = 0.2,
    size = 3
  ) +
  geom_text_repel(
    data = cv_sync_extremes %>% drop_na(Gene_name, sseqid),
    aes(label = Gene_name),
    size = 3,
    fontface = "bold",
    max.overlaps = 50
  ) +
  scale_fill_viridis_c(option = "plasma", name = "Sync.") +
  scale_color_viridis_c(option = "plasma", name = "Sync.") +
  theme(legend.position = "right")

Figure3d = plot1.auc.sync +
  plot2.auc.sync +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "bottom")
Figure3d
cowplot::save_plot(
  'figuresR/Figure.sync.kcv.AUC.pdf',
  Figure3d,
  dpi = 300,
  nrow = 1.4,
  ncol = 1.6
)

##################################################################
#Seasonal gene expression of the main genes identified by boruta
##################################################################

#all boruta genes
#################
genes.temporal.filtering = fin.ml.selection.meaning %>%
  dplyr::select(name, Gene_name, sseqid) %>%
  dplyr::rename(qseqid = name) %>%
  left_join(data.genes.leaves.monthly) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  )

#highlight geom ribbon
######################
highlight_2014 <- data.frame(
  xmin = "6.2014",
  xmax = "9.2014",
  ymin = -Inf,
  ymax = Inf
)

highlight_2017 <- data.frame(
  xmin = "6.2017",
  xmax = "9.2017",
  ymin = -Inf,
  ymax = Inf
)
highlight_2021 <- data.frame(
  xmin = "6.2021",
  xmax = "9.2021",
  ymin = -Inf,
  ymax = Inf
)

#for labelings
#seelct one genes (all have the same )
df.labels = genes.temporal.filtering %>% dplyr::filter(Gene_name == "FLC")
x_labels <- ifelse(
  df.labels$month == 6,
  paste0(df.labels$month, "-", df.labels$year),
  #paste0(df$month, "\n", df$year),
  paste0("", "\n")
  #as.character(df1$month)
)
x_levels <- with(df.labels, interaction(month, year, drop = TRUE))

label_df_oas.flo <- genes.temporal.filtering %>%
  drop_na(sseqid) %>%
  group_by(Gene_name) %>%
  filter(year == max(year)) %>%
  slice_tail(n = 1) %>% # safety in case of ties
  ungroup()

all.genes.figures.exp.time = ggplot(
  genes.temporal.filtering %>%
    drop_na(sseqid),
  aes(
    x = interaction(month, year),
    y = level.exp,
    group = TreeID
  )
) +
  facet_wrap(~Gene_name, scales = "free_y", ncol = 4) +
  geom_line(alpha = 0.8) +
  scale_x_discrete(
    limits = x_levels,
    labels = x_labels
  ) +
  theme_minimal(base_size = 12) +
  ggpubr::theme_pubr() +
  xlab("") +
  ylab("Gene expression (tpm)") +
  #hrbrthemes::theme_ipsum(base_size = 14, axis_title_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  geom_rect(
    data = highlight_2014,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey30",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2017,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey30",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2021,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey30",
    alpha = 0.2
  ) +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  )
all.genes.figures.exp.time

cowplot::save_plot(
  'figuresR/All.Boruta.Genes.pdf',
  all.genes.figures.exp.time,
  nrow = 3,
  ncol = 3
)

#now general expression over time
#focus on flowering genes and OAS cluster genes
###############################################
genes.temporal.filtering.oas.flo.sub = genes.temporal.filtering %>%
  drop_na(Gene_name, sseqid) %>%
  mutate(
    group.genes = case_when(
      Gene_name %in% c("FLC", "FT") ~ "Flowering genes",
      str_detect(Gene_name, "APR|PQ|LSU|SDI|TAA|SHM7") ~ "OAS cluster genes",
      TRUE ~ "Others"
    )
  )

genes.temporal.filtering.oas.flo.sub1 = genes.temporal.filtering.oas.flo.sub %>%
  filter(group.genes %in% c("Flowering genes", "OAS cluster genes")) %>%
  group_by(group.genes, Gene_name, month, year) %>%
  #summarise(
  #  mean.exp = mean(level.exp, na.rm = TRUE),
  #  sd.exp = sd(level.exp)
  #) %>%
  mutate(
    group.genes = factor(
      #relevel factor
      group.genes,
      levels = c("Flowering genes", "OAS cluster genes", "Others")
    )
  ) %>%
  mutate(
    Gene_name = factor(
      Gene_name,
      levels = c(
        "FT",
        "APR1",
        "APR_3",
        "LSU2",
        "SDI1_1",
        "FLC",
        "SDI2",
        "SHM7/MSA1",
        "TAA_1",
        "PQ",
        setdiff(
          unique(Gene_name),
          c(
            "FT",
            "APR1",
            "APR_3",
            "LSU2",
            "SDI1_1",
            "FLC",
            "SDI2",
            "SHM7/MSA1",
            "TAA_1",
            "PQ"
          )
        )
      )
    )
  ) %>%
  dplyr::filter(
    Gene_name %in%
      c(
        "FT",
        "APR1",
        "APR_3",
        "LSU2",
        "SDI1_1",
        "FLC",
        "SDI2",
        "SHM7/MSA1",
        "TAA_1",
        "PQ"
      )
  )

label_df_oas.flo.all <- genes.temporal.filtering.oas.flo.sub1 %>%
  group_by(Gene_name) %>%
  filter(year == max(year)) %>%
  slice_tail(n = 1) %>% # safety in case of ties
  ungroup()

oas.flo.exp = ggplot(
  genes.temporal.filtering.oas.flo.sub1,
  aes(
    x = interaction(month, year),
    y = level.exp,
    group = TreeID,
    col = TreeID
  )
) +
  geom_rect(
    data = highlight_2014,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey30",
    alpha = 0.1
  ) +
  geom_rect(
    data = highlight_2017,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey30",
    alpha = 0.1
  ) +
  geom_rect(
    data = highlight_2021,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey30",
    alpha = 0.1
  ) +
  geom_line() +
  facet_wrap(~Gene_name, scales = "free_y", nrow = 2) +
  theme_minimal() +
  xlab("") +
  ylab("Gene expression (tpm)") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  ) +
  scale_x_discrete(
    limits = x_levels,
    labels = x_labels
  ) +
  scale_color_manual(values = TreeID.color.pal) +
  scale_fill_manual(values = TreeID.color.pal)

oas.flo.exp
cowplot::save_plot(
  'figuresR/Figure.timesseries.OASFlo.pdf',
  oas.flo.exp,
  nrow = 1,
  ncol = 2
)
# #second option
# data.new.scaled.oas.flo = genes.temporal.filtering.oas.flo.sub %>%
#   filter(group.genes %in% c("Flowering genes", "OAS cluster genes")) %>%
#   group_by(group.genes, Gene_name) %>%
#   mutate(scale.expression = scale(level.exp)) %>%
#   ungroup() %>%
#   group_by(group.genes, Gene_name, month, year) %>%
#   summarise(
#     mean.exp = mean(scale.expression, na.rm = TRUE),
#     sd.exp = sd(scale.expression)
#   ) %>%
#   mutate(
#     group.genes = factor(
#       group.genes,
#       levels = c("Flowering genes", "OAS cluster genes", "Others")
#     )
#   )
#
# label_df_oas.flo.scaled <- data.new.scaled.oas.flo %>%
#   group_by(Gene_name, group.genes) %>%
#   filter(year == max(year)) %>%
#   slice_tail(n = 1) %>% # safety in case of ties
#   ungroup()
#
# oas.flo.plot.timeseries.std = ggplot(
#   data.new.scaled.oas.flo,
#   aes(
#     x = interaction(month, year),
#     y = mean.exp,
#     group = Gene_name,
#     col = Gene_name
#   )
# ) +
#   facet_grid(. ~ group.genes) +
#   geom_text_repel(
#     data = label_df_oas.flo.scaled,
#     aes(label = Gene_name),
#     hjust = 0,
#     direction = "y",
#     xlim = c(NA, Inf),
#     nudge_x = 3,
#     segment.color = "black",
#     size = 3,
#     show.legend = FALSE,
#     max.overlaps = Inf
#   ) +
#   geom_rect(
#     data = highlight_2014,
#     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#     inherit.aes = FALSE,
#     fill = "grey40",
#     alpha = 0.2
#   ) +
#   geom_rect(
#     data = highlight_2017,
#     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#     inherit.aes = FALSE,
#     fill = "grey40",
#     alpha = 0.2
#   ) +
#   geom_rect(
#     data = highlight_2021,
#     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#     inherit.aes = FALSE,
#     fill = "grey40",
#     alpha = 0.2
#   ) +
#   geom_line(linewidth = .6) +
#   theme(
#     axis.text.x = element_text(vjust = 0.5, size = 12, angle = 45),
#     axis.text.y = element_text(size = 12),
#     axis.title = element_text(size = 14),
#     strip.text = element_text(face = "bold", size = 12),
#     panel.grid.minor = element_blank(),
#     legend.position = "none",
#     strip.background = element_blank(),
#     plot.margin = margin(5.5, 60, 5.5, 5.5)
#   ) +
#   scale_x_discrete(
#     limits = x_levels,
#     labels = x_labels
#   ) +
#   xlab("Year") +
#   ylab("Gene expression level (dimensionless)") +
#   scale_color_manual(values = rep(genes.color.pal, length.out = 13)) +
#   scale_fill_manual(values = rep(genes.color.pal, length.out = 13)) +
#   coord_cartesian(clip = "off")
# oas.flo.plot.timeseries.std
# cowplot::save_plot(
#   'figuresR/Figure.timesseries.OASFlo.pdf',
#   oas.flo.plot.timeseries.std,
#   nrow = 1.1,
#   ncol = 1.4
# )

#PLOT FOR FLC AND FT

FLC.plot = ggplot(
  genes.temporal.filtering.oas.flo.sub1 %>% filter(Gene_name == "FLC"),
  aes(x = interaction(month, year), y = level.exp, group = TreeID, col = TreeID)
) +
  geom_line() +
  #scale_x_discrete(
  # labels = unique(interaction(df$month, factor(df$year), sep = "\n"))
  #) +
  scale_x_discrete(
    limits = x_levels,
    labels = x_labels
  ) +
  theme_bw() +
  xlab("") +
  ylab("Gene expression (tpm)") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  geom_rect(
    data = highlight_2014,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2017,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2021,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  scale_color_manual(values = TreeID.color.pal) +
  scale_fill_manual(values = TreeID.color.pal)
FLC.plot

FT.plot = ggplot(
  genes.temporal.filtering.oas.flo.sub1 %>% filter(Gene_name == "FT"),
  aes(x = interaction(month, year), y = level.exp, group = TreeID, col = TreeID)
) +
  geom_line() +
  #scale_x_discrete(
  #  labels = unique(interaction(df1$month, factor(df1$year), sep = "\n"))
  #) +
  scale_x_discrete(
    limits = x_levels,
    labels = x_labels
  ) +
  theme_bw() +
  xlab("") +
  ylab("Gene expression (tpm)") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  geom_rect(
    data = highlight_2014,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2017,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2021,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  scale_color_manual(values = TreeID.color.pal) +
  scale_fill_manual(values = TreeID.color.pal)
FT.plot
expressionFTFLC = (FT.plot) +
  (FLC.plot + ylab("")) +
  #plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 2)
expressionFTFLC
cowplot::save_plot(
  'figuresR/FigureFLCFT.pdf',
  expressionFTFLC,
  nrow = .8,
  ncol = 1.5
)

#additional time series for other important genes
#Bark Storage Proteins (BSPs) genes
###################################

bsp.list <- readxl::read_excel(
  "/Users/valentinjourne/Dropbox/F_crenata/ISO-seq/BSPcandidate.xlsx",
  sheet = 1,
  col_names = F
) %>%
  dplyr::rename(qseqid = 1)

genes.temporal.bsp = bsp.list %>%
  dplyr::select(qseqid, 2, 3) %>%
  left_join(data.genes.leaves.monthly) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  )

BSP.plot = ggplot(
  genes.temporal.bsp %>% drop_na(),
  aes(
    x = interaction(month, year),
    y = level.exp,
    group = TreeID,
    col = TreeID
  )
) +
  scale_x_discrete(
    limits = x_levels,
    labels = x_labels
  ) +
  geom_line() +
  facet_wrap(~qseqid, scales = "free_y", nrow = 1) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  ) +
  geom_rect(
    data = highlight_2014,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2017,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2021,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  xlab("") +
  ylab("Gene expression (tpm)") +
  scale_color_manual(values = TreeID.color.pal) +
  scale_fill_manual(values = TreeID.color.pal)
BSP.plot
cowplot::save_plot(
  'figuresR/FigureBSP.pdf',
  BSP.plot,
  nrow = .8,
  ncol = 1.5
)

#mainenant for slim 1

slim1.list = isoform.cleand %>%
  dplyr::rename(qseqid = 1) %>%
  left_join(blastparabidopsis %>% dplyr::rename(sseqid = gene)) %>%
  filter(Tair.up == "AT1G73730") %>% #tair ID for slim 1 or EIL3 is the same
  dplyr::select(qseqid, Tair.up) %>%
  left_join(data.genes.leaves.monthly) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  )

unique(slim1.list$qseqid)

slim1.plot = ggplot(
  slim1.list %>% drop_na(),
  aes(
    x = interaction(month, year),
    y = level.exp,
    group = TreeID,
    col = TreeID
  )
) +
  ggpubr::theme_pubr() +
  scale_x_discrete(
    limits = x_levels,
    labels = x_labels
  ) +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  ) +
  geom_rect(
    data = highlight_2014,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2017,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2021,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_line() +
  xlab("") +
  ylab("Gene expression (tpm)") +
  scale_color_manual(values = TreeID.color.pal) +
  scale_fill_manual(values = TreeID.color.pal)
slim1.plot
cowplot::save_plot(
  'figuresR/FigureSLIM1.pdf',
  slim1.plot,
  nrow = .8,
  ncol = .6
)


#FLC the other version of the genes 2
sub.flc.2version = data.genes.all %>%
  filter(Tissue == "L") %>%
  dplyr::rename(name = qseqid) %>%
  right_join(
    isoform.cleand %>%
      dplyr::select(name, sseqid) %>%
      filter(
        name %in%
          c(
            #"Facr_v2.5_s1cl064033",
            #"Facr_v2.5_s1cl064034",
            #"Facr_v2.5_s1cl067823",
            #"Facr_v2.5_s1cl067824"
            "Facr_v2.5_s1cl076318"
          )
      )
  ) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  ) %>%
  drop_na(IndexNo) #remove the missing ones

plot.FLC2 = ggplot(
  sub.flc.2version,
  aes(
    x = interaction(month, year),
    y = level.exp,
    group = TreeID,
    col = TreeID
  )
) +
  ggpubr::theme_pubr() +
  scale_x_discrete(
    limits = x_levels,
    labels = x_labels
  ) +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  ) +
  geom_rect(
    data = highlight_2014,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2017,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  geom_rect(
    data = highlight_2021,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey",
    alpha = 0.2
  ) +
  xlab("") +
  ylab("Gene expression (tpm)") +
  scale_color_manual(values = TreeID.color.pal) +
  scale_fill_manual(values = TreeID.color.pal) +
  geom_line()
plot.FLC2
cowplot::save_plot(
  'figuresR/FigureFLC2.pdf',
  plot.FLC2,
  nrow = .8,
  ncol = .6
)

######################################
genes.satake <- c(
  "AT5G48850",
  "AT1G04770",
  "AT4G08620",
  "AT1G78000",
  "AT5G10180",
  "AT5G19600",
  "AT5G13550",
  "AT3G12520",
  "AT5G26220",
  "AT5G18170",
  "AT5G07440",
  "AT5G08100",
  "AT4G17830"
)

subset.sulfur.exp = isoform.cleand %>%
  dplyr::rename(qseqid = 1) %>%
  left_join(blastparabidopsis %>% dplyr::rename(sseqid = gene)) %>%
  filter(Tair.up %in% c(genes.satake)) %>% #tair ID for slim 1 or EIL3 is the same
  dplyr::select(qseqid, Tair.up) %>%
  left_join(data.genes.leaves.monthly) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  ) %>%
  group_by(qseqid, Tair.up, month, year) %>%
  summarise(
    mean.exp = mean(level.exp, na.rm = TRUE),
    sd.exp = sd(level.exp)
  ) %>%
  drop_na(month)


ggplot(
  subset.sulfur.exp,
  aes(x = interaction(month, year), y = mean.exp, group = qseqid)
) +
  geom_line() +
  facet_wrap(. ~ Tair.up, ncol = 6, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))


###################################
#Amino acid relationship
###################################
#amino acid initial weird table format
AA.Phloem <- readxl::read_excel(
  "/Users/valentinjourne/Dropbox/F_crenata/data_amino_acid/Amino acid data_Fagus (2021-2023)_Miyazawa.xlsx",
  sheet = 4,
  col_names = FALSE
)
amino.acid.phloem = format.amino.acid.long(AA.Phloem) %>%
  mutate(Tissue = "PHLOEM") %>%
  mutate(Sample = str_remove(Sample, 'P'))

AA.Xylem <- readxl::read_excel(
  "/Users/valentinjourne/Dropbox/F_crenata/data_amino_acid/Amino acid data_Fagus (2021-2023)_Miyazawa.xlsx",
  sheet = 5,
  col_names = FALSE
)

amino.acid.xylem = format.amino.acid.long(AA.Xylem) %>%
  mutate(Tissue = "XYLEM") %>%
  mutate(Sample = str_remove(Sample, 'X'))

# aggregate.aa.tissue = bind_rows(amino.acid.phloem, amino.acid.xylem) %>%
#   mutate(Year = lubridate::year(Date)) %>%
#   filter(Year < 2023) %>%
#   group_by(Date, Year, Month, Compound, Tissue) %>%
#   summarise(
#     mean.value = mean(Value, na.rm = T),
#     median.value = median(Value, na.rm = T),
#     sd = sd(Value, na.rm = T)
#   ) %>%
#   ungroup() %>%
#   mutate(Year.chr = as.character(Year))

data.cleaned.aa.phloem = amino.acid.phloem %>%
  mutate(Year = lubridate::year(Date)) %>%
  filter(Compound %in% c("NH4+", "NO3-")) %>%
  mutate(Year.chr = as.character(Year)) %>%
  filter(Year < 2023)


#Analaysis with sulfur data, GSH
#################################
total.sufur <- readxl::read_excel(
  "/Users/valentinjourne/Dropbox/F_crenata/Maruyama/FW修正版/sulfur.xlsx",
  sheet = 1,
  col_names = T
)

#test on phloem data
####################
total.sufur.gsh.clean.phloem = total.sufur %>%
  dplyr::select(
    `Sample name`,
    `sulfate_p\r\n (μmol/gFW)`,
    `GSH_p (nmol/gFW)`,
    `totalS_p (mg/gFW)`
  ) %>%
  separate(
    col = `Sample name`,
    into = c("Sample", "Year"),
    sep = "([A-Z])_|_",
    remove = FALSE
  ) %>%
  mutate(
    Sample = as.integer(Sample),
    Year = as.integer(Year)
  ) %>%
  dplyr::rename(sulfate = 4, GSH = 5, totalS = 6) %>%
  mutate(Tissue = "PHLOEM") %>% #here my tissue I want
  filter(Year < 2023) %>% #two years only
  pivot_longer(
    cols = c(sulfate, GSH, totalS),
    names_to = "compound",
    values_to = "concentration"
  ) %>%
  mutate(
    Year = as.factor(Year),
    compound = factor(
      compound,
      levels = c("sulfate", "GSH", "totalS"),
      labels = c("Sulfate", "GSH", "Total S")
    )
  )

#function to do wilcoxon test stat and get the sig label for the plot
data.for.phloem.gsh.totalS = do.stat.element_gsh_sulfate_totalS(
  total.sufur.gsh.clean.phloem
)


phloem.gsh.sulfate = ggplot(
  data.for.phloem.gsh.totalS$summary_stats,
  aes(x = Year, y = mean, fill = Year)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    linewidth = 0.4,
    alpha = .7
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    linewidth = 0.7
  ) +
  geom_jitter(
    data = total.sufur.gsh.clean.phloem,
    aes(x = Year, y = concentration),
    width = 0.08,
    size = 2,
    alpha = 0.6,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = data.for.phloem.gsh.totalS$sig_labels,
    aes(x = 1.5, y = y_pos, label = sig_label),
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  facet_wrap(~compound, scales = "free_y") +
  #scale_fill_manual(values = c("2021" = "#5B8DB8", "2022" = "#E07B54")) + #c("#008B45FF", "#EE0000FF")
  scale_fill_manual(values = c("2021" = "#008B45FF", "2022" = "#EE0000FF")) + #c("#008B45FF", "#EE0000FF")
  labs(
    x = "Year",
    y = "Concentration (µmol/g FW)",
    fill = "Year"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none"
  )
phloem.gsh.sulfate
cowplot::save_plot(
  'figuresR/phloem.gsh.sulfatephloem.pdf',
  phloem.gsh.sulfate,
  nrow = .8,
  ncol = 1
)

#now the same with Xylem
########################
total.sufur.gsh.clean.xylen = total.sufur %>%
  dplyr::select(
    `Sample name`,
    `sulfate_x\r\n (mmol/L)`,
    `GSH_x（µM)`
  ) %>%
  separate(
    col = `Sample name`,
    into = c("Sample", "Year"),
    sep = "([A-Z])_|_",
    remove = FALSE
  ) %>%
  mutate(
    Sample = as.integer(Sample),
    Year = as.integer(Year)
  ) %>%
  dplyr::rename(sulfate = 4, GSH = 5) %>%
  mutate(Tissue = "XYLEM") %>% #here my tissue I want
  filter(Year < 2023) %>% #two years only
  pivot_longer(
    cols = c(sulfate, GSH),
    names_to = "compound",
    values_to = "concentration"
  ) %>%
  mutate(
    Year = as.factor(Year),
    compound = factor(
      compound,
      levels = c("sulfate", "GSH"),
      labels = c("Sulfate", "GSH")
    )
  )

data.for.xylem.gsh.totalS = do.stat.element_gsh_sulfate_totalS(
  total.sufur.gsh.clean.xylen
)

xylem.gsh.sulfate = ggplot(
  data.for.xylem.gsh.totalS$summary_stats,
  aes(x = Year, y = mean, fill = Year)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    linewidth = 0.4,
    alpha = .7
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    linewidth = 0.7
  ) +
  geom_jitter(
    data = total.sufur.gsh.clean.xylen,
    aes(x = Year, y = concentration),
    width = 0.08,
    size = 2,
    alpha = 0.6,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = data.for.xylem.gsh.totalS$sig_labels,
    aes(x = 1.5, y = y_pos, label = sig_label),
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  facet_wrap(~compound, scales = "free_y") +
  scale_fill_manual(values = c("2021" = "#008B45FF", "2022" = "#EE0000FF")) + #c("#008B45FF", "#EE0000FF")
  labs(
    x = "Year",
    y = "Concentration",
    fill = "Year"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none"
  )
xylem.gsh.sulfate
cowplot::save_plot(
  'figuresR/xylem.gsh.sulfatephloem.pdf',
  xylem.gsh.sulfate,
  nrow = .8,
  ncol = 1 / 2
)

#additional test now again with aa for NO3 and NH4
##################################################
data.t.test.N = bind_rows(amino.acid.phloem, amino.acid.xylem) %>%
  mutate(Year = lubridate::year(Date)) %>%
  filter(Year < 2023) %>%
  filter(Compound %in% c("NH4+", "NO3-")) %>%
  mutate(MastONOFF = if_else(Year == 2021, "Mast", "Non-Mast")) %>%
  mutate(
    Year = as.factor(Year),
    Sample = as.factor(Sample),
    Month = as.factor(Month)
  ) %>%
  mutate(
    log.value = log1p(Value),
    Value.noise = if_else(Value == 0, 0.001, Value)
  )


sub.data.t.test.NH4 = data.t.test.N %>%
  filter(Compound == "NH4+" & Tissue == "XYLEM")

#model.NH4.XYLEM <- glmmTMB(
#  Value ~ Year + (1 | Sample) + (1 | Month),
#  family = tweedie(link = "log"),
#  data = sub.data.t.test.NH4
#)
#summary(model.NH4.XYLEM)
#hist((sub.data.t.test.NH4$Value))
#twwedie too complex, but trends are the same

model.NH4.XYLEM <- glmmTMB(
  log.value ~ Year + (1 | Sample) + (1 | Month), #I tried with fixed month, results are same same
  data = sub.data.t.test.NH4
)
summary(model.NH4.XYLEM)
emmeans(model.NH4.XYLEM, ~Year)
pairs(emmeans(model.NH4.XYLEM, ~Year))
#simulationOutput <- simulateResiduals(fittedModel = model.NH4.XYLEM, plot = F)
#plot(simulationOutput)
#shapiro.test(sub.data.t.test.N$Value[sub.data.t.test.N$MastONOFF == "Non-Mast"])
#wilcox.test(Value ~ MastONOFF, data = sub.data.t.test.N)

sub.data.t.test.NH4.ph = data.t.test.N %>%
  filter(Compound == "NH4+" & Tissue != "XYLEM")
model.NH4.PHLOEM <- glmmTMB(
  log.value ~ Year + (1 | Sample) + (1 | Month), #I tried with fixed month, results are same same
  data = sub.data.t.test.NH4.ph
)
summary(model.NH4.PHLOEM)

sub.data.t.test.NO3.ph = data.t.test.N %>%
  filter(Compound == "NO3-" & Tissue != "XYLEM")
model.NO3.PHLOEM <- glmmTMB(
  log.value ~ Year + (1 | Sample) + (1 | Month), #I tried with fixed month, results are same same
  data = sub.data.t.test.NO3.ph
)
summary(model.NO3.PHLOEM)

sub.data.t.test.NO3.x = data.t.test.N %>%
  filter(Compound == "NO3-" & Tissue == "XYLEM")
model.NO3.XYLEM <- glmmTMB(
  log.value ~ Year + (1 | Sample) + (1 | Month), #I tried with fixed month, results are same same
  data = sub.data.t.test.NO3.x
)
summary(model.NO3.XYLEM)

#make a plot for NH4
sub.data.t.test.NH4.summary.stat = sub.data.t.test.NH4.ph %>%
  group_by(Year) %>%
  summarise(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    se = sd / sqrt(n),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

nh4.ph = ggplot(
  sub.data.t.test.NH4.summary.stat,
  aes(x = Year, y = mean, fill = Year)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    linewidth = 0.4,
    alpha = .7
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    linewidth = 0.7
  ) +
  scale_fill_manual(values = c("2021" = "#008B45FF", "2022" = "#EE0000FF")) + #c("#008B45FF", "#EE0000FF")
  labs(
    x = "Year",
    y = "Concentration NH4",
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none"
  )
nh4.ph
cowplot::save_plot(
  'figuresR/nh4.ph.pdf',
  nh4.ph,
  nrow = .8,
  ncol = 1 / 2
)

#make a plot for NO3
sub.data.t.test.N03.summary.stat = sub.data.t.test.NO3.ph %>%
  group_by(Year) %>%
  summarise(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    se = sd / sqrt(n),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

no3 = ggplot(
  sub.data.t.test.N03.summary.stat,
  aes(x = Year, y = mean, fill = Year)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    linewidth = 0.4,
    alpha = .7
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    linewidth = 0.7
  ) +
  scale_fill_manual(values = c("2021" = "#008B45FF", "2022" = "#EE0000FF")) + #c("#008B45FF", "#EE0000FF")
  labs(
    x = "Year",
    y = "Concentration NO3",
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none"
  )
no3
cowplot::save_plot(
  'figuresR/no3.ph.pdf',
  no3,
  nrow = .8,
  ncol = 1 / 2
)


#now only single for total S
totalS = ggplot(
  data.for.phloem.gsh.totalS$summary_stats %>% filter(compound == "Total S"),
  aes(x = Year, y = mean, fill = Year)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    linewidth = 0.4,
    alpha = .7
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    linewidth = 0.7
  ) +
  #scale_fill_manual(values = c("2021" = "#5B8DB8", "2022" = "#E07B54")) + #c("#008B45FF", "#EE0000FF")
  scale_fill_manual(values = c("2021" = "#008B45FF", "2022" = "#EE0000FF")) + #c("#008B45FF", "#EE0000FF")
  labs(
    x = "Year",
    y = "Concentration Total S",
    fill = "Year"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none"
  )
totalS
cowplot::save_plot(
  'figuresR/totalS.pdf',
  totalS,
  nrow = .8,
  ncol = 1 / 2
)
################################################################################################################################
################################################################################################################
################################################################################################################
########################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################################################
########################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################################################

final.highcv.tair.function.go = list.genes.cv %>%
  left_join(GO.gene.list) %>%
  left_join(
    isoform.cleand %>%
      dplyr::select(name, sseqid, ProteinName) %>%
      dplyr::rename(qseqid = name)
  ) %>%
  right_join(summary_table %>% dplyr::rename(GO_ID = go_id)) %>%
  distinct()


cc = final.highcv.tair.function.go %>% distinct()

unique(final.highcv.tair.function.go$sseqid)
# write_csv(
#   final.highcv.tair.function.go %>% dplyr::select(-combined.metric),
#   here::here("final.highkcv.tair.function.go.csv")
# )

#evm.model.Scaf_5.31 #FT
#evm.model.Scaf_3.2022 #FLC

FLC = data.genes.all %>%
  left_join(
    isoform.cleand %>%
      dplyr::select(name, sseqid) %>%
      distinct() %>%
      dplyr::rename(qseqid = name)
  ) %>%
  left_join(
    blastparabidopsis %>%
      dplyr::select(Tair.up, gene) %>%
      dplyr::rename(sseqid = gene)
  ) %>%
  #filter(Tair.up == "AT5G10140") %>%
  filter(qseqid == "Facr_v2.5_s1cl026553") %>%
  filter(month != 10) %>%
  filter(Tissue == "L")

ggplot(
  FLC,
  aes(x = interaction(month, year), y = level.exp, group = TreeID)
) +
  geom_line()

subsrt.immuno.genes = final.highcv.tair.function.go %>%
  filter(
    TAIR.ID %in%
      c(
        "AT3G23240",
        "AT2G47730",
        "AT3G25830",
        "AT1G74590",
        "AT4G16740",
        "AT5G23960"
      )
  ) %>%
  dplyr::select(qseqid, TAIR.ID, CV) %>%
  distinct() %>%
  left_join(data.genes.leaves.monthly) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  )


subsrt.immuno.genes %>%
  mutate(newgene.name = paste0(TAIR.ID, "; CV=", round(CV, 3))) %>%
  ggplot(aes(
    x = timepoint,
    y = (level.exp),
    group = TreeID,
    col = TreeID
  )) +
  facet_wrap(~newgene.name, scales = "free_y") +
  geom_line(alpha = 0.6) +
  theme_minimal(base_size = 12) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 10, angle = 90),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Gene expression")


##############################################################################################################
#########################################################################################################################
####################################################################################################################################
#########################################################################################################################
####################################################################################################################################
##############################################################################################################
#########################################################################################################################
####################################################################################################################################
#########################################################################################################################

meaningGenBVOC = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/RNAseqdata/F.crenata_Flowering_BVOC_genes_summary_geneID.xlsx',
  sheet = 4
)
meaningGenBVOC.function = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/RNAseqdata/F.crenata_Flowering_BVOC_genes_summary_geneID.xlsx',
  sheet = 3
) %>%
  dplyr::select(Gene_name, At_gene_id, F.crenata, Pathway_name) %>%
  drop_na(Gene_name) %>%
  dplyr::rename(
    Gene_name_Original = Gene_name,
    AT_gene_id = At_gene_id
  )

meaningFloGenLFTC.function %>% dplyr::select(1, 2)

flowering.genes.variation = data.genes.leaves.monthly %>%
  right_join(
    isoform.cleand %>%
      dplyr::select(name, sseqid) %>%
      dplyr::rename(qseqid = name)
  ) %>%
  left_join(GO.gene.list %>% dplyr::select(qseqid, TAIR.ID) %>% distinct()) %>%
  right_join(
    meaningFloGenLFTC.function %>% dplyr::rename(TAIR.ID = AT_gene_id)
  ) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  ) %>%
  drop_na(qseqid) %>%
  left_join(
    cv_sync_toghether %>%
      dplyr::select(qseqid, kCV, mean_synchrony, Category_kCV, Category_sync)
  ) %>%
  mutate(
    label = paste0(
      "kCV = ",
      round(kCV, 2),
      "\n",
      "Sync = ",
      round(mean_synchrony, 2)
    )
  )

flowering.genes.variation.selection = flowering.genes.variation %>%
  filter(Category_sync == "High_sync") #High_sync


schemes.to.keep = "Photoperiod" #"Flo|flo"
label_df <- flowering.genes.variation.selection %>%
  #filter(str_detect(schemes, schemes.to.keep)) %>%
  group_by(Gene_name_Original) %>%
  summarise(
    kCV = unique(kCV),
    mean_synchrony = unique(mean_synchrony),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "kCV = ",
      round(kCV, 2),
      "\n",
      "Sync = ",
      round(mean_synchrony, 2)
    )
  )

flowering.genes.variation.selection %>%
  #filter(str_detect(schemes, schemes.to.keep)) %>%
  ggplot(aes(
    x = timepoint,
    y = level.exp,
    group = TreeID
  )) +
  facet_wrap(~Gene_name_Original, scales = "free") +
  geom_line(alpha = 0.6) +
  theme_minimal(base_size = 12) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  geom_vline(
    xintercept = which(levels(genes.temporal.filtering$timepoint) == "6-2014"),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = which(levels(genes.temporal.filtering$timepoint) == "6-2017"),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = which(levels(genes.temporal.filtering$timepoint) == "6-2021"),
    linetype = "dashed",
    color = "red"
  ) +
  geom_text(
    data = label_df,
    aes(label = label),
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.1,
    color = "red",
    size = 3,
    inherit.aes = FALSE
  ) +
  xlab("") +
  ylab("Gene expression")

meaningGenBVOC = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/RNAseqdata/F.crenata_Flowering_BVOC_genes_summary_geneID.xlsx',
  sheet = 4
)
meaningGenBVOC.function = read_excel(
  '/Users/valentinjourne/Dropbox/F_crenata/RNAseqdata/F.crenata_Flowering_BVOC_genes_summary_geneID.xlsx',
  sheet = 3
) %>%
  dplyr::select(Gene_name, At_gene_id, F.crenata, Pathway_name) %>%
  drop_na(Gene_name) %>%
  dplyr::rename(
    Gene_name_Original = Gene_name,
    AT_gene_id = At_gene_id
  )


BVOC.genes.variation = data.genes.leaves.monthly %>%
  right_join(
    isoform.cleand %>%
      dplyr::select(name, sseqid) %>%
      dplyr::rename(qseqid = name)
  ) %>%
  left_join(GO.gene.list %>% dplyr::select(qseqid, TAIR.ID) %>% distinct()) %>%
  right_join(
    meaningGenBVOC.function %>% dplyr::rename(TAIR.ID = AT_gene_id)
  ) %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  ) %>%
  drop_na(qseqid) %>%
  left_join(
    cv_sync_toghether %>%
      dplyr::select(qseqid, kCV, mean_synchrony, Category_kCV, Category_sync)
  ) %>%
  mutate(
    label = paste0(
      "kCV = ",
      round(kCV, 2),
      "\n",
      "Sync = ",
      round(mean_synchrony, 2)
    )
  )


BVOC.genes.variation %>%
  ggplot(aes(
    x = timepoint,
    y = level.exp,
    group = TreeID
  )) +
  facet_wrap(~Gene_name_Original, scales = "free") +
  geom_line(alpha = 0.6) +
  theme_minimal(base_size = 12) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  xlab("") +
  geom_vline(
    xintercept = which(levels(genes.temporal.filtering$timepoint) == "6-2014"),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = which(levels(genes.temporal.filtering$timepoint) == "6-2017"),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = which(levels(genes.temporal.filtering$timepoint) == "6-2021"),
    linetype = "dashed",
    color = "red"
  ) +
  ylab("Gene expression")
# Create unique timepoint label
sub.akikosul <- sub.akikosul %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  )

# Create axis labels: show year only at first occurrence per year
axis_labels_df <- sub.akikosul %>%
  distinct(timepoint, year, month) %>%
  arrange(year, month) %>%
  mutate(
    label = ifelse(
      !duplicated(year),
      paste0(month, "-", year),
      as.character(month)
    )
  )

axis_labels <- setNames(axis_labels_df$label, axis_labels_df$timepoint)

# Now plot
ggplot(
  sub.akikosul,
  aes(
    x = timepoint,
    y = level.exp,
    group = TreeID
  )
) +
  facet_wrap(~sseqid, scales = "free_y", ncol = 1) +
  geom_line(alpha = 0.8) +
  scale_x_discrete(labels = axis_labels) +
  theme_minimal(base_size = 12) +
  ggpubr::theme_pubr() +
  xlab("") +
  ylab("Gene expression") +
  hrbrthemes::theme_ipsum(base_size = 14, axis_title_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  geom_vline(
    xintercept = which(levels(sub.akikosul$timepoint) == "6-2014"),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = which(levels(sub.akikosul$timepoint) == "6-2017"),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = which(levels(sub.akikosul$timepoint) == "6-2021"),
    linetype = "dashed",
    color = "red"
  )

data.genes.all %>%
  #filter(qseqid == "Facr_v2.5_s1cl008229") %>%
  filter(qseqid == "Facr_v2.5_s1cl008231") %>%
  mutate(
    timepoint = paste(month, year, sep = "-"),
    timepoint = factor(
      timepoint,
      levels = unique(timepoint[order(year, month)])
    )
  ) %>%
  ggplot(aes(x = timepoint, y = level.exp, group = TreeID)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

######################################
#lets make netowkr genes
# dyn.load(here("dynGENIE3-master", "dynGENIE3_R_C_wrapper", "dynGENIE3.so"))
# getLoadedDLLs()
# library(doRNG)
# library(doParallel)
# library(reshape2)
# source("dynGENIE3-master/dynGENIE3_R_C_wrapper/dynGENIE3.R")
# file.exists("dynGENIE3-master/dynGENIE3_R_C_wrapper/dynGENIE3.R")
# summary_all_genes
# setwd("dynGENIE3-master/dynGENIE3_R_C_wrapper/")
# time_vec <- as.numeric(colnames(summary_all_genes))
# names(time_vec) <- as.character(time_vec)
# TS.data.mastsample = data_subset_norm.correlation
#
# res <- dynGENIE3(list(TS.data.mastsample), list(time_vec))
# link.list <- get.link.list(res$weight.matrix)
#
# top_edges <- as.data.frame(as.table(res$weight.matrix)) %>%
#   dplyr::rename(regulator = Var1, target = Var2, weight = Freq) %>%
#   filter(regulator != target) %>%
#   group_by(target) %>%
#   ungroup()
# quant <- quantile(top_edges$weight, 0.95) # or 0.95
# filtered_edges <- top_edges %>% filter(weight >= quant)
# # Create graph
# gene_network <- igraph::graph_from_data_frame(filtered_edges, directed = TRUE)
# plot(gene_network)
