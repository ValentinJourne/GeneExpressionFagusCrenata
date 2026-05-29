#####################################
#Gene expression analysis in Fagus crenata
#loading package and options
#####################################

#creator : Valentin Journ\'e; Kyushu University
#contact : journe.valentin@gmail.com
#file created on 2024-09-26
#updated on 2026 - 06

#####################################
#loading package
#####################################
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

#####################################
#loading options
#####################################
#set theme for ggplot
theme_set(theme_cowplot())
themesize = theme(
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 16),
  plot.margin = margin(0, 0, 0, 0)
)
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

#base_cols_treeid <- ggsci::pal_iterm("Rose Pine")(6)
#TreeID.color.pal = colorRampPalette(base_cols_treeid)(8)
TreeID.color.pal = rev(ggsci::pal_aaas()(8))
#ggsci::pal_aaas()(8)
#[1] "#3B4992FF" "#EE0000FF" "#008B45FF" "#631879FF" "#008280FF" "#BB0021FF" "#5F559BFF"
#[8] "#A20056FF"
genes.color <- ggsci::pal_aaas()(10)
genes.color.pal = colorRampPalette(genes.color)(13)


#RColorBrewer::brewer.pal(n = 12, name = "BuPu")
colors.pvalues <- colorRampPalette(c(
  "#C6DBEF",
  "#9ECAE1",
  "#6BAED6",
  "#4292C6",
  "#2171B5",
  "#08519C",
  "#08306B"
))(100)

less.vivid.color = colorRampPalette(c("#929CC2FF", "#F9D6ACFF", "#D4A6C9FF"))(
  10
)
