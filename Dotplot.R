features <- list(
  "cDC-CLEC9A" = c("CLEC9A", "CCSER1", "DNASE1L3"),
  "cDC-CD1C" = c("CD1C", "PPA1", "HLA-DPB1"),
  "pDC-LILRA4" = c("LILRA4", "PTPRS", "JCHAIN"),
  "Macrophage-MKI67" = c("MKI67", "H2AFZ", "TUBB"),
  "Macrophage-IL1B" = c("IL1B", "CXCL2", "NFKBIA"),
  "Macrophage-CD209" = c("CD209", "F13A1", "MSR1"),
  "cDC-LAMP3" = c("LAMP3", "CCL22", "CCR7"),
  "Macrophage-IDO1" = c("IDO1", "NAPSB", "CLNK"),
  "Macrophage-APOE" = c("APOE", "S100A6", "FBP1"),
  "Neutrophil" = c("FCGR3B", "CSF3R", "CXCL8"),
  "Macrophage-C1QA" = c("C1QA", "C1QB", "C1QC")
)
DotPlot(object = Myeloid, features = features) &
  theme_bw() &
  geom_point(shape = 21, aes(size = pct.exp), stroke = 1) &
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(color = 'black', angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(color = 'black', face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), 'cm'),
    panel.border = element_rect(color = "black", size = 1.2, linetype = "solid"),
    panel.spacing = unit(0.12, "cm"),
    legend.box.background = element_rect(colour = "black"),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.title = element_text(color = 'black', face = "bold", size = 9)
  ) &
  scale_color_gradientn(colours = colorRampPalette(c("#C2ADCE", "white", "#ECA8A9"))(100)) &
  labs(tag = "") &
  theme(
    plot.tag.position = c(0.3, 1.05),
    plot.tag = element_text(size = 12, face = "bold")
  ) &
  guides(
    size = guide_legend(title = "Proportion of\nexpressing cells"),
    colour = guide_colorbar(title = "Average\nexpression")
  )