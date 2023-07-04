rm(list = ls())
library(tidyverse)
library(ggrepel)
load(file = "data/rdata/deg_GSE82107_deg_oavsctrl.Rdata")
DEG_all <- DEGs_ebayes %>%
    na.omit() %>%
    mutate(
        symbol = rownames(DEGs_ebayes),
        logFC,
        Pvalue = P.Value
    ) %>%
    dplyr::select(
        symbol,
        logFC,
        Pvalue
    ) %>%
    mutate(
        direction = factor(
            ifelse(Pvalue < 0.05 & abs(logFC) > 0.5,
                ifelse(logFC > 0.5, "Up", "Down"), "NS"
            ),
            levels = c("Up", "Down", "NS")
        )
    )
pdf(
    file = "res/pic/final/Volcano_plot_gse82107.pdf", # nolint
    width = 16,
    height = 8
)
par(mfrow = c(1, 1))
ggplot(
    data = DEG_all,
    aes(
        x = logFC,
        y = -log10(Pvalue),
        colour = direction
    )
) +
    geom_point(alpha = 0.6) +
    scale_color_manual(
        values = c("#DC143C", "#00008B", "#808080")
    ) +
    geom_text_repel(
        data = DEG_all %>%
            filter(
                Pvalue < 0.05,
                abs(logFC) > 1
            ),
        aes(
            label = symbol
        ),
        size = 3,
        segment.color = "black",
        show.legend = FALSE
    ) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    ylab(
        expression(
            -log[10]("Adjusted P Value")
        )
    ) +
    xlab(
        expression(
            log[2]("Fold Change")
        )
    ) +
    xlim(-5, 5) +
    ylim(0, 5) +
    geom_vline(
        xintercept = c(-0.5, 0.5),
        lty = 2,
        col = "black",
        lwd = 0.6
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        lty = 2,
        col = "black",
        lwd = 0.6
    )
dev.off()
