rm(list = ls())
library(tidyverse)
library(RColorBrewer)
load(file = "data/rdata/matrix_GSE82107_raw.Rdata")
Expr_gse82107 <- Expr_GSE82107_anno
gene <- "EGR1"
Targets <- Targets %>%
    rownames_to_column("sample_id")
plot_df <- Expr_gse82107[gene, ] %>%
    t() %>%
    as.data.frame() %>%
    set_names("value") %>%
    rownames_to_column("sample_id") %>%
    left_join(Targets, by = "sample_id") %>%
    dplyr::select(sample_id, group, value)
plot_df$group <- factor(
    plot_df$group,
    levels = c("OA", "CTRL")
)
p <- t.test(
    plot_df[which(plot_df$group == "OA"), "value"],
    plot_df[which(plot_df$group == "CTRL"), "value"]
)$p.value
jco <- c("#1e62ac", "#dc3b3b")
ggplot(
    data = plot_df,
    aes(x = group, y = value, fill = group)
) +
    scale_fill_manual(values = jco) +
    geom_boxplot(
        outlier.size = -1,
        color = "black",
        lwd = 0.6,
        alpha = 0.7
    ) +
    geom_point(
        shape = 21,
        size = 5,
        position = position_jitterdodge(),
        color = "black",
        alpha = 1
    ) +
    theme_classic() +
    ylab(paste0("Expression value of ", gene)) +
    xlab("") +
    annotate(
        geom = "text",
        cex = 4,
        x = 1.5,
        y = 15,
        label = paste0(
            "P ",
            ifelse(p < 0.001, "< 0.001",
                paste0("= ", round(p, 3))
            )
        ),
        color = "black"
    ) +
    theme(
        panel.border = element_rect(
            colour = "black",
            fill = NA, size = 0.2
        ),
        axis.ticks = element_line(
            size = 0.2,
            color = "black"
        ),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
    )
fn <- paste0("res/pic/final/boxplot_of_", gene, "_gse82107.pdf")
ggsave(fn,
    width = 4.5, height = 4
)
