# ============================================================
# Shared style and export for Figures 1-7 (10 Nov 2026, review): one theme (Arial, theme_bw, bold title),
#   gene symbols in italics (markdown through ggtext; Ensembl IDs stay upright), P-values written as
#   a × 10^b, and the submission files. save_fig() writes the analysis/ PDF and PNG each plotting script has
#   always written (plus an SVG) and, when FIG_DIR is set (step 141), FigureN.pdf (vector, fonts embedded),
#   FigureN.svg and FigureN.tif (600 dpi, LZW) for the journal.
#   Sourced by 69b, 81c, 87, 76e, 89b and 70b.
# ============================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggtext) })

FIG_FONT   <- "Arial"
COL_NR1H3  <- "#E69F00"
COL_ZNF408 <- "#0072B2"
COL_SIG    <- "firebrick"
COL_EST    <- "#2C3E50"
FIG_W      <- 7.48    # 190 mm, the full page width of the journal: Figures 1-6 are drawn at print size
MINUS      <- "−"

theme_paper <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = FIG_FONT) +
    theme(panel.grid.minor   = element_blank(),
          plot.title         = element_markdown(face = "bold", size = rel(1.15)),
          plot.subtitle      = element_markdown(size = rel(0.9), colour = "grey25", lineheight = 1.25),
          plot.title.position = "plot",
          axis.title.x       = element_markdown(),
          axis.title.y       = element_markdown(),
          axis.text.y        = element_markdown(colour = "grey15"),
          axis.text.x        = element_text(colour = "grey15"),
          legend.title       = element_markdown(),
          strip.background   = element_rect(fill = "grey92", colour = "grey60"),
          strip.text         = element_markdown(face = "bold", hjust = 0))
}

# gene symbol(s) in italics for markdown text; Ensembl IDs are not symbols and stay upright
gi <- function(g) ifelse(is.na(g) | grepl("^ENSG", g), g, paste0("*", g, "*"))
# plain-text labels (geom_text / geom_text_repel): fontface for a gene label
gface <- function(g, bold = FALSE) ifelse(grepl("^ENSG", g), if (bold) "bold" else "plain",
                                          if (bold) "bold.italic" else "italic")

# axis labels with a true minus sign (scale_*_continuous(labels = lab_minus))
lab_minus <- function(x) sub("^-", MINUS, format(x, trim = TRUE, drop0trailing = TRUE))

# number with a true minus sign
fmt_num <- function(x, d = 2) sub("^-", MINUS, formatC(x, format = "f", digits = d))

# P-value for markdown: 0.37, 0.066, 0.0010, or 6.9 × 10<sup>−6</sup>
fmt_p <- function(p, sig = 2, cut = 1e-3) vapply(p, function(x) {
  if (is.na(x)) return("NA")
  if (x >= cut) return(formatC(signif(x, sig), format = "fg", digits = sig, flag = "#"))
  e <- floor(log10(x)); m <- round(x / 10^e, sig - 1)
  if (m >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("%s × 10<sup>%s%d</sup>", formatC(m, format = "f", digits = sig - 1), MINUS, -e)
}, character(1))

# same for plain-text labels (geom_text): superscript digits
fmt_p_plain <- function(p, sig = 2, cut = 1e-3) {
  sup <- function(s) chartr("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹", s)
  out <- fmt_p(p, sig, cut)
  m <- regmatches(out, regexec("^(.*)<sup>(.*)</sup>$", out))
  vapply(seq_along(out), function(i) if (length(m[[i]]) == 3)
    paste0(m[[i]][2], sub(MINUS, "⁻", sup(m[[i]][3]), fixed = TRUE)) else out[i], character(1))
}

# Broken y-axis: the plot is drawn twice, the upper panel over [brk[2], top] and the lower over [0, brk[1]],
# with out-of-range data dropped (scale limits), so points and labels appear only in their own segment; the
# break is marked with diagonal strokes on both axis lines. rel = height of the upper panel relative to the lower.
broken_y <- function(p, brk, top, rel = 0.3, upper_breaks = waiver()) {
  marks <- function(at) {
    s <- unit(3.5, "pt")
    annotation_custom(grid::segmentsGrob(x0 = unit(c(0, 1), "npc") - s, x1 = unit(c(0, 1), "npc") + s,
                                         y0 = unit(c(at, at), "npc") - s, y1 = unit(c(at, at), "npc") + s,
                                         gp = grid::gpar(col = "grey20", lwd = 1.1)))
  }
  up <- p + scale_y_continuous(limits = c(brk[2], top), breaks = upper_breaks, expand = expansion(mult = c(0.04, 0.08))) +
    marks(0) + coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(5.5, 5.5, 2, 5.5))
  lo <- p + scale_y_continuous(limits = c(0, brk[1]), expand = expansion(mult = c(0.02, 0.02))) +
    marks(1) + coord_cartesian(clip = "off") + labs(title = NULL, subtitle = NULL) +
    theme(plot.margin = margin(2, 5.5, 5.5, 5.5))
  patchwork::wrap_plots(up, lo, ncol = 1, heights = c(rel, 1)) + patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# Vector PDF: R's cairo devices need XQuartz on macOS and the base pdf() device cannot set Unicode text in Arial,
# so each figure is written as SVG (svglite; text kept as text) and printed to PDF at the figure's exact size by
# headless Google Chrome/Chromium, which embeds subsets of the fonts. Without Chrome the SVG and PNG/TIFF are
# still written and the PDF is skipped with a message. CHROME = path to the browser if not in a standard place.
find_chrome <- function() {
  cand <- c(Sys.getenv("CHROME"), "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            Sys.which(c("google-chrome", "google-chrome-stable", "chromium", "chromium-browser")))
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (length(cand)) cand[[1]] else ""
}
svg_to_pdf <- function(svg, pdf, width, height) {
  chrome <- find_chrome()
  if (!nzchar(chrome)) { message("No Chrome/Chromium found: ", pdf, " not written (SVG kept)"); return(invisible(FALSE)) }
  x <- paste(readLines(svg, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  x <- substring(x, regexpr("<svg", x, fixed = TRUE))
  html <- tempfile(fileext = ".html")
  writeLines(c("<!doctype html><html><head><meta charset=\"utf-8\"><style>",
               sprintf("@page{size:%sin %sin;margin:0}html,body{margin:0;padding:0}svg{display:block;width:%sin;height:%sin}",
                       width, height, width, height), "</style></head><body>", x, "</body></html>"), html, useBytes = TRUE)
  out <- normalizePath(pdf, mustWork = FALSE)
  prof <- tempfile("chrome_profile_")     # own profile per call: concurrent or interactive Chrome sessions can block
  unlink(out)
  system2(chrome, c("--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--no-first-run",
                    "--no-default-browser-check", paste0("--user-data-dir=", shQuote(prof)),
                    "--virtual-time-budget=30000",   # large SVGs (Figure 1: ~9,500 points) are printed before they finish drawing without it
                    paste0("--print-to-pdf=", shQuote(out)), shQuote(paste0("file://", html))),
          stdout = FALSE, stderr = FALSE, timeout = 180)
  unlink(c(html, prof), recursive = TRUE)
  if (!file.exists(out)) message("Chrome did not write ", pdf)
  invisible(file.exists(out))
}

save_fig <- function(p, stem, width, height, fig = NULL, png_dpi = 300) {
  ggsave(paste0(stem, ".svg"), p, width = width, height = height, device = svglite::svglite)
  svg_to_pdf(paste0(stem, ".svg"), paste0(stem, ".pdf"), width, height)
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = png_dpi, device = ragg::agg_png)
  fd <- Sys.getenv("FIG_DIR")
  if (!is.null(fig) && nzchar(fd)) {
    dir.create(fd, showWarnings = FALSE, recursive = TRUE)
    file.copy(paste0(stem, ".svg"), file.path(fd, paste0(fig, ".svg")), overwrite = TRUE)
    svg_to_pdf(file.path(fd, paste0(fig, ".svg")), file.path(fd, paste0(fig, ".pdf")), width, height)
    ggsave(file.path(fd, paste0(fig, ".tif")), p, width = width, height = height, dpi = 600,
           device = ragg::agg_tiff, compression = "lzw")
  }
  invisible(stem)
}
