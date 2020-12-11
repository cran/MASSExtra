## ----presetup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "",
                      message = FALSE,
                      warning = FALSE)
old_plot_new <- getHook("plot.new")
old_before_plot_new <- getHook("before.plot.new")
setHook("plot.new", list(), "replace")
setHook("before.plot.new", rev(old_before_plot_new)[1], "replace")
old_par <- par(no.readonly = TRUE)
old_options <- options()

library(tidyverse)
old_theme <- theme_get()

## ----packages_etc-------------------------------------------------------------
suppressPackageStartupMessages({
  library(visreg)
  library(knitr)
  library(tidyverse)
  library(patchwork)
  library(MASSExtra)
})
options(knitr.kable.NA = "")
theme_set(theme_bw() + theme(plot.title = element_text(hjust = 0.5)))

## ----setup, fig.height = 5, fig.width = 10, out.width="100%", fig.cap="Box-cox, old and new displays"----
par(mfrow = c(1, 2))
mod0 <- lm(MPG.city ~ Weight, Cars93)
boxcox(mod0)  ## MASS
box_cox(mod0) ## MASSExtra tweak

## ---- fig.height=5, fig.width=10, out.width="100%", fig.cap="The Box-Cox transformation effect"----
p0 <- ggplot(Cars93) + aes(x = Weight) + geom_point(colour = "#2297E6")  + xlab("Weight (lbs)") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, size=0.7, colour = "black")
p1 <- p0 + aes(y = MPG.city) + ylab("Miles per gallon (MPG)") + ggtitle("Untransformed response")
p2 <- p0 + aes(y = bc(MPG.city, lambda(mod0))) + ggtitle("Transformed response") + 
  ylab(bquote(bc(MPG, .(round(lambda(mod0), 2)))))
p1 + p2

## ---- fig.width=8, fig.height=6, fig.align="center", out.width="75%"----------
big_model <- lm(medv ~ . + (rm + tax + lstat + dis)^2 + poly(dis, 2) + poly(rm, 2) +
                   poly(tax, 2) + poly(lstat, 2), Boston)
big_model %>% drop_term(k = "GIC") %>% plot() %>% kable(booktabs=TRUE, digits=3)

## ---- fig.width=8, fig.height=6, fig.align="center", out.width="75%"----------
base_model <- lm(medv ~ ., Boston)
gic_model <- step_GIC(base_model, scope = list(lower = ~1, upper = formula(big_model)))
drop_term(gic_model) %>% plot() %>% kable(booktabs = TRUE, digits = 3)

## ---- fig.width=14, fig.height=6, fig.align="center", out.width="100%", fig.cap="Two component visualisations for the Boston house price model"----
capture.output(suppressWarnings({
   g1 <- visreg(gic_model, "dis", plot = FALSE)
   g2 <- visreg(gic_model, "lstat", plot = FALSE)
   plot(g1, gg = TRUE) + plot(g2, gg = TRUE) 
})) -> void

## ---- fig.width=8, fig.height=6, fig.align="center", out.width="75%"----------
add_term(base_model, scope = formula(big_model), k = "gic") %>% 
   plot() %>% kable(booktabs = TRUE, digits = 3)

## -----------------------------------------------------------------------------
quine_full <- glm.nb(Days ~ Age*Eth*Sex*Lrn, data = quine)
drop_term(quine_full) %>% kable(booktabs = TRUE, digits = 4)
quine_gic <- step_GIC(quine_full)
drop_term(quine_gic) %>% kable(booktabs = TRUE, digits = 4)

## -----------------------------------------------------------------------------
mpg0 <- glm(MPG.city ~ Weight + Cylinders + EngineSize + Origin, 
            family = Gamma(link = "inverse"), data = Cars93)
drop_term(mpg0) %>% kable(booktabs = TRUE, digits = 3)
mpg_gic <- step_down(mpg0, k = "gic")      ## simple backward elimination, mainly used for GLMMs
drop_term(mpg_gic) %>% kable(booktabs = TRUE, digits = 3)

## ---- fig.width=10, fig.height=7, fig.align="center", out.width="80%", fig.cap="MPG model using a Gamma family with inverse link"----
ggplot(Cars93) + aes(x = Weight, y = MPG.city, colour = Cylinders) + geom_point() +
   geom_smooth(method = "glm", method.args = list(family = Gamma), formula = y ~ x) +
   ylab("Miles per gallon (city driving)") + scale_colour_brewer(palette = "Dark2")

## -----------------------------------------------------------------------------
dput(MASSExtra:::makepredictcall.normalise) ## exported only as an S3 method
dput(.normalise)                            ## an exported helper, but with a hidden name
dput(zs)                                    ## the scale() equivalent
dput(zq)                                    ## the 'boxplot scaling' (based on quantiles)

## ---- echo=FALSE--------------------------------------------------------------
setdiff(getNamespaceExports("MASS"), getNamespaceExports("MASSExtra")) %>% sort() %>%  noquote()

## ---- echo=FALSE--------------------------------------------------------------
intersect(getNamespaceExports("MASS"), getNamespaceExports("MASSExtra")) %>% sort() %>%  noquote()

## ---- echo=FALSE--------------------------------------------------------------
setdiff(getNamespaceExports("MASSExtra"), getNamespaceExports("MASS")) %>% 
   setdiff(getNamespaceExports("splines")) %>% 
   grep("^[.]__T", ., invert = TRUE, value = TRUE) %>% ## exclude S4 class objects
   sort() %>% noquote()

## ---- echo=FALSE--------------------------------------------------------------
intersect(getNamespaceExports("MASSExtra"), getNamespaceExports("splines")) %>% sort() %>%  noquote()

## ----cleanup, include=FALSE---------------------------------------------------
setHook("plot.new", old_plot_new, "replace")
setHook("before.plot.new", old_before_plot_new, "replace")
par(old_par)
options(old_options)
theme_set(old_theme)

