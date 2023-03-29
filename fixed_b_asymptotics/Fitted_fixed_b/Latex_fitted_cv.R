library(xtable)

fitted_cv <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/fixed_b_asymptotics/Fitted_fixed_b/fitted_CV.csv"
the_fits <- read.csv(fitted_cv)
the_fits <- the_fits[, -1]
the_fits <- the_fits[the_fits$d==1, ]
the_fits <- the_fits[,-6]

the_kernel <- "Bartlett"
latex_table <- the_fits[the_fits$kernel == "Bartlett",]
latex_table <- latex_table[,-6]
lug_index <- latex_table$is_lugsail
alpha <- latex_table$alpha[!lug_index]
latex_table <- latex_table[,-c(6,7)]
latex_table <- cbind(alpha,
                     latex_table[!lug_index,], 
                     latex_table[lug_index,])

latex_table <- xtable(latex_table)


colnames(latex_table) <- c("$\\\alpha$", rep(c("$a_0$", "$a_1$", "$a_2$", "$a_3$", "$R^2$"), 2))
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- c(" \\hline &  & \\multicolumn{5}{c|}{\textit{Mother}}& \\multicolumn{2}{c|}{\textit{Zero}}  \\\\\n")

rownames(latex_table) <- c("Bartlett", "", "  ", "   ")

italic <- function(x){
  paste0('{\\emph{ ', x, '}}')
}

align(latex_table) <- "|cc||ccccc||ccccc|"

print(latex_table, sanitize.rownames.function = italic,
      add.to.row = addtorow,
      hline.after = c(-1,  4), 
      include.colnames = T)
