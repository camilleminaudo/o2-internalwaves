



dt <- read.csv("C:/Projects/myGit/o2-internalwaves/data/Gergal_example_Camille.csv")


head(dt)


dt$time = as.POSIXct(paste0(dt$pr_dt," ", dt$time_GMT), tz = "GMT")


names(dt)

df_wide <- dt[,c("time",
                 "DO_1m", "DO_9m", "DO_19m",
                 "wt_1m", "wt_9m", "wt_19m",
                 "PAR_1m", "PAR_9m", "PAR_19m",
                 "Viento_ms",
                 "Zmax_hydro")]
names(df_wide) <- c("time",
                    "DO_1", "DO_9", "DO_19",
                    "T_1", "T_9", "T_19",
                    "PAR_1", "PAR_9", "PAR_19",
                    "wind",
                    "z_th")

# df_wide has the correct format for running the do_budget_metabolism.R script
write.csv(x = df_wide, file = "C:/Projects/myGit/o2-internalwaves/data/Gergal_example_wide.csv", row.names = F)




  df_1 <- dt[,c("time","DO_1m","wt_1m")]
names(df_1) <- c("time","DO","T")
df_1$depth = 1

df_9 <- dt[,c("time","DO_9m","wt_9m")]
names(df_9) <- c("time","DO","T")
df_9$depth = 9

df_19 <- dt[,c("time","DO_19m","wt_19m")]
names(df_19) <- c("time","DO","T")
df_19$depth = 19

df <- rbind(df_1, df_9, df_19)

summary(df)


# df has the correct format for running directly the internal wave correction time domain script.
write.csv(x = df, file = "C:/Projects/myGit/o2-internalwaves/data/Gergal_example_Camille_formatted.csv", row.names = F)


