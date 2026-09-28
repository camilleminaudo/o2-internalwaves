



dt <- read.csv("C:/Projects/myGit/o2-internalwaves/data/Gergal_example_Camille.csv")


head(dt)


dt$time = as.POSIXct(paste0(dt$pr_dt," ", dt$time_GMT), tz = "GMT")


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

write.csv(x = df, file = "C:/Projects/myGit/o2-internalwaves/data/Gergal_example_Camille_formatted.csv", row.names = F)


