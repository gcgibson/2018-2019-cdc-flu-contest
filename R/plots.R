library(ggplot2)
library(lme4)
library(nnet)


lag_df <- read.csv("./data/lag_df")
lag_df0 <- read.csv("./data/lag_df0")

lm_fit <- lm(X0~Incidence +season_week + Region, data=lag_df0)
summary(lm_fit)
lm_fit_hierarchical <- lmer(X0~Incidence +season_week + (1|Region), data=lag_df0)
summary(lm_fit_hierarchical)
nn_fit <- nnet(X0~Incidence +season_week + Region, subset_lag_df, size=12, maxit=500, linout=T, decay=0.01)


library(ggplot2)

ggplot(lag_df,aes(x=1:length(X0),y=X0,col=Region)) +geom_point() + theme(axis.title.x=element_blank(),
                                                                            axis.text.x=element_blank(),
                                                                            axis.ticks.x=element_blank()) 
ggplot(lag_df,aes(x=1:length(X15),y=X15,col=Region)) +geom_point() + theme(axis.title.x=element_blank(),
                                                                         axis.text.x=element_blank(),
                                                                         axis.ticks.x=element_blank()) 

ggplot(lag_df,aes(x=1:length(X0),y=X51,col=Region)) +geom_point() + theme(axis.title.x=element_blank(),
                                                                         axis.text.x=element_blank(),
                                                                         axis.ticks.x=element_blank()) 


ggplot(lag_df,aes(x=season_week,y=X0,col=Region)) +geom_point() +theme_bw()


ggplot(lag_df[lag_df$Region=="National",],aes(x=season_week,y=X0)) +geom_point(color='maroon',alpha=.5) +theme_bw()
ggplot(lag_df0[lag_df0$Region=="National",],aes(x=Incidence,y=X0)) +geom_point(color='maroon',alpha=.5) +theme_bw()
ggplot(lag_df0,aes(x=Region,y=X0)) +geom_point(color='maroon',alpha=.5) +theme_bw()


predict_df <- data.frame(season_week = 1:53,Region=rep("National",53),week=1:53,X=1:53,Incidence=rep(3,53))
lm.predict <- predict(lm_fit,newdata = predict_df)

p <- ggplot(data.frame(season_week=1:53,X0=lm.predict), aes(x=season_week,y=X0)) 
p <- p+ geom_point(data=lag_df0,aes(x=season_week,y=X0),alpha=.1,col="cornflowerblue") + theme_bw()+geom_line()
p

lm.predict <- predict(lm_fit_hierarchical,newdata = predict_df)

p <- ggplot(data.frame(season_week=1:53,X0=lm.predict), aes(x=season_week,y=X0)) 
p <- p+ geom_point(data=lag_df0,aes(x=season_week,y=X0),alpha=.1,col="cornflowerblue") + theme_bw()+geom_line()
p

lm.predict <- predict(nn_fit,newdata = predict_df)

p <- ggplot(data.frame(season_week=1:53,X0=lm.predict), aes(x=season_week,y=X0)) 
p <- p+ geom_point(data=lag_df0,aes(x=season_week,y=X0),alpha=.1,col="cornflowerblue") + theme_bw()+geom_line()
p
