library(tidyverse)
library(vegan)
library(GGally)
library(FactoMineR)
library(factoextra)
library(viridis)
library(corrplot)
library(car)
library(mgcv)
library(mgcViz)
library(foreach)

main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20, hjust = 0.5),
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=28))

read.table("data/data_clean/Metadata_quadra_CAM.txt", header  = T)%>%
  filter(!str_detect(Host_code, "22_CAM_15...."))-> metadata_quad
read.table("data/data_clean/Metadata_grid_CAM.txt", header  = T)%>%
  filter(!str_detect(Grid_code, "22_CAM_15"))-> metadata_grid

metadata_quad%>%
  arrange(Host_code)%>%
  select(!c(Ecosystem, Collection_date, Locality, Country))%>%
  right_join(metadata_grid, by =  "Grid_code") -> all_metadata

str(all_metadata)
##### Richess and biomass #####
##### FAMD and glm ##### 

##### more convignente setup for data
Y = all_metadata%>%
  select(Biomass, Plant_richness)
Y_log = all_metadata%>%
  reframe(log_biomass = log1p(Biomass),
         log_plant_richness = log1p(Plant_richness))

X =  all_metadata%>%
  dplyr::select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                                   CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                depth_oxy, Pature, Fauche, natural_landscape, cultivated, artificial, non_emitting,wetland)

##### FAMD to extract env linear combination
X %>%
  select(-c(natural_landscape, cultivated, artificial, non_emitting,wetland))%>%
  mutate(across(c(depth_oxy,Fauche, Pature), as.factor),
         across(where(is.numeric), ~decostand(., method = "standardize"))) -> metadata_grid_standard

FAMD_env <- FAMD(metadata_grid_standard, graph = F)
env_coord_FAMD = FAMD_env$ind$coord
fviz_contrib(FAMD_env, axes = c(1,2,3,4,5), choice  = "var")
X11()
fviz_contrib(FAMD_env, axes = c(1), choice  = "var")
X11()
fviz_contrib(FAMD_env, axes = c(2), choice  = "var")
X11()
fviz_contrib(FAMD_env, axes = c(3), choice  = "var")
X11()
fviz_contrib(FAMD_env, axes = c(4), choice  = "var")
X11()
fviz_contrib(FAMD_env, axes = c(5), choice  = "var")
dev.off()

corrplot(cor(X%>%select(where(is.numeric))),method = "shade", hclust.method = "ward.D2")
FAMD_env$var
fviz_famd_var(FAMD_env,
              choice = "quanti.var", # show points only (nbut not "text")
              col.var = "cos2", # color by groups
             axes = c(1,2),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

fviz_famd_ind(FAMD_env,
              col.ind = "cos2", # color by groups
              axes = c(1,2),
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

cluster_3 = hcut(Y, hc_metric = "euclidean", hc_method = "ward.D2", k =3)
cluster3_quad_df = data.frame(cluster_3$data, clust = as.factor(cluster_3$cluster))


##### Richness and biomass a first overview #####
df_glm = cbind(Y_log,
    env_coord_FAMD)

ggpairs(df_glm #, ggplot2::aes(colour=cluster3_quad_df$clust)
)

X_selected = X%>%
  select(clay, lime_tot, depth_oxy, pHwater, MO, Cond)%>%
  mutate(depth_oxy =  factor(X_selected$depth_oxy, level = c( "0-9","10-19", "20-29", "30-39", ">40")))
  #mutate(across(where(is.numeric), ~decostand(., method = "standardize")))%>%
  #cbind(Y_log)



##### biomass glm ######
cbind(Biomass = Y$Biomass,X_selected%>%select(!depth_oxy))%>%
  pivot_longer(-Biomass)%>%
  ggplot()+
  facet_wrap(~name, scale = "free")+
  geom_point(aes(value, Biomass))+
  #geom_smooth(aes(value, Biomass))+
  main_theme
# full_biom = glm(log_biomass~Dim.1+Dim.2+Dim.3+Dim.4+Dim.5, data = df_glm, family = Gamma(link = "log"))
# summary(full_biom)
# hist(full_biom$residuals)
# (full_biom$null.deviance-full_biom$deviance)/full_biom$null.deviance

selected_biom_glm = glm(Y$Biomass ~ depth_oxy + clay  + lime_tot + pHwater +
                          MO + Cond, data = X_selected%>%
                          mutate(across(where(is.numeric), ~decostand(., method = "standardize"))),
                        family = Gamma(link = "log"))
summary(selected_biom_glm)

(selected_biom_glm$null.deviance-selected_biom_glm$deviance)/selected_biom_glm$null.deviance

hist(selected_biom_glm$residuals)


step_biom_glm = step(selected_biom_glm, direction = "both")

step_biom_glm_summary = summary(step_biom_glm)
step_biom_glm_summary
step_biom_glm_summary$coefficients%>%
  round(3)%>%
  knitr::kable()

(step_biom_glm$null.deviance-step_biom_glm$deviance)/step_biom_glm$null.deviance
vif(step_biom_glm)



crPlots(step_biom_glm, ylab = "Partial response/residuals (Biom)",
        col.lines=c("red", "blue"))

mod = varpart(Y = Y$Biomass, ~lime_tot, ~pHwater,~
          MO,~ Cond, data = X_selected)
plot(mod, bg = c("hotpink","skyblue"))

cbind(Biomass = Y$Biomass,X_selected%>%select(!depth_oxy),
      fit = fitted(step_biom_glm))%>%
  pivot_longer(-c(Biomass,fit))%>%
  ggplot()+
  facet_wrap(~name, scale = "free")+
  geom_point(aes(value, Biomass))+
  geom_line(aes(value, fit))+
  main_theme

# Convert crPlots object to ggplot
str(cr_plot_data)

##### Biomass gam #####

selected_biom_gam1 <- gam(Y$Biomass ~  s(lime_tot) + s(clay) + depth_oxy +
                               s(pHwater) + s(MO) + s(Cond), data = X_selected,family = Gamma("log"))

selected_biom_gam2 <- gam(Y$Biomass ~ s(lime_tot) + s(clay) +
                           s(pHwater) + s(MO) + s(Cond), data = X_selected,family = Gamma("log"))

selected_biom_gam3 <- gam(Y$Biomass ~  s(lime_tot) + 
                            s(pHwater) + s(MO) + s(Cond), data = X_selected,family = Gamma("log"))

AIC(step_biom_glm, selected_biom_gam1, selected_biom_gam2, selected_biom_gam3,selected_biom_gam4)%>%
  knitr::kable()

par(mfrow = c(2,2))
gam.check(selected_biom_gam2)
summary(selected_biom_gam2)
temp = summary(selected_biom_gam2)
str(temp)
temp$s.table%>%
  round(3)%>%
  knitr::kable()
k.check(selected_biom_gam2)

par(mfrow = c(3,2))
p_obj_biom <- plot(selected_biom_gam2, residuals = TRUE,all.terms = T)


sm_biom_df <- foreach(i = 1:length(p_obj_biom), .combine=rbind )%do%{
  as.data.frame(p_obj_biom[[i]][c("x", "se", "fit","xlab")])
  
}
data_biom_df <- foreach(i = 1:length(p_obj_biom), .combine=rbind )%do%{
  temp = as.data.frame(p_obj_biom[[i]][c("raw", "p.resid", "xlab")])
  temp = cbind(temp, Y)
  temp
} 
## plot

ggplot(sm_biom_df, aes(x = x, y = fit)) +
  facet_wrap("xlab", scale = "free")+
  geom_rug(data = data_biom_df, mapping = aes(x = raw, y = NULL),
           sides = "b") +
  geom_smooth(data =data_biom_df, aes(x = raw, y = p.resid),method = "lm",se=F, col = "red")+
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 0.1) +
  geom_point(data = data_biom_df, mapping = aes(x = raw, y = p.resid, size = Plant_richness ),
             alpha = 0.5) +
  geom_line(col = "blue") +
  
  labs(x ="Env values", y = "Partial response/resid (Biomass)", col = "richn", size = "biom")+
  scale_color_viridis()+
  main_theme


##### Richness model #####
# full_richn = glm(log_plant_richness~Dim.1+Dim.2+Dim.3+Dim.4+Dim.5, data = df_glm, family = Gamma(link = "log"))
# summary(full_richn)
# hist(full_richn$residuals)
# (full_richn$null.deviance-full_richn$deviance)/full_richn$null.deviance
# Biomasse correlated with DIm 1, 2 and 3
# Richness corralated with dim 4 and 5

###
selected_rich_glm = glm(Y$Plant_richness ~ depth_oxy + clay  + lime_tot + pHwater + MO + Cond, 
                        data = X_selected, family = Gamma(link = "log"))
step(selected_rich_glm)
summary(selected_rich_glm)
summary(selected_rich_glm)$coefficients%>%
  round(3)%>%
  knitr::kable()
hist(selected_rich_glm$residuals, breaks = 40)

(selected_rich_glm$null.deviance-selected_rich_glm$deviance)/selected_rich_glm$null.deviance

crPlots(selected_rich_glm,col.lines=c("red", "blue"))
summary(glm(Y$Plant_richness ~ depth_oxy + clay  + lime_tot + pHwater + MO + Cond, 
                        data = X_selected%>%
              mutate(across(where(is.numeric), ~decostand(., method = "range"))), family = Gamma(link = "log")))
vif(selected_rich_glm)
#####
##### Ratio biom/rich #####

Y%>%mutate(
  ratio_RB = Plant_richness/Biomass) -> Y
ggplot()+
  geom_point(data = Y, aes(Biomass, Plant_richness, col = log1p(ratio_RB)), cex= 2)+
  scale_color_viridis()+
  labs(y = "Plant richness", col = "log(1+Rich/biom)")+
  main_theme
hist(Y$ratio_RB, breaks = 40)

Y_log%>%mutate(
  log_ratio_RB = log_plant_richness/log_biomass) -> Y_log
ggplot()+
  geom_point(data = Y_log, aes(log_biomass, log_plant_richness, col = log_ratio_RB), cex= 2)+
  scale_color_viridis()+
  main_theme
hist(Y_log$log_ratio_RB)


ggpairs(cbind(Y_log,X_selected),
        lower = list(continuous = "points", combo = "box_no_facet", discrete = "count", na = "na"),
        upper = list(continuous = "cor")#ggplot2::aes(colour=CAH_env_FAMD$data.clust$clust)
      # Define color based on correlation strength
)

cbind(ratio_RB = Y$ratio_RB,X_selected%>%select(!depth_oxy))%>%
  pivot_longer(-ratio_RB)%>%
  ggplot()+
  facet_wrap(~name, scale = "free")+
  geom_point(aes(value, ratio_RB))+
  geom_smooth(aes(value, ratio_RB))

cbind(ratio_RB = Y$Biomass,X_selected)%>%
  ggplot()+
  geom_jitter(aes(depth_oxy,ratio_RB), width = 0.1)

##### Ratio biom/rich models #####

selected_ratioBR_glm = glm(Y$ratio_RB ~ depth_oxy + clay  + lime_tot + pHwater + MO + Cond, data = X_selected, family = Gamma(link = "log"))
summary(selected_ratioBR_glm)
step_ratioBR_glm = step(selected_ratioBR_glm, direction ="both")
summary(step_ratioBR_glm)
par(mfrow = c(2,2))
plot(step_ratioBR_glm)
par(mfrow = c(1,1))
hist(step_ratioBR_glm$residuals,breaks = 40)

#####
ID_test_batch = sample(1:369, round(369*0.2))
test_batch_X = X_selected[ID_test_batch,]
test_batch_Y = Y[ID_test_batch,]
training_batch_X = X_selected[-ID_test_batch,]
training_batch_Y = Y[-ID_test_batch,]


selected_ratioBR_glm_2 = glm(log1p(training_batch_Y$ratio_RB) ~ depth_oxy + clay  +
                             lime_tot + pHwater + MO + Cond,
                           data = training_batch_X, family =gaussian(link = "identity"))

summary(selected_ratioBR_glm_2)
step_ratioBR_glm_2 = step(selected_ratioBR_glm, direction ="both")
summary(step_ratioBR_glm)
preticed_ratio = predict.glm(step_ratioBR_glm, newdata = test_batch_X) 

plot(log1p(test_batch_Y$ratio_RB), preticed_ratio, xlim = c(-2.5,2.5),  ylim = c(-0.5,2.5))
abline(0,1)

##### gam #####

selected_ratioBR_gam <- gam(Y$ratio_RB ~  depth_oxy + clay  + lime_tot + 
                              pHwater + MO + Cond, data = X_selected, family = Gamma("log"))
par(mfrow = c(3,2))
summary(selected_ratioBR_gam)

selected_ratioBR_gam2 <- gam(Y$ratio_RB ~  s(lime_tot) + s(clay) +
                              s(pHwater) + s(MO) + s(Cond), data = X_selected,family = Gamma("log"))

par(mfrow = c(2,2))
gam.check(selected_ratioBR_gam2)

summary(selected_ratioBR_gam2)
par(mfrow = c(3,2))
str(plot(selected_ratioBR_gam2, all.terms = T))
AIC(selected_ratioBR_gam,selected_ratioBR_glm, selected_ratioBR_gam2)
k.check(selected_ratioBR_gam2)

par(mfrow = c(3,2))
p_obj <- plot(selected_ratioBR_gam2, residuals = TRUE,all.terms = T)

X_selected$clay
sm_df <- foreach(i = 1:length(p_obj), .combine=rbind )%do%{
  as.data.frame(p_obj[[i]][c("x", "se", "fit","xlab")])
  
}
data_df <- foreach(i = 1:length(p_obj), .combine=rbind )%do%{
  temp = as.data.frame(p_obj[[i]][c("raw", "p.resid", "xlab")])
  temp = cbind(temp, Y)
  temp
} 
## plot

ggplot(sm_df, aes(x = x, y = fit)) +
  facet_wrap("xlab", scale = "free")+
  geom_rug(data = data_df, mapping = aes(x = raw, y = NULL),
           sides = "b") +
  
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 0.3) +
  geom_point(data = data_df, mapping = aes(x = raw, y = p.resid, col =Plant_richness, size = Biomass ),
            alpha = 0.7) +
  geom_line() +
  labs(x ="Env values", y = "Partial response (richn per biom)", col = "richn", size = "biom")+
  scale_color_viridis()+
  main_theme

# plot(selected_ratioBR_gam2)
# b <- getViz(selected_ratioBR_gam2)
# o <- plot( sm(b, 1) )
# o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
#   l_ciLine(mul = 5, colour = "blue", linetype = 2) +
#   l_points(shape = 19, size = 1, alpha = 0.1) +
#   labs(x = "Clay", y = "biomass per species")
#   main_theme
# b <- getViz(selected_ratioBR_gam2)
# str(b)
# plot(selected_ratioBR_gam2$fitted.values)
# print(plot(b) + l_fitLine(colour ="blue") + l_ciLine() + l_rug() + l_points(), pages =1)
# 
# pl <- plot(b, allTerms = T) + l_points() +
#   l_ciLine(colour = 2) + l_ciBar() + theme_get() + labs(title = NULL)
# 


##### test with interaction #####
selected_ratioBR_gam3 <- gam(log1p(Y$ratio_RB) ~  as.factor(depth_oxy) + s(clay,Cond)  + s(lime_tot) +
                               s(pHwater) + s(MO) , data = X_selected)
par(mfrow = c(1,1))
vis.gam(selected_ratioBR_gam3, theta = 40, n.grid = 50, lwd = 0.4)
AIC(selected_ratioBR_gam2,selected_ratioBR_gam3)
k.check(selected_ratioBR_gam3)
plot(selected_ratioBR_gam3, page = 1, scheme = 2)
summary(selected_ratioBR_gam3)

##### 
selected_ratioBR_gam2_test <- gam(training_batch_Y$ratio_RB ~ as.factor(depth_oxy) + s(clay)  + s(lime_tot) +
                               s(pHwater) + s(MO) + s(Cond), data = training_batch_X,
                               family = Gamma("log"))

summary(selected_ratioBR_gam2_test)

preticed_ratio_gam = predict.gam(selected_ratioBR_gam2_test, newdata = test_batch_X) 

preticed_ratio_gam
par(mfrow = c(1,1))
plot(test_batch_Y$ratio_RB, exp(preticed_ratio_gam))
abline(0,1, lty = 2, col ="red")
mean(abs(log1p(test_batch_Y$ratio_RB) - preticed_ratio_gam))
par(mfrow = c(2,2))

gam.check(selected_ratioBR_gam2_test)

model_p = predict_gam(selected_ratioBR_gam2)

model_p 
##### Ratio biom/rich  log models #####

selected_ratioBR_glm_log = glm(Y_log$log_ratio_RB ~ depth_oxy + clay  + lime_tot + pHwater + MO + Cond, data = X_selected, family =gaussian(link = "identity"))
summary(selected_ratioBR_glm_log)
step_ratioBR_glm_log = step(selected_ratioBR_glm_log, direction ="both")
summary(step_ratioBR_glm_log)
par(mfrow = c(2,2))
plot(step_ratioBR_glm_log)
par(mfrow = c(1,1))
hist(step_ratioBR_glm_log$residuals)

#####
ID_test_batch = sample(1:369, round(369*0.2))
test_batch_X = X_selected[ID_test_batch,]
test_batch_Y_log = Y_log[ID_test_batch,]
training_batch_X = X_selected[-ID_test_batch,]
training_batch_Y_log = Y_log[-ID_test_batch,]


selected_ratioBR_glm_log_2 = glm(training_batch_Y_log$log_ratio_RB ~ depth_oxy + clay  +
                               lime_tot + pHwater + MO + Cond,
                             data = training_batch_X, family =gaussian(link = "identity"))

summary(selected_ratioBR_glm_log_2)
step_ratioBR_glm_log_2 = step(selected_ratioBR_glm_log_2, direction ="both")
summary(step_ratioBR_glm_log_2)
preticed_ratio = predict.glm(step_ratioBR_glm_log_2, newdata = test_batch_X) 

plot(log1p(test_batch_Y_log$log_ratio_RB), preticed_ratio, xlim = c(-0.1,2.5),  ylim = c(-0.5,2.5))
abline(0,1)

##### gam log #####

selected_ratioBR_gam_log <- gam(Y_log$log_ratio_RB ~  depth_oxy + clay  + lime_tot + 
                              pHwater + MO + Cond, data = X_selected)
selected_ratioBR_gam2_log <- gam(Y_log$log_ratio_RB ~ as.factor(depth_oxy) + s(clay)  + s(lime_tot) +
                               s(pHwater) + s(MO) + s(Cond), data = X_selected)
summary(selected_ratioBR_gam2_log)
par(mfrow = c(3,2))
plot(selected_ratioBR_gam2_log, all.terms = T)
AIC(selected_ratioBR_gam_log,selected_ratioBR_glm_log, selected_ratioBR_gam2_log)

selected_ratioBR_gam3_log <- gam(Y_log$log_ratio_RB ~  as.factor(depth_oxy) + s(clay,Cond)  + s(lime_tot) +
                               s(pHwater) + s(MO) , data = X_selected)
par(mfrow = c(1,1))
vis.gam(selected_ratioBR_gam3, theta = 40, n.grid = 50, lwd = 0.4)
AIC(selected_ratioBR_gam2_log,selected_ratioBR_gam3_log)
k.check(selected_ratioBR_gam3)
plot(selected_ratioBR_gam3, page = 1, scheme = 2)
summary(selected_ratioBR_gam3)

##### 
selected_ratioBR_gam2_log_test <- gam(training_batch_Y_log$log_ratio_RB ~ as.factor(depth_oxy) + s(clay)  + s(lime_tot) +
                                    s(pHwater) + s(MO) + s(Cond), data = training_batch_X)

summary(selected_ratioBR_gam2_log_test)

preticed_ratio_gam_log = predict.gam(selected_ratioBR_gam2_log_test, newdata = test_batch_X) 

preticed_ratio_gam_log

plot(test_batch_Y_log$log_ratio_RB, preticed_ratio_gam_log, xlim = c(0,3),  ylim = c(0,3))
abline(0,1, lty = 2, col ="red")
mean(abs(test_batch_Y$ratio_RB - preticed_ratio_gam))
par(mfrow = c(2,2))
gam.check(selected_ratioBR_gam2_test)



##### RDA Famd #####
rda_axis = rda(Y~., as.data.frame(env_coord_FAMD))
summary(rda_axis)
anova.cca(rda_axis, permutations = 10000,by ="term")
ordiplot(rda_axis,scaling = 1)

df_glm%>%
  pivot_longer(-c(log_biomass,log_plant_richness), values_to = "value_env", names_to = "axis_env")%>%
  group_by(axis_env)%>%
  mutate(value_env = exp(decostand(value_env, method = "range")+1))%>%
  ggplot()+
  facet_wrap(~axis_env)+
  geom_point(aes(log_biomass, log_plant_richness, size =value_env, col =value_env), alpha = 0.3)+
  scale_colour_continuous(type = "viridis")+
  main_theme

##### RDA selected env #####

FAMD_env_selected <- FAMD(X_selected, graph = F)
plot(FAMD_env_selected)
CAH_env_FAMD = HCPC(FAMD_env_selected, metric = "euclidean", method = "ward", nb.clust =5)

str(X_selected%>%select(!depth_oxy))
ggpairs(cbind(Y,X_selected%>%select(!depth_oxy)),
        lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.1), combo = "box_no_facet", discrete = "count", na = "na"),
        upper = list(continuous = wrap("smooth", alpha = 0.3, size=0.1), combo = "box_no_facet", discrete = "count", na = "na")
        #,ggplot2::aes(colour=CAH_env_FAMD$data.clust$clust)
)

rda_selected = rda(Y~. , as.data.frame(X_selected))

rda_summary = summary(rda_selected)

anova.cca(rda_selected, permutations = 10000,by ="term")
ordiplot(rda_selected,scaling = 1)

rda_ordistep <- ordiR2step(rda(Y ~ 1, as.data.frame(X_selected)),
                           scope = formula(rda_selected),
                           direction = "both",
                           R2scope = F,
                           pstep = 10000,
                           trace = FALSE)
anova.cca(rda_ordistep, permutations = 10000,by ="term")
ordiplot(rda_ordistep,scaling = 2)
c=100 ; b = 2; a = -20
x= seq(-2,2, length.out = 100)
plot(x, c + b*(x) + a*(x)**2 )
##### VGAM

selected_richn = glm(Y$log_plant_richness~., data = X_selected, family =gaussian(link = "identity") )
plot(selected_richn)
hist(selected_richn$residuals)
summary(selected_richn)
ggpairs(df_glm #, ggplot2::aes(colour=cluster3_quad_df$clust)
)


library(VGAM)
test_vglm = vglm(as.matrix(Y)~., data = X_selected, family = uninormal())

anova.vglm(test_vglm)
summary(test_vglm)
drop1.vglm(test_vglm)
preds_vglm = predictvglm(test_vglm)
fit = fittedvlm(test_vglm)
plot(X_selected$MO, rnorm(preds_vglm[,1],exp(exp(preds_vglm[,2]))))
points(X_selected$MO, Y[,1], col = "red")
constraints(test_vglm)
model.matrix(test_vglm, type = "lm")  # LM model matrix
model.matrix(test_vglm) 
preds <-predict(test_vglm)
head(preds_vglm)
head(preds)
plot(as.matrix(Y), bg= "black", alpha = 0.2, pch =21, cex =1)
points(rnorm(preds[,1],exp(preds[,2])) , rnorm(preds[,3],exp(preds[,4])),bg ="blue", pch =21, cex =1)
preds_data = data.frame(rnorm(preds[,1],exp(preds[,2])) , rnorm(preds[,3],exp(preds[,4])))
error = (as.matrix(Y)[,1] - preds_data[,1]) + (as.matrix(Y)[,2] - preds_data[,2])

