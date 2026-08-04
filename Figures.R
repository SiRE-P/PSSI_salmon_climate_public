###Main script to process model fit data and produce all manuscript figures###

library(ncdf4)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidybayes)
library(GGally)
library(ggpubr)
library(ggrepel)
library(scales)

#Model fit####
fit <- readRDS(file="./output/sockeye_climate_model_full.rds")

#in levelx reduce Babine Lake (spawning channel) to just Babine Lake
levelx=c("Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Chilko Lake","Chilliwack Lake","Francois Lake","Fraser Lake","Shuswap Lake","Quesnel Lake","Babine Lake","Tahltan Lake","Tatsamenie Lake")
levelxs=c("Osoyoos","Wenatchee","Great Central","Sproat","Chilko","Chilliwack","Francois","Fraser","Shuswap","Quesnel","Babine","Tahltan","Tatsamenie")
mcol = c("Osoyoos Lake"="hotpink","Wenatchee Lake"="pink2","Great Central Lake"="springgreen4","Sproat Lake"="springgreen2","Chilko Lake"="skyblue4","Chilliwack Lake"="skyblue1","Francois Lake"="blue1","Fraser Lake"="blue4","Shuswap Lake"="cyan4","Quesnel Lake"="cyan2","Babine Lake (spawning channel)"="gold1","Tahltan Lake"="purple4","Tatsamenie Lake"="orange3")
mcolx = c("Osoyoos Lake"="hotpink","Wenatchee Lake"="pink2","Great Central Lake"="springgreen4","Sproat Lake"="springgreen2","Chilko Lake"="skyblue4","Chilliwack Lake"="skyblue1","Francois Lake"="blue1","Fraser Lake"="blue4","Shuswap Lake"="cyan4","Quesnel Lake"="cyan2","Babine Lake"="gold1","Tahltan Lake"="purple4","Tatsamenie Lake"="orange3")
cus <- c("Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Chilko Lake","Chilliwack Lake","Francois Lake","Fraser Lake","Shuswap Lake","Quesnel Lake","Babine Lake","Tahltan Lake","Tatsamenie Lake")
cu_reg <- data.frame(cu_name=c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake","Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Babine Lake","Babine Lake (spawning channel)","Tahltan Lake","Tahltan Lake (enhanced)","Tatsamenie Lake","Tatsamenie Lake (enhanced)","Chutine Lake"), major_watershed=c("Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Columbia","Columbia","Somass","Somass","Skeena","Skeena","Stikine","Stikine","Taku","Taku","Taku"))
cuname_select <- c("Chilko Lake", "Great Central Lake","Sproat Lake","Tahltan Lake","Tatsamenie Lake","Quesnel Lake","Osoyoos Lake","Babine Lake (spawning channel)","Fraser Lake","Chilliwack Lake","Shuswap Lake","Francois Lake","Wenatchee Lake")
domainx<-c("South","Fraser","North")
names(domainx)<-c(1,2,3)

cuid_datx <- cuid_dat
cuid_datx$cu_name <- as.character(cuid_datx$cu_name)
cuid_datx$cu_name[cuid_datx$cu_name=="Babine Lake (spawning channel)"] = "Babine Lake"
cuid_datx$cu_name <- factor(cuid_datx$cu_name, levels=levelx)


#Map####

#coastline
coast = read_sf("./ne_10m_ocean.shp")
coast2 = st_transform(coast, crs="EPSG:4326")

#invert coastline from marine to land
sf::sf_use_s2(FALSE)
coast2 <- coast2 |> 
  st_make_valid() |>
  st_union() |>
  st_make_valid()

#Create a proper world polygon (dateline-safe)
world <- st_polygon(list(rbind(
  c(-180, -90),
  c(-180,  90),
  c( 180,  90),
  c( 180, -90),
  c(-180, -90)
))) |> 
  st_sfc(crs = 4326)

#Land = world minus ocean polygon
land <- st_difference(world, coast2)

#shelf data at 500m
shelf = read_sf("./shelf_500m.shp")

#convert
#library(terra)?
shelf <- shelf%>%mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2],
    shallow = ifelse(bathy > -500, 1, NA)  # adjust sign if needed
  )

#marine polygons
polsf = read_sf("./goa_polygons_borders2.shp")
polsf$Region = factor(polsf$Region, levels=c('Taku','Stikine','Skeena','Fraser','Somass','Columbia'))
polsf = st_transform(polsf, crs="EPSG:4326")

#freshwater polygons
#Lakes
lw3 = read_sf("./lakes_polygons_cu_plus.shp")
lw3 = st_transform(lw3, crs="EPSG:4326")
lw3$Region = sub(" River","",lw3$watershed)
lw3 <- lw3[which(lw3$cu%in%cus),]
lw3 <- st_make_valid(lw3)

#Select CUs
lw4 <- lw3[which(lw3$cu%in%cus),]

#centerpoints for labels
lw5 <- lw4[which(lw4$GNSNM1%in%c("Osoyoos Lake","Great Central Lake","Sproat Lake","Chilko Lake","Chilliwack Lake","François Lake","Fraser Lake","Shuswap Lake","Quesnel Lake","Babine Lake","Tahltan Lake","Tatsamenie Lake")),]
#lw5$cp <- st_point_on_surface(lw5$geometry)
coords <- st_coordinates(st_centroid(lw5))
lw5 <- as.data.frame(cbind(lw5,coords))
#add Wenatchee centerpoints for label
lw5 <- lw5%>%select(cu,X,Y)%>%add_row(cu="Wenatchee Lake",X=-120.7852,Y=47.8270)
lw5$cux <- sub('Lake','L.',lw5$cu)
lw5$cu <- factor(lw5$cu,levels=levelx)
lw5 <- lw5[order(lw5$cu),]
lw5 <- lw5%>%left_join(cuid_datx,join_by(cu==cu_name))

#Rivers
rw3 = read_sf("./rivers_polygons_cu_plus.shp")
rw3 = st_transform(rw3, crs="EPSG:4326")
rw3 <- rw3[which(rw3$cu%in%cus),]
rw3$Region = sub(" River","",rw3$watershed)
rw3 <- rw3 %>% group_by(cu,watershed,Region) %>% summarise(geometry = st_union(geometry), .groups = "drop")

#plot with selected marine Polygons
ggplot()+
  #geom_sf(data=shelf, colour='lightgrey')+
  geom_raster(data=shelf|>filter(shallow==1), aes(x=x,y=y),fill='lightgrey')+
  geom_sf(data=land, fill="white", colour="lightgrey")+
  geom_sf(data=polsf[which(polsf$Region%in%c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")),],linetype="solid",linewidth=1.5,show.legend=FALSE,aes(fill=Region,colour=Region))+
  #geom_sf(data=polsf[which(polsf$Region=="QCS"),],fill=NA,linetype="dotted",linewidth=1.5,show.legend=FALSE,colour='goldenrod3')+
  #geom_sf(data=polsf[which(polsf$Region=="Hecate"),],fill=NA,linetype="dotted",linewidth=1.5,show.legend=FALSE,colour='seagreen3')+
  #geom_sf(data=polsf[which(polsf$Region=="Kodiak"),],fill=NA,linetype="dotted",linewidth=1.5,show.legend=FALSE,colour='seagreen')+
  scale_fill_manual(values=c("Columbia"="darkorchid1", "Somass"="darkorchid4", "Fraser"="goldenrod3", "Skeena"="seagreen1", "Stikine"="seagreen3", "Taku"="seagreen4"))+
  geom_sf(data=lw4,fill='blue',linewidth=0.3,show.legend="point",colour='blue')+ #geom_point(aes(x=-120.7852,y=47.8270),show.legend=FALSE,shape=16,size=1.5,colour='blue')+
  geom_sf(data=rw3,fill=NA,linewidth=1,show.legend="point",size=3,aes(colour=Region))+ #geom_sf(data=we,fill=NA,linewidth=1,show.legend=FALSE,size=3,colour='darkorchid1')+
  #scale_color_manual(values=c("Columbia"="hotpink", "Somass"="springgreen4", "QCS"="red", "Hecate"="maroon4", "Fraser"="skyblue3", "Skeena"="gold1", "Stikine"="purple4", "Taku"="orange3", "Kodiak"="orchid4", "OpenOcean"="navy"))+
  scale_color_manual(values=c("Columbia"="darkorchid1", "Somass"="darkorchid4", "Fraser"="goldenrod3", "Skeena"="seagreen1", "Stikine"="seagreen3", "Taku"="seagreen4"))+
  coord_sf(xlim=c(-155,-118),ylim=c(45,60),expand=FALSE)+
  #geom_label_repel(data=lw5, size=3, nudge_x=-2, direction="both", label.size = 0, box.padding=0.1, aes(label=cux,x=X,y=Y))+
  #geom_text_repel(data=lw5, size=4, nudge_x=0.5,nudge_y=0.5, direction="both", aes(label=cu_id,x=X,y=Y))+
  geom_text_repel(data=lw5, size=4, nudge_x=2,nudge_y=1, direction="both", aes(label=pop_name_short,x=X,y=Y))+
  labs(title=NULL,colour="Watersheds",fill="Watersheds",x=NULL,y=NULL)+theme_classic(base_size=14)+theme(legend.position=c(0.1,0.3)) 
  
#Env Data####
env_dat <- read.csv("./data/env_data.csv")

#relabel "Babine Lake" to "... (spawning channel)" to match smo_dat naming
#env_dat$cu[env_dat$cu=="Babine Lake"] = "Babine Lake (spawning channel)"

env_select2 <- env_dat |>
  filter(cu %in% cus) |>
  select(year, cu_name = cu,
         FTsre = FTsreMean,
         FTwre = FTwreMean, ,
         #TEMPcs2 = TEMPcsMean2,
         #TEMPcs = TEMPcsMean,
         #MLDcs2 = MLDcsMean2,
         #MLDcs = MLDcsMean,
         TEMPsoo = TEMPsooMean,
         TEMPwoo = TEMPwooMean,
         SSTcs = SSTcsMean,
         SSTsoo = SSTsooMean,
         SSTwoo = SSTwooMean,
         FTum = FTumMean,
         FDum = FDumMean) |>
  mutate(cu_name = factor(cu_name, levels = levelx)) |>
  mutate(major_watershed = case_when(cu_name %in% c("Osoyoos Lake","Wenatchee Lake") ~ "Columbia",
                                     cu_name %in% c("Great Central Lake","Sproat Lake") ~ "Somass",
                                     cu_name %in% c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake") ~ "Fraser",
                                     cu_name == "Babine Lake" ~ "Skeena",
                                     cu_name == "Tahltan Lake" ~ "Stikine",
                                     cu_name == "Tatsamenie Lake" ~ "Taku")) |>
  mutate(major_watershed = factor(major_watershed, levels = c("Columbia","Somass","Fraser","Skeena","Stikine","Taku"))) |>
  mutate(domain = case_when(major_watershed %in% c("Columbia", "Somass") ~ "South",
                            major_watershed == "Fraser" ~ "Fraser",
                            major_watershed  %in% c("Stikine", "Taku", "Skeena") ~ "North")) |>
  mutate(domain = factor(domain, levels = c("South", "Fraser", "North")))

#Set to NA 2024 winter rearing temps for all as they lack the jan-march 2025 data to integrate
env_select2[which(env_select2$year==2024),'FTwre'] = NA
#env_select[which(env_select$year==2024),'GCMTwre'] = NA

#Set to NA summer rearing temps for Osoyoos and Wenatchee lake 2024 because they are based on unreliable observation data and seem unreasonably high
env_select2[which(env_select2$cu_name%in%c("Osoyoos Lake","Wenatchee Lake")&env_select2$year==2024),'FTsre'] = NA
#env_select[which(env_select$cu_name%in%c("Osoyoos Lake","Wenatchee Lake")&env_select$year==2024),'GCMTsre'] = NA

#Coastal Data
coastal_dat2 <- read.csv("./data/env_data_region.csv")
coastal_dat2 <- coastal_dat2 |> filter(!region=="OpenOcean")|>
  mutate(region = factor(region, levels=c("Columbia","Fraser","Somass","QCS","Hecate","Skeena","Stikine","Taku","Kodiak","OpenOcean")))|>
  mutate(domain = case_when(region %in% c("Columbia", "Somass") ~ "South",
                            region == "Fraser" ~ "Fraser",
                            region %in% c("QCS","Hecate") ~ "Hecate",
                            region  %in% c("Stikine", "Taku", "Skeena") ~ "North",
                            region == "Kodiak" ~ "Kodiak")) |>
  mutate(domain = factor(domain, levels = c("South","Fraser","Hecate","North","Kodiak")))

#Raw data overview
# summary(coastal_dat2)
# summary(coastal_dat2%>%filter(domain=="South"))
# summary(coastal_dat2%>%filter(domain=="Fraser"))
# summary(coastal_dat2%>%filter(domain=="North"))

env_long <- env_select2 |>
  gather(key = variable, value = value, -cu_name, -year, -major_watershed, -domain) |>
  group_by(variable) |>
  mutate(value_stnd = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |> #Global stdn
  group_by(variable, cu_name, major_watershed) |>
  mutate(value_stnd_cu = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |> #CU stdn
  ungroup()
 
# env_long2 <- env_long |> mutate(variable=factor(variable, levels=c("FTwre","FTsre","TEMPcs","MLDcs","TEMPwo","TEMPso","SSTwo","SSTso","FTum","FDum"))) |>
#   filter(year<2024)
# 
# env_long2$cu_namex <- factor(str_remove(as.character(env_long2$cu_name)," Lake"),levels=cuid_dat$pop_name_short)
# env_long2$variablex <- as.character(env_long2$variable)
# env_long2 <- env_long2%>%mutate(variablex=recode(variablex,
#                                                  'FTwre'='Temp. Winter Rearing',
#                                                  'FTsre'='Temp. Summer Rearing',
#                                                  'TEMPcs'='Temp. Coast',
#                                                  'MLDcs'='MLD Coast',
#                                                  'TEMPwo'='Temp. Winter Open Ocean',
#                                                  'TEMPso'='Temp. Summer Open Ocean',
#                                                  'SSTwo'='SST Winter Open Ocean',
#                                                  'SSTso'='SST Summer Open Ocean',
#                                                  'FTum'='Temp. Return Migration',
#                                                  'FDum'='Discharge Return Migration'))%>%
#   mutate(variablex=factor(variablex,levels=c('Temp. Winter Rearing','Temp. Summer Rearing','Temp. Coast','MLD Coast','Temp. Winter Open Ocean','Temp. Summer Open Ocean','SST Winter Open Ocean','SST Summer Open Ocean','Temp. Return Migration','Discharge Return Migration')))
# 
# ggplot(data=env_long2,aes(x = year, y = cu_namex, fill = value_stnd_cu))+
#   geom_tile()+
#   scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="lightgrey")+
#   facet_wrap(~variablex,ncol=2)+
#   labs(x = "", y = "",fill="Scaled Value")+
#   scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
#   theme_bw(base_size = 10)+
#   theme(legend.position='bottom')+
#   coord_cartesian(expand = FALSE)

#Posterior env. data
library(patchwork)
major_watershed<-unique(cuid_datx$major_watershed)
domain<-unique(cuid_datx$domain)

FTwre<-spread_draws(fit, FTwre_complete[y,p])%>%group_by(y,p)%>%summarize(value=mean(FTwre_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs))%>%ungroup()%>%group_by(pop)%>%mutate(value_scaled=(value-mean(value,na.rm=TRUE))/sd(value,na.rm=TRUE))%>%select(pop, year, value,value_scaled)
FTsre<-spread_draws(fit, FTsre_complete[y,p])%>%group_by(y,p)%>%summarize(value=mean(FTsre_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs))%>%ungroup()%>%group_by(pop)%>%mutate(value_scaled=(value-mean(value,na.rm=TRUE))/sd(value,na.rm=TRUE))%>%select(pop, year, value,value_scaled)
FTum<-spread_draws(fit, FTum_complete[y,p])%>%group_by(y,p)%>%summarize(value=mean(FTum_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs))%>%ungroup()%>%select(pop, year, value)
FDum<-spread_draws(fit, FDum_complete[y,p])%>%group_by(y,p)%>%summarize(value=mean(FDum_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs))%>%ungroup()%>%select(pop, year, value)
#TEMPcs<-spread_draws(fit, TEMPcs_complete[y,r])%>%group_by(y,r)%>%summarize(value=mean(TEMPcs_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),watershed=factor(major_watershed[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(watershed, year, value)
#MLDcs<-spread_draws(fit, MLDcs_complete[y,r])%>%group_by(y,r)%>%summarize(value=mean(MLDcs_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),watershed=factor(major_watershed[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(watershed, year, value)
TEMPcs<-spread_draws(fit, TEMPstnd[y,r])%>%group_by(y,r)%>%summarize(value=mean(TEMPstnd,na.rm=TRUE))%>%mutate(year=y+(year_start-1),watershed=factor(major_watershed[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(watershed, year, value)
MLDcs<-spread_draws(fit, MLDstnd[y,r])%>%group_by(y,r)%>%summarize(value=mean(MLDstnd,na.rm=TRUE))%>%mutate(year=y+(year_start-1),watershed=factor(major_watershed[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(watershed, year, value)
TEMPwoo<-spread_draws(fit, TEMPwoo_complete[y])%>%group_by(y)%>%summarize(value=mean(TEMPwoo_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),openocean="Open Ocean")%>%ungroup()%>%select(openocean, year, value)
TEMPsoo<-spread_draws(fit, TEMPsoo_complete[y])%>%group_by(y)%>%summarize(value=mean(TEMPsoo_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),openocean="Open Ocean")%>%ungroup()%>%select(openocean, year, value)

##Add projection data, and empty space to seperate historic and projection
FTwrep <- read.csv("./data/FTwre_projected245_2041_2070.csv")
FTwrep <- FTwrep|>select(-X) |> gather(key = pop, value = value, Babine.Lake:Wenatchee.Lake) |> 
  #mutate(value_scaled = (value - mean(env_select2$FTwre, na.rm = TRUE)) / sd(env_select2$FTwre, na.rm = TRUE)) |>
  mutate(pop = str_trim(str_remove(pop, "\\b.[Ll]ake\\b"))) |> mutate(pop = ifelse(pop == "Great.Central", "Great Central", pop)) |>
  group_by(pop) |> summarize(value = mean(value, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), pop = factor(pop,levels=levelxs)) 
FTwre2 <- rbind(FTwre|>select(pop,value,year)|> mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(pop=FTwrep$pop,value=NA,year="2025-2040"),FTwrep) |> mutate(value_scaled=(value-mean(value,na.rm=TRUE))/sd(value,na.rm=TRUE)) |> select(pop, year, value,value_scaled)

FTsrep <- read.csv("./data/FTsre_projected245_2041_2070.csv")
FTsrep <- FTsrep|>select(-X) |> gather(key = pop, value = value, Babine.Lake:Wenatchee.Lake) |> 
  mutate(pop = str_trim(str_remove(pop, "\\b.[Ll]ake\\b"))) |> mutate(pop = ifelse(pop == "Great.Central", "Great Central", pop)) |>
  group_by(pop) |> summarize(value = mean(value, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), pop = factor(pop,levels=levelxs)) 
FTsre2 <- rbind(FTsre|>select(pop,value,year)|>mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(pop=FTsrep$pop,value=NA,year="2025-2040"),FTsrep) |> mutate(value_scaled=(value-mean(value,na.rm=TRUE))/sd(value,na.rm=TRUE)) |> select(pop, year, value,value_scaled)

TEMPcsp <- read.csv("./data/TEMPcs_projected245_2041_2070.csv")
TEMPcsp <- TEMPcsp|>select(-X) |> gather(key = watershed, value = value, Fraser:Taku) |> 
  group_by(watershed) |> summarize(value = mean(value, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), watershed = factor(watershed,levels=major_watershed)) 
TEMPcs2 <- rbind(TEMPcs|>select(watershed,value,year)|>mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(watershed=TEMPcsp$watershed,value=NA,year="2025-2040"),TEMPcsp)

MLDcsp <- read.csv("./data/MLDcs_projected245_2041_2070.csv")
MLDcsp <- MLDcsp|>select(-X) |> gather(key = watershed, value = value, Fraser:Taku) |> 
  group_by(watershed) |> summarize(value = mean(value, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), watershed = factor(watershed,levels=major_watershed)) 
MLDcs2 <- rbind(MLDcs|>select(watershed,value,year)|>mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(watershed=MLDcsp$watershed,value=NA,year="2025-2040"),MLDcsp)

TEMPoop <- read.csv("./data/TEMP_OpenOcean_projected245_2041_2070.csv")
TEMPwoop <- TEMPoop|>select(-year, TEMPwoo) |> summarize(value = mean(TEMPwoo, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), openocean="Open Ocean") 
TEMPwoo2 <- rbind(TEMPwoo|>select(value,year,openocean)|>mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(openocean=TEMPwoop$openocean,value=NA,year="2025-2040"),TEMPwoop)

TEMPsoop <- TEMPoop|>select(-year, TEMPsoo) |> summarize(value = mean(TEMPsoo, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), openocean="Open Ocean") 
TEMPsoo2 <- rbind(TEMPsoo|>select(value,year,openocean)|>mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(openocean=TEMPsoop$openocean,value=NA,year="2025-2040"),TEMPsoop)

FTump <- read.csv("./data/FTum_projected245_2041_2070.csv")
FTump <- FTump|>select(-X) |> gather(key = pop, value = value, Osoyoos.Lake:Tatsamenie.Lake) |> 
  mutate(pop = str_trim(str_remove(pop, "\\b.[Ll]ake\\b"))) |> mutate(pop = ifelse(pop == "Great.Central", "Great Central", pop)) |>
  group_by(pop) |> summarize(value = mean(value, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), pop = factor(pop,levels=levelxs)) 
FTum2 <- rbind(FTum|>select(pop,value,year)|>mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(pop=FTump$pop,value=NA,year="2025-2040"),FTump)

FDump <- read.csv("./data/FDum_projected245_2041_2070.csv")
FDump <- FDump|>select(-X) |> gather(key = pop, value = value, Osoyoos.Lake:Tatsamenie.Lake) |> 
  mutate(pop = str_trim(str_remove(pop, "\\b.[Ll]ake\\b"))) |> mutate(pop = ifelse(pop == "Great.Central", "Great Central", pop)) |>
  group_by(pop) |> summarize(value = mean(value, na.rm = TRUE)) |> mutate(year = factor("2041-2070"), pop = factor(pop,levels=levelxs)) 
FDum2 <- rbind(FDum|>select(pop,value,year)|>mutate(year=factor(year, levels=as.character(1981:2024))),data.frame(pop=FDump$pop,value=NA,year="2025-2040"),FDump)

ggplot(data=FTwre2,aes(x = year, y = pop, fill = value_scaled))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Winter Freshwater Temperature",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE) +
  
  ggplot(data=FTsre2,aes(x = year, y = pop, fill = value_scaled))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Summer Freshwater Temperature",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE)+
  
  ggplot(data=TEMPcs2,aes(x = year, y = watershed, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Coastal Temperature",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE)+
  
  ggplot(data=MLDcs2,aes(x = year, y = watershed, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Coastal Mixed Layer Depth",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE)+
  
  ggplot(data=TEMPwoo2,aes(x = year, y = openocean, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Winter Open Ocean Temperature",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE)+
  
  ggplot(data=TEMPsoo2,aes(x = year, y = openocean, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Summer Open Ocean Temperature",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE)+
  
  ggplot(data=FTum2,aes(x = year, y = pop, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Return Migration Temperature",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE)+
  
  ggplot(data=FDum2,aes(x = year, y = pop, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue", high = "red",na.value="white")+
  labs(title="Return Migration Discharge",x=NULL, y=NULL,fill="Scaled Value")+
  scale_x_discrete(breaks=c("1985","1990","2000","2010","2020","2041-2070"))+
  theme_bw(base_size = 12)+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  geom_vline(xintercept = "2025-2040",linetype=2,lwd=1)+
  coord_cartesian(expand = FALSE)+
  
  plot_layout(ncol=2, axis_titles='collect_x',axes='collect_y',guides='collect')&theme(legend.position='none')

##Freshwater imputations as Supplementary Figures
mcoly = c("Osoyoos"="hotpink","Wenatchee"="pink2","Great Central"="springgreen4","Sproat"="springgreen2","Chilko"="skyblue4","Chilliwack"="skyblue1","Francois"="blue1","Fraser"="blue4","Shuswap"="cyan4","Quesnel"="cyan2","Babine"="gold1","Tahltan"="purple4","Tatsamenie"="orange3")

spread_draws(fit, FTwre_complete[y,p])|>
  mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs)) |>
  mutate(domain = case_when(pop %in% c("Great Central","Sproat", "Osoyoos", "Wenatchee") ~ "South",
                            pop %in% c("Fraser","Francois","Chilko","Chilliwack","Cultus","Shuswap","Quesnel") ~ "Fraser", 
                            pop  %in% c("Babine", "Tahltan", "Tatsamenie") ~ "North")) |>
  mutate(FTwre = FTwre_complete * sd(env_select2$FTwre, na.rm = TRUE) + mean(env_select2$FTwre, na.rm = TRUE)) |>
  ggplot(aes(x = year, y = FTwre, colour=pop))+
  scale_colour_manual(values=mcoly)+
  #ggplot(aes(x = year, y = FTwre, colour=domain))+
  #scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3'))+
  stat_lineribbon(.width=0.5,fill="lightgrey")+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Temp. Winter Rearing",x=NULL)+
  theme(legend.position='bottom') #+ facet_wrap(~pop, scales = "free_y")

spread_draws(fit, FTsre_complete[y,p])|>
  mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs)) |>
  mutate(domain = case_when(pop %in% c("Great Central","Sproat", "Osoyoos", "Wenatchee") ~ "South",
                            pop %in% c("Fraser","Francois","Chilko","Chilliwack","Cultus","Shuswap","Quesnel") ~ "Fraser", 
                            pop  %in% c("Babine", "Tahltan", "Tatsamenie") ~ "North")) |>
  mutate(FTsre = FTsre_complete * sd(env_select2$FTsre, na.rm = TRUE) + mean(env_select2$FTsre, na.rm = TRUE)) |>
  ggplot(aes(x = year, y = FTsre, colour=pop))+
  scale_colour_manual(values=mcoly)+
  #ggplot(aes(x = year, y = FTsre, colour=domain))+
  #scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3'))+
  stat_lineribbon(.width=0.5,fill="lightgrey")+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Temp. Summer Rearing",x=NULL)+
  theme(legend.position='bottom') #+facet_wrap(~pop, scales = "free_y")

FTum_unscaled <- env_select2 |> mutate(pop = factor(str_trim(str_remove(cu_name, "\\b[Ll]ake\\b")),levels=levelxs)) |> select(year, pop, FTum_unscaled=FTum) 
spread_draws(fit, FTum_complete[y,p])|>
  mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs)) |>
  mutate(domain = case_when(pop %in% c("Great Central","Sproat", "Osoyoos", "Wenatchee") ~ "South",
                            pop %in% c("Fraser","Francois","Chilko","Chilliwack","Cultus","Shuswap","Quesnel") ~ "Fraser", 
                            pop  %in% c("Babine", "Tahltan", "Tatsamenie") ~ "North")) |>
  left_join(FTum_unscaled) |>
  group_by(pop)|>
  mutate(FTum = FTum_complete * sd(FTum_unscaled, na.rm = TRUE) + mean(FTum_unscaled, na.rm = TRUE)) |> ungroup() |>
  ggplot(aes(x = year, y = FTum, colour=pop))+
  scale_colour_manual(values=mcoly)+
  #ggplot(aes(x = year, y = FTum, colour=domain))+
  #scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3'))+
  stat_lineribbon(.width=0.5,fill="lightgrey")+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Temp. Return Migration",x=NULL)+
  theme(legend.position='bottom') #+facet_wrap(~pop, scales = "free_y")

FDum_unscaled <- env_select2 |> mutate(pop = factor(str_trim(str_remove(cu_name, "\\b[Ll]ake\\b")),levels=levelxs)) |> select(year, pop, FDum_unscaled=FDum) 
spread_draws(fit, FDum_complete[y,p])|>
  mutate(year=y+(year_start-1),pop=factor(cuid_datx$pop_name_short[p],levels=levelxs)) |>
  mutate(domain = case_when(pop %in% c("Great Central","Sproat", "Osoyoos", "Wenatchee") ~ "South",
                            pop %in% c("Fraser","Francois","Chilko","Chilliwack","Cultus","Shuswap","Quesnel") ~ "Fraser", 
                            pop  %in% c("Babine", "Tahltan", "Tatsamenie") ~ "North")) |>
  left_join(FDum_unscaled) |>
  group_by(pop)|>
  mutate(FDum = FDum_complete * sd(FDum_unscaled, na.rm = TRUE) + mean(FDum_unscaled, na.rm = TRUE)) |> ungroup() |>
  ggplot(aes(x = year, y = FDum, colour=pop))+
  scale_colour_manual(values=mcoly)+
  #ggplot(aes(x = year, y = FDum, colour=domain))+
  #scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3'))+
  stat_lineribbon(.width=0.5,fill="lightgrey")+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Disc. Return Migration",x=NULL)+
  theme(legend.position='bottom') #+facet_wrap(~pop, scales = "free_y")

##Marine imputations as Supplementary Figures
regionx = c("Columbia","Fraser","Somass","QCS","Hecate","Skeena","Stikine","Taku","Kodiak","OpenOcean")

TEMPcs_unscaled <- coastal_dat2 |> mutate(region = factor(region,levels=c("Columbia","Fraser","Somass","QCS","Hecate","Skeena","Stikine","Taku","Kodiak"))) |> select(year, region, TEMPcs_unscaled=TEMP) 
spread_draws(fit, TEMP_complete[y,r])|>
  mutate(year=y+(year_start-1), region = factor(regionx[r],levels=c("Columbia","Fraser","Somass","QCS","Hecate","Skeena","Stikine","Taku","Kodiak","OpenOcean"))) |>
  mutate(TEMP = TEMP_complete * sd(TEMPcs_unscaled$TEMPcs_unscaled, na.rm = TRUE) + mean(TEMPcs_unscaled$TEMPcs_unscaled, na.rm = TRUE)) |> 
  ggplot(aes(x = year, y = TEMP, colour=region))+
  scale_colour_manual(values=c("Columbia"="hotpink", "Somass"="springgreen4", "QCS"="red", "Hecate"="maroon4", "Fraser"="skyblue3", "Skeena"="gold1", "Stikine"="purple4", "Taku"="orange3", "Kodiak"="orchid4"))+
  stat_lineribbon(.width=0.5,fill="lightgrey")+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Coastal Temperature",x=NULL)+
  theme(legend.position='bottom') 

MLDcs_unscaled <- coastal_dat2 |> mutate(region = factor(region,levels=c("Columbia","Fraser","Somass","QCS","Hecate","Skeena","Stikine","Taku","Kodiak"))) |> select(year, region, MLDcs_unscaled=MLD) 
spread_draws(fit, MLD_complete[y,r])|>
  mutate(year=y+(year_start-1), region = factor(regionx[r],levels=c("Columbia","Fraser","Somass","QCS","Hecate","Skeena","Stikine","Taku","Kodiak","OpenOcean"))) |>
  mutate(MLD = MLD_complete * sd(MLDcs_unscaled$MLDcs_unscaled, na.rm = TRUE) + mean(MLDcs_unscaled$MLDcs_unscaled, na.rm = TRUE)) |> 
  ggplot(aes(x = year, y = MLD, colour=region))+
  scale_colour_manual(values=c("Columbia"="hotpink", "Somass"="springgreen4", "QCS"="red", "Hecate"="maroon4", "Fraser"="skyblue3", "Skeena"="gold1", "Stikine"="purple4", "Taku"="orange3", "Kodiak"="orchid4"))+
  stat_lineribbon(.width=0.5,fill="lightgrey")+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Coastal Mixed Layer Depth",x=NULL)+
  theme(legend.position='bottom') 

spread_draws(fit, TEMPwoo_complete[y])%>%mutate(year=y+(year_start-1), TEMP = TEMPwoo_complete * sd(env_select2$TEMPwoo, na.rm = TRUE) + mean(env_select2$TEMPwoo, na.rm = TRUE),season="Winter")%>%select(-TEMPwoo_complete)|>
  ggplot(aes(x = year, y = TEMP))+
  stat_lineribbon(.width=0.5,fill="lightgrey",colour='dodgerblue')+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Winter Open Ocean Temperatures",x=NULL)

spread_draws(fit, TEMPsoo_complete[y])%>%mutate(year=y+(year_start-1), TEMP = TEMPsoo_complete * sd(env_select2$TEMPsoo, na.rm = TRUE) + mean(env_select2$TEMPsoo, na.rm = TRUE),season="Summer")%>%select(-TEMPsoo_complete)|>
  ggplot(aes(x = year, y = TEMP))+
  stat_lineribbon(.width=0.5,fill="lightgrey",colour='red')+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Summer Open Ocean Temperatures",x=NULL)

#Predicted_vs_Observed####

#Adult Observed
se_adult <- read.csv("./data/adult_table_public.csv", header = TRUE)
adult_data <- se_adult%>%filter(!en_route_mort=="included") |>
  filter(cu_name %in% cuname_select, adult_return_year>=1981) |>
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |>
  select(cu_name, year=adult_return_year, spawners=total_escape, returns_1.1=total_returns_1.1, returns_1.2=total_returns_1.2, returns_1.3=total_returns_1.3, returns_2.1=total_returns_2.1, returns_2.2=total_returns_2.2, returns_2.3=total_returns_2.3)|>
  mutate(returns_total=rowSums(across(starts_with('returns_')),na.rm=TRUE))
adult_data$cu_name <- as.character(adult_data$cu_name)
adult_data$cu_name[adult_data$cu=="Babine Lake (spawning channel)"] = "Babine Lake"
adult_data$cu_name <- factor(adult_data$cu_name, levels=levelx)
adult_data$returns_total[which(adult_data$returns_total==0)] = NA
adult_data$domain <- "observed"

#Predicted
predictions <- spread_draws(fit, returns_age[y, fw_age, mar_age, cu_id]) |> 
  group_by(y,cu_id,.chain,.iteration,.draw)|>
  summarise(returns_total=sum(returns_age))|> 
  left_join(spread_draws(fit, spawners[y,cu_id]))|>
  ungroup()|>
  mutate(year=y+(year_start-1))|>
  left_join(cuid_datx)|>
  mutate(cu_name=factor(cu_name,levels=c("Tatsamenie Lake","Tahltan Lake","Babine Lake","Quesnel Lake","Shuswap Lake","Fraser Lake","Francois Lake","Chilliwack Lake","Chilko Lake","Sproat Lake","Great Central Lake","Wenatchee Lake","Osoyoos Lake")))

ggplot()+
  stat_pointinterval(data=predictions,aes(x = year, y = returns_total,colour=domain),shape=19,.width=0.66,linewidth=1,size=7)+
  #stat_pointinterval(data=predictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = returns_total,colour=domain),shape=19,.width=0.66,linewidth=2,size=9)+
  #stat_pointinterval(data=predictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = returns_total,colour=domain),shape=NA,.width=0.66,linewidth=NA,size=9)+
  scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3','observed'='grey29'))+
  geom_point(data=adult_data,aes(x=year,y=returns_total,colour=domain),shape=18,size=3)+
  #geom_point(data=adult_data%>%filter(cu_name=="Sproat Lake"),aes(x=year,y=returns_total,colour=domain),shape=18,size=5)+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme_bw(base_size=10)+
  #theme_bw(base_size=14)+
  labs(title=NULL,x=NULL,y="Returns",colour=NULL)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme(legend.position='bottom')
  #theme(legend.position='none')

ggplot()+
  stat_pointinterval(data=predictions,aes(x = year, y = spawners,colour=domain),shape=19,.width=0.66,linewidth=1,size=7)+
  #stat_pointinterval(data=predictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = spawners,colour=domain),shape=19,.width=0.66,linewidth=2,size=9)+
  #stat_pointinterval(data=predictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = spawners,colour=domain),shape=NA,.width=0.66,linewidth=NA,size=9)+
  scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3','observed'='grey29'))+
  geom_point(data=adult_data,aes(x=year,y=spawners,colour=domain),shape=18,size=3)+
  #geom_point(data=adult_data%>%filter(cu_name=="Sproat Lake"),aes(x=year,y=spawners,colour=domain),shape=18,size=5)+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme_bw(base_size=10)+
  #theme_bw(base_size=14)+
  labs(title=NULL,x=NULL,y="Spawners",colour=NULL)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme(legend.position='bottom')
  #theme(legend.position='none')

#Juvenile Observed
se_juv <- read.csv("./data/juvenile_table_public.csv", header = TRUE)
juv_data <- se_juv |> 
  filter(cu_name %in% cuname_select) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  select(cu_name, smolt_migration_year, smolt_abundance_age0, smolt_abundance_age1, smolt_abundance_age2, smolt_abundance_age3, smolt_abundance_total, presmolt_abundance) |> 
  filter(smolt_migration_year >= year_start) |> 
  group_by(cu_name, smolt_migration_year) |> 
  mutate(smolt_abundance_total = ifelse(is.na(smolt_abundance_total), sum(smolt_abundance_age0, smolt_abundance_age1, smolt_abundance_age2, smolt_abundance_age3, na.rm = TRUE), smolt_abundance_total)) |> 
  mutate(smolt_abundance_total = ifelse(smolt_abundance_total == 0, NA, smolt_abundance_total)) |> 
  mutate(smolt_abundance_total = ifelse(!is.na(smolt_abundance_total) &
                                          !is.na(presmolt_abundance) &
                                          smolt_abundance_total == presmolt_abundance,
                                        NA, smolt_abundance_total)) |>  
  filter(!is.na(smolt_abundance_total) | !is.na(presmolt_abundance)) |> 
  arrange(cu_name, smolt_migration_year) |> 
  rename(smolts = smolt_abundance_total, presmolts = presmolt_abundance) |> 
  mutate(smolts = smolts - rowSums(across(c(smolt_abundance_age0, smolt_abundance_age3)), na.rm = TRUE)) |> #remove age0s and 3s from totals
  select(cu_name, smolt_migration_year, smolts, presmolts) |> 
  left_join(juv_select_table) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  mutate(juv_abund = ifelse(type == "smolts", smolts, presmolts)) |> 
  select(cu_name, smolt_migration_year, juv_abund, type) |> 
  filter(!is.na(juv_abund))
juv_data$cu_name <- as.character(juv_data$cu_name)
juv_data$cu_name[juv_data$cu_name=="Babine Lake (spawning channel)"] = "Babine Lake"
juv_data$cu_name <- factor(juv_data$cu_name, levels=levelx)
juv_data$juv_abund[which(juv_data$juv_abund==0)] = NA
juv_data$domain <- "observed"

#Predicted
spredictions <- spread_draws(fit, smolts_age[y, fw_age, cu_id]) |> 
  group_by(y,cu_id,.chain,.iteration,.draw)|>
  summarise(smolts_total=sum(smolts_age))|> 
  left_join(spread_draws(fit, spawners[y,cu_id]))|>
  ungroup()|>
  mutate(year=y+(year_start-1))|>
  left_join(cuid_datx)|>
  mutate(cu_name=factor(cu_name,levels=c("Tatsamenie Lake","Tahltan Lake","Babine Lake","Quesnel Lake","Shuswap Lake","Fraser Lake","Francois Lake","Chilliwack Lake","Chilko Lake","Sproat Lake","Great Central Lake","Wenatchee Lake","Osoyoos Lake")))

ggplot()+
  stat_pointinterval(data=spredictions,aes(x = year, y = smolts_total,colour=domain),shape=19,.width=0.66,linewidth=1,size=7)+
  #stat_pointinterval(data=spredictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = smolts_total,colour=domain),shape=19,.width=0.66,linewidth=2,size=9)+
  #stat_pointinterval(data=spredictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = smolts_total,colour=domain),shape=NA,.width=0.66,linewidth=NA,size=9)+
  scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3','observed'='grey29'))+
  geom_point(data=juv_data,aes(x=smolt_migration_year,y=juv_abund,colour=domain),shape=18,size=3)+
  #geom_point(data=juv_data%>%filter(cu_name=="Sproat Lake"),aes(x=smolt_migration_year,y=juv_abund,colour=domain),shape=18,size=5)+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme_bw(base_size=10)+
  #theme_bw(base_size=14)+
  labs(title=NULL,x=NULL,y="Smolts",colour=NULL)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme(legend.position='bottom')
  #theme(legend.position='none')

#Parameters####

#mu parameters
b_FTwre_hmu<-spread_draws(fit, b_FTwre_hyper_mu)%>%mutate(d="Hyper",variable="FTwre")%>%select(.chain,.iteration,.draw,d,value=b_FTwre_hyper_mu,variable)
b_FTwre_mu<-spread_draws(fit, b_FTwre_mu[d,w])%>%mutate(d=domainx[d],variable="FTwre")%>%select(.chain,.iteration,.draw,d,value=b_FTwre_mu,variable)
b_FTsre_hmu<-spread_draws(fit, b_FTsre_hyper_mu)%>%mutate(d="Hyper",variable="FTsre")%>%select(.chain,.iteration,.draw,d,value=b_FTsre_hyper_mu,variable)
b_FTsre_mu<-spread_draws(fit, b_FTsre_mu[d,s])%>%mutate(d=domainx[d],variable="FTsre")%>%select(.chain,.iteration,.draw,d,value=b_FTsre_mu,variable)
#b_FTwre_mu<-spread_draws(fit, b_FTwre_mu[w])%>%mutate(d="Hyper",variable="FTwre")%>%select(.chain,.iteration,.draw,d,value=b_FTwre_mu,variable)
#b_FTsre_mu<-spread_draws(fit, b_FTsre_mu[s])%>%mutate(d="Hyper",variable="FTsre")%>%select(.chain,.iteration,.draw,d,value=b_FTsre_mu,variable)
b_TEMPcs_hmu<-spread_draws(fit, b_TEMPcs_hyper_mu)%>%mutate(d="Hyper",variable="TEMPcs")%>%select(.chain,.iteration,.draw,d,value=b_TEMPcs_hyper_mu,variable)
b_TEMPcs_mu<-spread_draws(fit, b_TEMPcs_mu[d])%>%mutate(d=domainx[d],variable="TEMPcs")%>%select(.chain,.iteration,.draw,d,value=b_TEMPcs_mu,variable)
b_MLDcs_hmu<-spread_draws(fit, b_MLDcs_hyper_mu)%>%mutate(d="Hyper",variable="MLDcs")%>%select(.chain,.iteration,.draw,d,value=b_MLDcs_hyper_mu,variable)
b_MLDcs_mu<-spread_draws(fit, b_MLDcs_mu[d])%>%mutate(d=domainx[d],variable="MLDcs")%>%select(.chain,.iteration,.draw,d,value=b_MLDcs_mu,variable)
b_TEMPwoo_hmu<-spread_draws(fit, b_TEMPwoo_hyper_mu)%>%mutate(d="Hyper",variable="TEMPwoo")%>%select(.chain,.iteration,.draw,d,value=b_TEMPwoo_hyper_mu,variable)
b_TEMPwoo_mu<-spread_draws(fit, b_TEMPwoo_mu[d,w])%>%mutate(d=domainx[d],variable="TEMPwoo")%>%select(.chain,.iteration,.draw,d,value=b_TEMPwoo_mu,variable)
b_TEMPsoo_hmu<-spread_draws(fit, b_TEMPsoo_hyper_mu)%>%mutate(d="Hyper",variable="TEMPsoo")%>%select(.chain,.iteration,.draw,d,value=b_TEMPsoo_hyper_mu,variable)
b_TEMPsoo_mu<-spread_draws(fit, b_TEMPsoo_mu[d,s])%>%mutate(d=domainx[d],variable="TEMPsoo")%>%select(.chain,.iteration,.draw,d,value=b_TEMPsoo_mu,variable)
b_FTum_hmu<-spread_draws(fit, b_FTum_hyper_mu)%>%mutate(d="Hyper",variable="FTum")%>%select(.chain,.iteration,.draw,d,value=b_FTum_hyper_mu,variable)
b_FTum_mu<-spread_draws(fit, b_FTum_mu[d])%>%mutate(d=domainx[d],variable="FTum")%>%select(.chain,.iteration,.draw,d,value=b_FTum_mu,variable)
b_FDum_hmu<-spread_draws(fit, b_FDum_hyper_mu)%>%mutate(d="Hyper",variable="FDum")%>%select(.chain,.iteration,.draw,d,value=b_FDum_hyper_mu,variable)
b_FDum_mu<-spread_draws(fit, b_FDum_mu[d])%>%mutate(d=domainx[d],variable="FDum")%>%select(.chain,.iteration,.draw,d,value=b_FDum_mu,variable)

mdf<-data.frame(rbind(b_FTwre_hmu[,c(".chain",".iteration",".draw","d","value","variable")],b_FTwre_mu[,c(".chain",".iteration",".draw","d","value","variable")],
                      b_FTsre_hmu[,c(".chain",".iteration",".draw","d","value","variable")],b_FTsre_mu[,c(".chain",".iteration",".draw","d","value","variable")],
                      b_TEMPcs_hmu,b_TEMPcs_mu,
                      b_MLDcs_hmu,b_MLDcs_mu,
                      b_TEMPwoo_hmu[,c(".chain",".iteration",".draw","d","value","variable")],b_TEMPwoo_mu[,c(".chain",".iteration",".draw","d","value","variable")],
                      b_TEMPsoo_hmu[,c(".chain",".iteration",".draw","d","value","variable")],b_TEMPsoo_mu[,c(".chain",".iteration",".draw","d","value","variable")],
                      b_FTum_hmu,b_FTum_mu,
                      b_FDum_hmu,b_FDum_mu))

mdf$d <- factor(mdf$d, levels=c("Hyper","South","Fraser","North"))
mdf$variable <- factor(mdf$variable, levels=c("FTwre","FTsre","TEMPcs","MLDcs","TEMPwoo","TEMPsoo","FTum","FDum"))
mdf <- mdf%>%mutate(v = case_when(variable == "FTwre" ~ "Winter Freshwater Temperature",
                                  variable == "FTsre" ~ "Summer Freshwater Temperature",
                                  variable == "TEMPcs" ~ "Coastal Temperature",
                                  variable == "MLDcs" ~ "Coastal Mixed Layer Depth",
                                  variable == "TEMPwoo" ~ "Winter Open Ocean Temperature",
                                  variable == "TEMPsoo" ~ "Summer Open Ocean Temperature",
                                  variable == "FTum" ~ "Return Migration Temperature",
                                  variable == "FDum" ~ "Return Migration Discharge"))
mdf$v = factor(mdf$v,levels=c("Winter Freshwater Temperature","Summer Freshwater Temperature","Coastal Temperature","Coastal Mixed Layer Depth","Winter Open Ocean Temperature","Summer Open Ocean Temperature","Return Migration Temperature","Return Migration Discharge"))

#pop parameters
b_FTwre<-spread_draws(fit, b_FTwre[p,w])
b_FTsre<-spread_draws(fit, b_FTsre[p,s])
b_TEMPcs<-spread_draws(fit, b_TEMPcs[p])
b_MLDcs<-spread_draws(fit, b_MLDcs[p])
b_TEMPwoo<-spread_draws(fit, b_TEMPwoo[p,w])
b_TEMPsoo<-spread_draws(fit, b_TEMPsoo[p,s])
b_FTum<-spread_draws(fit, b_FTum[p])
b_FDum<-spread_draws(fit, b_FDum[p])

b_FTwre2<-b_FTwre%>%mutate(variable="FTwre")%>%select(.chain,.iteration,.draw,p,value=b_FTwre,variable)
b_FTsre2<-b_FTsre%>%mutate(variable="FTsre")%>%select(.chain,.iteration,.draw,p,value=b_FTsre,variable)
b_TEMPcs2<-b_TEMPcs%>%mutate(variable="TEMPcs")%>%select(.chain,.iteration,.draw,p,value=b_TEMPcs,variable)
b_MLDcs2<-b_MLDcs%>%mutate(variable="MLDcs")%>%select(.chain,.iteration,.draw,p,value=b_MLDcs,variable)
b_TEMPsoo2<-b_TEMPsoo%>%mutate(variable="TEMPsoo")%>%select(.chain,.iteration,.draw,p,value=b_TEMPsoo,variable)
b_TEMPwoo2<-b_TEMPwoo%>%mutate(variable="TEMPwoo")%>%select(.chain,.iteration,.draw,p,value=b_TEMPwoo,variable)
b_FTum2<-b_FTum%>%mutate(variable="FTum")%>%select(.chain,.iteration,.draw,p,value=b_FTum,variable)
b_FDum2<-b_FDum%>%mutate(variable="FDum")%>%select(.chain,.iteration,.draw,p,value=b_FDum,variable)

pdf<-data.frame(rbind(b_FTwre2[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_FTsre2[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_TEMPcs2[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_MLDcs2[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_TEMPwoo2[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_TEMPsoo2[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_FTum2[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_FDum2[,c(".chain",".iteration",".draw","p","value","variable")]))


pdf$p <- cuid_datx$cu_name[pdf$p]
pdf <- pdf%>%mutate(d = case_when(p %in% c("Great Central Lake","Sproat Lake", "Osoyoos Lake", "Wenatchee Lake") ~ "South",
                                  p %in% c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake") ~ "Fraser", 
                                  p  %in% c("Babine Lake", "Tahltan Lake", "Tatsamenie Lake") ~ "North"))
pdf$d <- factor(pdf$d, levels=c("Hyper","South","Fraser","North"))
pdf$variable <- factor(pdf$variable, levels=c("FTwre","FTsre","TEMPcs","MLDcs","TEMPwoo","TEMPsoo","FTum","FDum"))
pdf <- pdf%>%mutate(v = case_when(variable == "FTwre" ~ "Winter Freshwater Temperature",
                                  variable == "FTsre" ~ "Summer Freshwater Temperature",
                                  variable == "TEMPcs" ~ "Coastal Temperature",
                                  variable == "MLDcs" ~ "Coastal Mixed Layer Depth",
                                  variable == "TEMPwoo" ~ "Winter Open Ocean Temperature",
                                  variable == "TEMPsoo" ~ "Summer Open Ocean Temperature",
                                  variable == "FTum" ~ "Return Migration Temperature",
                                  variable == "FDum" ~ "Return Migration Discharge"))
pdf$v = factor(pdf$v,levels=c("Winter Freshwater Temperature","Summer Freshwater Temperature","Coastal Temperature","Coastal Mixed Layer Depth","Winter Open Ocean Temperature","Summer Open Ocean Temperature","Return Migration Temperature","Return Migration Discharge"))

ggplot()+
  geom_density(data=pdf,aes(x = value, group=p, colour = d),linewidth=0.5)+
  geom_density(data=mdf,aes(x = value, colour = d),linewidth=1.2)+
  scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3','Hyper'='grey29'))+
  theme_minimal(base_size=12)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title=NULL,colour=NULL,y="Posterior density",x="Parameter effect")+
  theme(legend.position='bottom')+
  geom_vline(xintercept=0, linetype= 2)+
  scale_x_continuous(limits=c(-1.5,1.5),breaks=c(-1,0,1))+
  scale_y_continuous(breaks=NULL)+
  facet_wrap(v ~ ., scales = "free_y",ncol=2)

#Effect Sizes####

###Pull out all variables with dimensions
#post<-extract(fit)
#post <- spread_draws(fit, a_smolts_log[fw_age, p], a_returns_logit[idx,p], b_smolts[p], b_returns[p], b_FTwre[p,w], b_FTsre[p,s], b_TEMPcs[p], b_MLDcs[p], b_TEMPwoo[p,w], b_TEMPsoo[p,w], b_FTum[p], b_FDum[p], FTwre_complete[y,p], FTsre_complete[y,p], FTum_complete[y,p], FDum_complete[y,p], TEMPcs_complete[y,r], MLDcs_complete[y,r], TEMPwoo_complete[y],TEMPsoo_complete[y])
#spread_draws(fit, a_smolts_log[fw_age, p])|>filter(fw_age==1,p==1)
#median(spread_draws(fit, FTwre_complete[y,p])|>filter(p==1)|>pull(FTwre_complete),na.rm=TRUE)
#a_smolts_log<-spread_draws(fit, a_smolts_log[f,p]) #removed in v4.6
#a_fry<-spread_draws(fit,a_fry[by,p])
a_fry_log_year_mu<-spread_draws(fit,a_fry_log_year_mu[p])
holdover_logit_mu<-spread_draws(fit,holdover_logit_mu[p])
age2_survival_logit_mu<-spread_draws(fit,age2_survival_logit_mu[p])
#b_smolts<-spread_draws(fit, b_smolts[p]) #changed from b_smolts[p] in v4.6
b_fry<-spread_draws(fit, b_fry[p]) #changed from b_smolts[p] in v4.6
a_returns_logit<-spread_draws(fit, a_returns_logit[idx,p])
#b_returns<-spread_draws(fit, b_returns[p]) #removed density dependent b_returns for the marine life stage!!!

b_FTwre<-spread_draws(fit, b_FTwre[p,w])
b_FTsre<-spread_draws(fit, b_FTsre[p,s])
b_TEMPcs<-spread_draws(fit, b_TEMPcs[p])
b_MLDcs<-spread_draws(fit, b_MLDcs[p])
b_TEMPwoo<-spread_draws(fit, b_TEMPwoo[p,w])
b_TEMPsoo<-spread_draws(fit, b_TEMPsoo[p,s])
b_FTum<-spread_draws(fit, b_FTum[p])
b_FDum<-spread_draws(fit, b_FDum[p])

FTwre<-spread_draws(fit, FTwre_complete[y,p])
FTsre<-spread_draws(fit, FTsre_complete[y,p])
FTum<-spread_draws(fit, FTum_complete[y,p])
FDum<-spread_draws(fit, FDum_complete[y,p])
TEMPcs<-spread_draws(fit, TEMPstnd[y,r])
MLDcs<-spread_draws(fit, MLDstnd[y,r])
TEMPwoo<-spread_draws(fit, TEMPwoo_complete[y])
TEMPsoo<-spread_draws(fit, TEMPsoo_complete[y])

#MLDstnd_a<-spread_draws(fit, MLDstnd_TEMPstnd_latent_a[r])
#MLDstnd_b<-spread_draws(fit, MLDstnd_TEMPstnd_b[r])
FTum_a<-spread_draws(fit, FTum_FDum_latent_a[p])
FTum_b<-spread_draws(fit, FTum_FDum_b[p])

###INCLUDE MODELED FISHING PROPORTION IN THE CALCULATION?

#esc_abd <- exp(model_data$e_priors_ln)
spawner_abd<-all_totals[!is.na(all_totals$spawners),c('cu_id','spawners')]%>% group_by(cu_id) %>%
  summarize(smean = mean(spawners, na.rm = TRUE)) %>%
  pivot_wider(names_from = "cu_id", values_from = "smean") %>%
  ungroup() %>%
  as.matrix()%>%
  as.vector()

#add domain to env_dat
# env_dat<-env_dat%>%filter(cu %in% cus)%>%mutate(domain = case_when(cu %in% c("Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake") ~ "South",
#                                              cu %in% c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake") ~ "Fraser", 
#                                              cu %in% c("Babine Lake","Tahltan Lake","Tatsamenie Lake") ~ "North")) |> mutate(domain = factor(domain, levels = c("South", "Fraser", "North")))

#Number of simulated env data points 
n<-50

#TEMPcs####
TEMP_unscaled <- coastal_dat%>%filter(!region=="OpenOcean")%>%select(region,year,TEMP)  
TEMPcs_complete <- spread_draws(fit, TEMPcs_complete[y,r])%>%group_by(y,r)%>%summarize(value=mean(TEMPcs_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),major_watershed=factor(unique(cuid_dat$major_watershed)[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(major_watershed, year, value)  
TEMPstnd <- spread_draws(fit, TEMPstnd[y,r])%>%group_by(y,r)%>%summarize(value=mean(TEMPstnd,na.rm=TRUE))%>%mutate(year=y+(year_start-1),major_watershed=factor(unique(cuid_dat$major_watershed)[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(major_watershed, year, value)  

TEMPcs_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  #TEMPcs_unscaled <- env_dat%>%filter(domain==cuid_datx$domain[pop])%>%select(TEMPcsMean)
  #TEMPcs_unscaled <- env_select%>%select(TEMPcs)
  TEMPcs_completex <- as.data.frame(TEMPcs_complete%>%filter(major_watershed==cuid_datx$major_watershed[pop])%>%select(value))
  #maxx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='TEMPcs')%>%summarise(maxv = max(value_stnd, na.rm = TRUE))%>%.[[1]]
  #minx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='TEMPcs')%>%summarise(minv = min(value_stnd, na.rm = TRUE))%>%.[[1]]
  maxx=TEMPstnd%>%filter(major_watershed==cuid_datx$major_watershed[pop])%>%summarise(maxv = max(value, na.rm = TRUE))%>%.[[1]]
  minx=TEMPstnd%>%filter(major_watershed==cuid_datx$major_watershed[pop])%>%summarise(minv = min(value, na.rm = TRUE))%>%.[[1]]
  TEMPstnd_x <- seq(minx, maxx, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    #MLD_x <- MLDstnd_a|>filter(r==cuid_datx$major_watershed_id[pop])|>pull(MLDstnd_TEMPstnd_latent_a) + MLDstnd_b|>filter(r==cuid_datx$major_watershed_id[pop])|>pull(MLDstnd_TEMPstnd_b) * TEMPstnd_x[x]
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * 
        FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * 
                               FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * 
                               FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre))
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1 = smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smolts_2 <- smolts}
      for(mar in 1:model_data$N_mar_ages){
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * TEMPstnd_x[x]+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * median(MLDcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(MLDstnd),na.rm=TRUE)+
                          #b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * MLD_x+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * median(FTum|>filter(p==pop)|>pull(FTum_complete),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    TEMPcs_contrast.df <- bind_rows(TEMPcs_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], TEMPcs = (TEMPstnd_x[x] * sd(TEMPcs_completex$value, na.rm = TRUE) + mean(TEMPcs_completex$value, na.rm = TRUE))*sd(TEMP_unscaled$TEMP,na.rm=TRUE)+mean(TEMP_unscaled$TEMP,na.rm=TRUE), x = TEMPstnd_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
TEMPcs_contrast.df %>% 
  ggplot(aes(x = TEMPcs, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  facet_wrap(~domain)+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Coastal Temperature (Celsius)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2) 

#sample plot
ggplot(data=TEMPcs_contrast.df%>%filter(pop%in%c('Osoyoos Lake','Fraser Lake','Tatsamenie Lake')), aes(x = TEMPcs, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=c('Osoyoos Lake'='darkorchid3','Fraser Lake'='goldenrod3','Tatsamenie Lake'='seagreen3'))+
  #facet_wrap(~domain)+
  theme_bw(base_size=14)+
  theme(legend.position = 'bottom')+
  labs(title=NULL,colour=NULL,x="Coastal Temperature (Celsius)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change 
#stndx<-env_long %>% filter(variable == "TEMPcs") %>% group_by(cu_name) %>% summarise(stnd_mean = mean(value_stnd,na.rm=TRUE), stnd_sd = sd(value_stnd,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
stndx<-TEMPstnd %>% group_by(major_watershed) %>% summarise(stnd_mean = mean(value,na.rm=TRUE), stnd_sd = sd(value,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
TEMPcs_contrast.df2<-data.frame()

for(popx in 1:model_data$N_pop){
  TEMPcs_contrast.x<-TEMPcs_contrast.df%>%filter(pop==cuid_datx$cu_name[popx])%>%select(pop,replace,draw,x,n)
  uni = unique(TEMPcs_contrast.x$x)
  #minx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))))
  #maxx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))))
  minx<-which(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdm15)%>%.[[1]]))))
  maxx<-which(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdp15)%>%.[[1]]))))
  TEMPcs_contrast.x<-TEMPcs_contrast.x%>%filter(n %in% c(minx,maxx))%>%
    select(pop,replace,draw,n)%>%
    pivot_wider(names_from=n,values_from=replace)%>%
    mutate(per=.[[4]] / .[[3]] * 100)%>%
    select(pop, draw,per)
  TEMPcs_contrast.df2<-rbind(TEMPcs_contrast.df2,TEMPcs_contrast.x)
}
TEMPcs_contrast.df2<-TEMPcs_contrast.df2%>%left_join(TEMPcs_contrast.df%>%filter(n==1))
TEMPcs_contrast.df2$domain<-factor(cuid_datx$domain[TEMPcs_contrast.df2$pop],levels=c("South","Fraser","North"))
TEMPcs_contrast.df2$popx=factor(str_remove(as.character(TEMPcs_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)

#MLDcs####
#Optionally with the indirect of MLD>Temp
MLD_unscaled <- coastal_dat%>%filter(!region=="OpenOcean")%>%select(region,year,MLD)  
MLDcs_complete <- spread_draws(fit, MLDcs_complete[y,r])%>%group_by(y,r)%>%summarize(value=mean(MLDcs_complete,na.rm=TRUE))%>%mutate(year=y+(year_start-1),major_watershed=factor(unique(cuid_dat$major_watershed)[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(major_watershed, year, value)  
MLDstnd <- spread_draws(fit, MLDstnd[y,r])%>%group_by(y,r)%>%summarize(value=mean(MLDstnd,na.rm=TRUE))%>%mutate(year=y+(year_start-1),major_watershed=factor(unique(cuid_dat$major_watershed)[r],levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku")))%>%ungroup()%>%select(major_watershed, year, value)  

MLDcs_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  #MLDcs_unscaled <- env_dat%>%filter(domain==cuid_datx$domain[pop])%>%select(MLDcsMean2)
  #MLDcs_unscaled <- env_select%>%filter(major_watershed==cuid_datx$major_watershed[pop])%>%select(MLDcs)
  #maxx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='MLDcs')%>%summarise(maxv = max(value_stnd, na.rm = TRUE))%>%.[[1]]
  #minx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='MLDcs')%>%summarise(minv = min(value_stnd, na.rm = TRUE))%>%.[[1]]
  MLDcs_completex <- as.data.frame(MLDcs_complete%>%filter(major_watershed==cuid_datx$major_watershed[pop])%>%select(value))
  maxx=MLDstnd%>%filter(major_watershed==cuid_datx$major_watershed[pop])%>%summarise(maxv = max(value, na.rm = TRUE))%>%.[[1]]
  minx=MLDstnd%>%filter(major_watershed==cuid_datx$major_watershed[pop])%>%summarise(minv = min(value, na.rm = TRUE))%>%.[[1]]
  MLDstnd_x <- seq(minx, maxx, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    #TEMP_x <- MLDstnd_a|>filter(r==cuid_datx$major_watershed_id[pop])|>pull(MLDstnd_TEMPstnd_latent_a) + MLDstnd_b|>filter(r==cuid_datx$major_watershed_id[pop])|>pull(MLDstnd_TEMPstnd_b) * MLDstnd_x[x]
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * 
        FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * 
                               FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * 
                               FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre))
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1=smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smolts_2=smolts}
      for(mar in 1:model_data$N_mar_ages){
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * MLDstnd_x[x]+
                          #b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * TEMP_x+
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * median(TEMPcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(TEMPstnd),na.rm=TRUE)+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * median(FTum|>filter(p==pop)|>pull(FTum_complete),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    MLDcs_contrast.df <- bind_rows(MLDcs_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], MLDcs = (MLDstnd_x[x] * sd(MLDcs_completex$value, na.rm = TRUE) + mean(MLDcs_completex$value, na.rm = TRUE))*sd(MLD_unscaled$MLD,na.rm=TRUE)+mean(MLD_unscaled$MLD,na.rm=TRUE), x = MLDstnd_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
MLDcs_contrast.df %>% 
  ggplot(aes(x = MLDcs, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  facet_wrap(~domain,scale="free_x")+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Coastal Mixed Layer Depth (Meters)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change 
#stndx<-env_long %>% filter(variable == "MLDcs") %>% group_by(cu_name) %>% summarise(stnd_mean = mean(value_stnd,na.rm=TRUE), stnd_sd = sd(value_stnd,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
stndx<-MLDstnd %>% group_by(major_watershed) %>% summarise(stnd_mean = mean(value,na.rm=TRUE), stnd_sd = sd(value,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
MLDcs_contrast.df2<-data.frame()

for(popx in 1:model_data$N_pop){
  MLDcs_contrast.x<-MLDcs_contrast.df%>%filter(pop==cuid_datx$cu_name[popx])%>%select(pop,replace,draw,x,n)
  uni = unique(MLDcs_contrast.x$x)
  #minx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))))
  #maxx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))))
  minx<-which(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdm15)%>%.[[1]]))))
  maxx<-which(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(major_watershed==cuid_datx$major_watershed[popx])%>%select(sdp15)%>%.[[1]]))))
  MLDcs_contrast.x<-MLDcs_contrast.x%>%filter(n %in% c(minx,maxx))%>%
    select(pop,replace,draw,n)%>%
    pivot_wider(names_from=n,values_from=replace)%>%
    mutate(per=.[[4]] / .[[3]] * 100)%>%
    select(pop, draw,per)
  MLDcs_contrast.df2<-rbind(MLDcs_contrast.df2,MLDcs_contrast.x)
}
MLDcs_contrast.df2<-MLDcs_contrast.df2%>%left_join(MLDcs_contrast.df%>%filter(n==1))
MLDcs_contrast.df2$domain<-factor(cuid_datx$domain[MLDcs_contrast.df2$pop],levels=c("South","Fraser","North"))
MLDcs_contrast.df2$popx=factor(str_remove(as.character(MLDcs_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)

#TEMPwoo####
TEMPwoo_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  TEMPwoo_unscaled <- env_select2$TEMPwoo
  maxx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='TEMPwoo')%>%summarise(maxv = max(value_stnd, na.rm = TRUE))%>%.[[1]]
  minx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='TEMPwoo')%>%summarise(minv = min(value_stnd, na.rm = TRUE))%>%.[[1]]
  TEMPwoo_x <- seq(minx, maxx, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * 
        FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * 
                               FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * 
                               FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre))
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1=smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smolts_2=smolts}
      for(mar in 1:model_data$N_mar_ages){
        if(mar==1){b_TEMPwoox = b_TEMPwoo|>filter(p==pop,w==1)|>pull(b_TEMPwoo) * TEMPwoo_x[x]}
        if(mar==2){b_TEMPwoox = b_TEMPwoo|>filter(p==pop,w==1)|>pull(b_TEMPwoo) * TEMPwoo_x[x] + b_TEMPwoo|>filter(p==pop,w==2)|>pull(b_TEMPwoo) * TEMPwoo_x[x]}
        if(mar==3){b_TEMPwoox = b_TEMPwoo|>filter(p==pop,w==1)|>pull(b_TEMPwoo) * TEMPwoo_x[x] + b_TEMPwoo|>filter(p==pop,w==2)|>pull(b_TEMPwoo) * TEMPwoo_x[x] + b_TEMPwoo|>filter(p==pop,w==3)|>pull(b_TEMPwoo) * TEMPwoo_x[x]}
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_TEMPwoox +
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * median(TEMPcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(TEMPstnd),na.rm=TRUE)+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * median(MLDcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(MLDstnd),na.rm=TRUE)+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * median(FTum|>filter(p==pop)|>pull(FTum_complete),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    TEMPwoo_contrast.df <- bind_rows(TEMPwoo_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], TEMPwoo = TEMPwoo_x[x] * sd(TEMPwoo_unscaled, na.rm = TRUE) + mean(TEMPwoo_unscaled, na.rm = TRUE), x = TEMPwoo_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
TEMPwoo_contrast.df %>% 
  ggplot(aes(x = TEMPwoo, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  facet_wrap(~domain,scale="free_y")+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Winter Open Ocean Temperature (Celsius)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change 
stndx<-env_long %>% filter(variable == "TEMPwoo") %>% group_by(cu_name) %>% summarise(stnd_mean = mean(value_stnd,na.rm=TRUE), stnd_sd = sd(value_stnd,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
TEMPwoo_contrast.df2<-data.frame()

for(popx in 1:model_data$N_pop){
  TEMPwoo_contrast.x<-TEMPwoo_contrast.df%>%filter(pop==cuid_datx$cu_name[popx])%>%select(pop,replace,draw,x,n)
  uni = unique(TEMPwoo_contrast.x$x)
  minx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))))
  maxx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))))
  TEMPwoo_contrast.x<-TEMPwoo_contrast.x%>%filter(n %in% c(minx,maxx))%>%
    select(pop,replace,draw,n)%>%
    pivot_wider(names_from=n,values_from=replace)%>%
    mutate(per=.[[4]] / .[[3]] * 100)%>%
    select(pop, draw,per)
  TEMPwoo_contrast.df2<-rbind(TEMPwoo_contrast.df2,TEMPwoo_contrast.x)
}
TEMPwoo_contrast.df2<-TEMPwoo_contrast.df2%>%left_join(TEMPwoo_contrast.df%>%filter(n==1))
TEMPwoo_contrast.df2$domain<-factor(cuid_datx$domain[TEMPwoo_contrast.df2$pop],levels=c("South","Fraser","North"))
TEMPwoo_contrast.df2$popx=factor(str_remove(as.character(TEMPwoo_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)

#TEMPsoo####
TEMPsoo_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  TEMPsoo_unscaled <- env_select2$TEMPsoo
  maxx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='TEMPsoo')%>%summarise(maxv = max(value_stnd, na.rm = TRUE))%>%.[[1]]
  minx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='TEMPsoo')%>%summarise(minv = min(value_stnd, na.rm = TRUE))%>%.[[1]]
  TEMPsoo_x <- seq(minx, maxx, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * 
        FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * 
                               FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * 
                               FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre))
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1=smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smolts_2=smolts}
      for(mar in 1:model_data$N_mar_ages){
        if(mar==1){b_TEMPsoox = b_TEMPsoo|>filter(p==pop,s==1)|>pull(b_TEMPsoo) * TEMPsoo_x[x]}
        if(mar==2){b_TEMPsoox = b_TEMPsoo|>filter(p==pop,s==1)|>pull(b_TEMPsoo) * TEMPsoo_x[x] + b_TEMPsoo|>filter(p==pop,s==2)|>pull(b_TEMPsoo) * TEMPsoo_x[x]}
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_TEMPsoox +
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * median(TEMPcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(TEMPstnd),na.rm=TRUE)+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * median(MLDcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(MLDstnd),na.rm=TRUE)+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * median(FTum|>filter(p==pop)|>pull(FTum_complete),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    TEMPsoo_contrast.df <- bind_rows(TEMPsoo_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], TEMPsoo = TEMPsoo_x[x] * sd(TEMPsoo_unscaled, na.rm = TRUE) + mean(TEMPsoo_unscaled, na.rm = TRUE), x = TEMPsoo_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
TEMPsoo_contrast.df %>% 
  ggplot(aes(x = TEMPsoo, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  facet_wrap(~domain,scale="free_y")+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Summer Open Ocean Temperature (Celsius)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change 
stndx<-env_long %>% filter(variable == "TEMPsoo") %>% group_by(cu_name) %>% summarise(stnd_mean = mean(value_stnd,na.rm=TRUE), stnd_sd = sd(value_stnd,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
TEMPsoo_contrast.df2<-data.frame()

for(popx in 1:model_data$N_pop){
  TEMPsoo_contrast.x<-TEMPsoo_contrast.df%>%filter(pop==cuid_datx$cu_name[popx])%>%select(pop,replace,draw,x,n)
  uni = unique(TEMPsoo_contrast.x$x)
  minx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))))
  maxx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))))
  TEMPsoo_contrast.x<-TEMPsoo_contrast.x%>%filter(n %in% c(minx,maxx))%>%
    select(pop,replace,draw,n)%>%
    pivot_wider(names_from=n,values_from=replace)%>%
    mutate(per=.[[4]] / .[[3]] * 100)%>%
    select(pop, draw,per)
  TEMPsoo_contrast.df2<-rbind(TEMPsoo_contrast.df2,TEMPsoo_contrast.x)
}
TEMPsoo_contrast.df2<-TEMPsoo_contrast.df2%>%left_join(TEMPsoo_contrast.df%>%filter(n==1))
TEMPsoo_contrast.df2$domain<-factor(cuid_datx$domain[TEMPsoo_contrast.df2$pop],levels=c("South","Fraser","North"))
TEMPsoo_contrast.df2$popx=factor(str_remove(as.character(TEMPsoo_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)

#FTum####
FTum_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  FTum_unscaled <- env_select2%>%filter(cu_name==cuid_datx$cu_name[pop])%>%select(FTum)
  maxx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='FTum')%>%summarise(maxv = max(value_stnd, na.rm = TRUE))%>%.[[1]]
  minx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='FTum')%>%summarise(minv = min(value_stnd, na.rm = TRUE))%>%.[[1]]
  FTum_x <- seq(minx, maxx, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * 
        FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * 
                               FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * 
                               FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre))
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1=smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smolts_2=smolts}
      for(mar in 1:model_data$N_mar_ages){
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * FTum_x[x]+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * median(MLDcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(MLDstnd),na.rm=TRUE)+
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * median(TEMPcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(TEMPstnd),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    FTum_contrast.df <- bind_rows(FTum_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], FTum = FTum_x[x] * sd(FTum_unscaled[,1], na.rm = TRUE) + mean(FTum_unscaled[,1], na.rm = TRUE), x = FTum_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
FTum_contrast.df %>% 
  ggplot(aes(x = FTum, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  facet_wrap(~domain)+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Return Migration Temperature (Celsius)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change 
stndx<-env_long %>% filter(variable == "FTum") %>% group_by(cu_name) %>% summarise(stnd_mean = mean(value_stnd,na.rm=TRUE), stnd_sd = sd(value_stnd,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
FTum_contrast.df2<-data.frame()

for(popx in 1:model_data$N_pop){
  FTum_contrast.x<-FTum_contrast.df%>%filter(pop==cuid_datx$cu_name[popx])%>%select(pop,replace,draw,x,n)
  uni = unique(FTum_contrast.x$x)
  minx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))))
  maxx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))))
  FTum_contrast.x<-FTum_contrast.x%>%filter(n %in% c(minx,maxx))%>%
    select(pop,replace,draw,n)%>%
    pivot_wider(names_from=n,values_from=replace)%>%
    mutate(per=.[[4]] / .[[3]] * 100)%>%
    select(pop, draw,per)
  FTum_contrast.df2<-rbind(FTum_contrast.df2,FTum_contrast.x)
}
FTum_contrast.df2<-FTum_contrast.df2%>%left_join(FTum_contrast.df%>%filter(n==1))
FTum_contrast.df2$domain<-factor(cuid_datx$domain[FTum_contrast.df2$pop],levels=c("South","Fraser","North"))
FTum_contrast.df2$popx=factor(str_remove(as.character(FTum_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)

#FDum####
#maxx=max(model_data$FDum)
#FDum_x <- seq(-maxx, maxx, length = n)
FDum_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  FDum_unscaled <- env_select2%>%filter(cu_name==cuid_datx$cu_name[pop])%>%select(FDum)
  FDum_x <- seq(-1.5, 1.5, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    FTum_xd = FTum_a|>filter(p==cuid_datx$major_watershed_id[pop])|>pull(FTum_FDum_latent_a)+FTum_b|>filter(p==cuid_datx$major_watershed_id[pop])|>pull(FTum_FDum_b)*FDum_x
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * 
        FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * 
                               FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * 
                               FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre))
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1=smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smolts_2=smolts}
      for(mar in 1:model_data$N_mar_ages){
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_FDum|>filter(p==pop)|>pull(b_FDum) * FDum_x[x]+
                          #b_FTum|>filter(p==pop)|>pull(b_FTum) * median(FTum|>filter(p==pop)|>pull(FTum_complete),na.rm=TRUE)+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * FTum_xd+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * median(MLDcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(MLDstnd),na.rm=TRUE)+
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * median(TEMPcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(TEMPstnd),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    FDum_contrast.df <- bind_rows(FDum_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], FDum = FDum_x[x] * sd(FDum_unscaled[,1], na.rm = TRUE) + mean(FDum_unscaled[,1], na.rm = TRUE), x = FDum_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
FDum_contrast.df %>% 
  ggplot(aes(x = FDum, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  facet_wrap(~domain, scale='free_x')+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Return Migration Discharge (m^3 s^-1)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change, different for FDum because it is normalized by cu and not globally 
minx<-which(abs(FDum_x-(-1.5))==min(abs(FDum_x-(-1.5)))) #-1.5sd
maxx<-which(abs(FDum_x-(1.5))==min(abs(FDum_x-(1.5)))) #+1.5sd
FDum_contrast.df2 <- FDum_contrast.df%>%select(pop,replace,draw,n)%>%pivot_wider(names_from=n,values_from=replace)%>%select(1,2,minx+2,maxx+2)
#FDum_x <- FDum_contrast.df%>%select(pop,replace,draw,n)%>%pivot_wider(names_from=n,values_from=replace)
#FDum_perx <- FDum_perx%>%mutate(across(3:52,~.x/FDum_perx%>%pull(3)*100))
#FDum_contrast.df2 <- FDum_perx%>%pivot_longer(cols=c(-pop,-draw),names_to='n',values_to='per')%>%mutate(n=as.numeric(n))%>%left_join(FDum_contrast.df)
FDum_contrast.df2 <- FDum_contrast.df2%>%mutate(per=.[[4]] / .[[3]] * 100)
FDum_contrast.df2<-FDum_contrast.df2%>%left_join(FDum_contrast.df%>%filter(n==1))
FDum_contrast.df2$domain<-factor(cuid_datx$domain[FDum_contrast.df2$pop],levels=c("South","Fraser","North"))
FDum_contrast.df2$popx=factor(str_remove(as.character(FDum_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)

#FTwre####
FTwre_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  FTwre_unscaled <- env_select2$FTwre
  maxx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='FTwre')%>%summarise(maxv = max(value_stnd, na.rm = TRUE))%>%.[[1]]
  minx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='FTwre')%>%summarise(minv = min(value_stnd, na.rm = TRUE))%>%.[[1]]
  FTwre_x <- seq(minx, maxx, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * FTwre_x[x]+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * 
        FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * FTwre_x[x]
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * 
                               FTsre|>filter(p==pop) |> group_by(.draw) |> summarise(FTsre = median(FTsre_complete))|>pull(FTsre) +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * FTwre_x[x])
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1=smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smotls_2=smolts}
      for(mar in 1:model_data$N_mar_ages){
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * median(FTum|>filter(p==pop)|>pull(FTum_complete),na.rm=TRUE)+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * median(MLDcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(MLDstnd),na.rm=TRUE)+
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * median(TEMPcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(TEMPstnd),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    FTwre_contrast.df <- bind_rows(FTwre_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], FTwre = FTwre_x[x] * sd(FTwre_unscaled, na.rm = TRUE) + mean(FTwre_unscaled, na.rm = TRUE), x = FTwre_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
FTwre_contrast.df %>% 
  ggplot(aes(x = FTwre, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  scale_y_log10()+
  facet_wrap(~domain,scale='free')+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Winter Freshwater Temperature (Celsius)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change 
stndx<-env_long %>% filter(variable == "FTwre") %>% group_by(cu_name) %>% summarise(stnd_mean = mean(value_stnd,na.rm=TRUE), stnd_sd = sd(value_stnd,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
FTwre_contrast.df2<-data.frame()

for(popx in 1:model_data$N_pop){
  FTwre_contrast.x<-FTwre_contrast.df%>%filter(pop==cuid_datx$cu_name[popx])%>%select(pop,replace,draw,x,n)
  uni = unique(FTwre_contrast.x$x)
  minx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))))
  maxx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))))
  FTwre_contrast.x<-FTwre_contrast.x%>%filter(n %in% c(minx,maxx))%>%
    select(pop,replace,draw,n)%>%
    pivot_wider(names_from=n,values_from=replace)%>%
    mutate(per=.[[4]] / .[[3]] * 100)%>%
    select(pop, draw,per)
  FTwre_contrast.df2<-rbind(FTwre_contrast.df2,FTwre_contrast.x)
}
FTwre_contrast.df2<-FTwre_contrast.df2%>%left_join(FTwre_contrast.df%>%filter(n==1))
FTwre_contrast.df2$domain<-factor(cuid_datx$domain[FTwre_contrast.df2$pop],levels=c("South","Fraser","North"))
FTwre_contrast.df2$popx=factor(str_remove(as.character(FTwre_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)

#FTsre####
FTsre_contrast.df <- data.frame()
for(pop in 1:model_data$N_pop){
  print(pop)
  FTsre_unscaled <- env_select2$FTsre
  maxx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='FTsre')%>%summarise(maxv = max(value_stnd, na.rm = TRUE))%>%.[[1]]
  minx=env_long%>%filter(cu_name==cuid_datx$cu_name[pop])%>%filter(variable=='FTsre')%>%summarise(minv = min(value_stnd, na.rm = TRUE))%>%.[[1]]
  FTsre_x <- seq(minx, maxx, length = n)
  for(x in 1:n){
    returns <- 0
    smolts <- 0
    for(fw in 1:model_data$N_fw_ages){
      log_fw = a_fry_log_year_mu |> filter(p == pop) |> pull(a_fry_log_year_mu) +
        b_FTwre|>filter(p==pop,w==1)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)+
        b_FTsre|>filter(p==pop,s==1)|>pull(b_FTsre) * FTsre_x[x] +
        b_FTwre|>filter(p==pop,w==2)|>pull(b_FTwre) * 
        FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre)
      
      a_fry = exp(log_fw)
      
      holdover = plogis(holdover_logit_mu |> filter(p == pop) |> pull(holdover_logit_mu))
      
      age2_survival = plogis(age2_survival_logit_mu |> filter(p == pop) |> pull(age2_survival_logit_mu)+
                               b_FTsre|>filter(p==pop,s==2)|>pull(b_FTsre) * FTsre_x[x] +
                               b_FTwre|>filter(p==pop,w==3)|>pull(b_FTwre) * 
                               FTwre|>filter(p==pop) |> group_by(.draw) |> summarise(FTwre = median(FTwre_complete))|>pull(FTwre))
      
      fry = (a_fry *  spawner_abd[pop]) / (1 + b_fry |> filter(p == pop) |> pull(b_fry) * ( spawner_abd[pop] / 1e6));
      if(fw==1){smolts = fry * (1 - holdover);smolts_1=smolts}
      if(fw==2){smolts = fry * holdover * age2_survival;smolts_2=smolts}
      for(mar in 1:model_data$N_mar_ages){
        alpha <- plogis(a_returns_logit|>filter(idx==((fw-1)*model_data$N_mar_ages+mar),p==pop)|>pull(a_returns_logit)+
                          b_FTum|>filter(p==pop)|>pull(b_FTum) * median(FTum|>filter(p==pop)|>pull(FTum_complete),na.rm=TRUE)+
                          b_MLDcs|>filter(p==pop)|>pull(b_MLDcs) * median(MLDcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(MLDstnd),na.rm=TRUE)+
                          b_TEMPcs|>filter(p==pop)|>pull(b_TEMPcs) * median(TEMPcs|>filter(r==cuid_datx[which(cuid_datx$cu_id==pop),"major_watershed_id"])|>pull(TEMPstnd),na.rm=TRUE))
        #returns <- returns + (alpha * smolts)/(1 + (b_returns|>filter(p==pop)|>pull(b_returns) * smolts/1e6))
        returns <- returns + (alpha * smolts)
      }
    }
    FTsre_contrast.df <- bind_rows(FTsre_contrast.df, data.frame(pop = factor(cuid_datx$cu_name[pop],levels=levelx), domain = factor(cuid_datx$domain[pop],levels=c("South","Fraser","North")), smolts = smolts_1+smolts_2, returns = returns, replace=returns/spawner_abd[pop], FTsre = FTsre_x[x] * sd(FTsre_unscaled, na.rm = TRUE) + mean(FTsre_unscaled, na.rm = TRUE), x = FTsre_x[x], n=x, draw  = 1:length(returns)))
  }
}
#plot
FTsre_contrast.df %>% 
  ggplot(aes(x = FTsre, y = replace, group = pop, color = pop))+
  stat_lineribbon(.width = c(0.5))+
  scale_fill_manual(values = "grey90", guide = "none")+
  scale_color_manual(values=mcolx)+
  scale_y_log10()+
  facet_wrap(~domain,scale='free')+
  theme_bw(base_size=12)+
  labs(title=NULL,colour="Populations",x="Summer Freshwater Temperature (Celsius)",y="Productivity")+
  geom_hline(yintercept = 1, lty = 2)

#percentage change 
stndx<-env_long %>% filter(variable == "FTsre") %>% group_by(cu_name) %>% summarise(stnd_mean = mean(value_stnd,na.rm=TRUE), stnd_sd = sd(value_stnd,na.rm=TRUE))%>%mutate(sdm1=stnd_mean-stnd_sd, sdp1=stnd_mean+stnd_sd, sdm15=stnd_mean-(1.5*stnd_sd), sdp15=stnd_mean+(1.5*stnd_sd))
FTsre_contrast.df2<-data.frame()

for(popx in 1:model_data$N_pop){
  FTsre_contrast.x<-FTsre_contrast.df%>%filter(pop==cuid_datx$cu_name[popx])%>%select(pop,replace,draw,x,n)
  uni = unique(FTsre_contrast.x$x)
  minx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdm15)%>%.[[1]]))))
  maxx<-which(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))==min(abs(uni-(stndx%>%filter(cu_name==cuid_datx$cu_name[popx])%>%select(sdp15)%>%.[[1]]))))
  FTsre_contrast.x<-FTsre_contrast.x%>%filter(n %in% c(minx,maxx))%>%
    select(pop,replace,draw,n)%>%
    pivot_wider(names_from=n,values_from=replace)%>%
    mutate(per=.[[4]] / .[[3]] * 100)%>%
    select(pop, draw,per)
  FTsre_contrast.df2<-rbind(FTsre_contrast.df2,FTsre_contrast.x)
}
FTsre_contrast.df2<-FTsre_contrast.df2%>%left_join(FTsre_contrast.df%>%filter(n==1))
FTsre_contrast.df2$domain<-factor(cuid_datx$domain[FTsre_contrast.df2$pop],levels=c("South","Fraser","North"))
FTsre_contrast.df2$popx=factor(str_remove(as.character(FTsre_contrast.df2$pop)," Lake"),levels=cuid_datx$pop_name_short)


###Combined Effects
library(patchwork)

axis_range <- c(2,3000)
minorb <- c(seq(2,10,1),seq(20,90,10),seq(100,900,100),seq(1000,3000,1000))
majorb <- c(5,10,50,100,500,1000)

ggplot(FTwre_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title="Winter Freshwater Temperature",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0)) +
  #scale_y_log10(limits=c(20,500), minor_breaks=minorb, breaks=c(20,30,40,50,100,500), expand=expansion(mult=c(0.05,0.05)))

  ggplot(FTsre_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title="Summer Freshwater Temperature",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0)) +
  #scale_y_log10(limits=c(30,300), minor_breaks=minorb, breaks=c(30,40,50,100,200), expand=expansion(mult=c(0.05,0.05)))
  
  ggplot(TEMPcs_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title="Coastal Temperature",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0)) +
  #scale_y_log10(limits=c(5,1000), minor_breaks=minorb, breaks=c(10,20,30,40,50,100,200,300,500), expand=expansion(mult=c(0.05,0.05)))
  
  ggplot(MLDcs_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title="Coastal Mixed Layer Depth",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0)) +
  #scale_y_log10(limits=c(6,600), minor_breaks=minorb, breaks=c(10,20,30,40,50,100,200,300), expand=expansion(mult=c(0.05,0.05)))
  
  ggplot(TEMPwoo_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title= "Winter Open Ocean Temperature",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0)) +
  #scale_y_log10(limits=c(1,1000), minor_breaks=minorb, breaks=c(10,20,30,40,50,100,200,500), expand=expansion(mult=c(0.05,0.05)))
  
  ggplot(TEMPsoo_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title="Summer Open Ocean Temperature",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0)) +
  #scale_y_log10(limits=c(1,2000), minor_breaks=minorb, breaks=c(10,20,30,40,50,100,200,500,1000), expand=expansion(mult=c(0.05,0.05)))
  
  ggplot(FTum_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title="Return Migration Temperature",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0)) +
  #scale_y_log10(limits=c(40,200), minor_breaks=minorb, breaks=c(40,50,100,200), expand=expansion(mult=c(0.05,0.05)))
  
  ggplot(FDum_contrast.df2, aes(x = popx, y = per, fill = domain))+
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  theme_bw(base_size=12)+
  #theme_bw(base_size=14)+
  labs(title="Return Migration Discharge",fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0))+ 
  #scale_y_log10(limits=c(5,500), minor_breaks=minorb, breaks=c(10,20,30,40,50,100,200,300,500), expand=expansion(mult=c(0.05,0.05)))
  
  plot_layout(ncol=2, axis_titles='collect_x',axes='collect_y',guides='collect')&theme(legend.position='bottom')


#Capacities for freshwater model####
a_fry_log_year_mu <- spread_draws(fit,a_fry_log_year_mu[p])
a_fry <- spread_draws(fit,a_fry[y,p])
a_fry <- a_fry|>select(-y)|>group_by(p,.chain,.iteration,.draw)|>summarize(a_fry = mean(a_fry))|>ungroup()
b_fry <- spread_draws(fit, b_fry[p]) #changed from b_smolts[p] in v4.6

#fry mu capacity
c_fry_mu <- a_fry_log_year_mu|>left_join(b_fry)
c_fry_mu <- c_fry_mu|>mutate(smolt_capacity = (exp(a_fry_log_year_mu) - 1) / (b_fry / 1e6))
c_fry_mu$cu_name <- cuid_datx$cu_name[c_fry_mu$p]

ggplot()+geom_violin(data=c_fry_mu, aes(x=cu_name, y=smolt_capacity),fill='dodgerblue1', position="dodge", alpha=0.5)+
  scale_y_log10()+
  labs(x=NULL, y="Capacity")+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))

c_fry_x <- c_fry_mu|>select(cu_name,smolt_capacity)|>group_by(cu_name)|>summarize(smolt_capacity = mean(smolt_capacity,na.rm=TRUE))

ggplot()+
  stat_pointinterval(data=spredictions,aes(x = year, y = smolts_total,colour=domain),shape=19,.width=0.66,linewidth=1,size=7)+
  #stat_pointinterval(data=spredictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = smolts_total,colour=domain),shape=19,.width=0.66,linewidth=2,size=9)+
  #stat_pointinterval(data=spredictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = smolts_total,colour=domain),shape=NA,.width=0.66,linewidth=NA,size=9)+
  scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3','observed'='grey29'))+
  geom_point(data=juv_data,aes(x=smolt_migration_year,y=juv_abund,colour=domain),shape=18,size=3)+
  #geom_point(data=juv_data%>%filter(cu_name=="Sproat Lake"),aes(x=smolt_migration_year,y=juv_abund,colour=domain),shape=18,size=5)+
  geom_hline(data=c_fry_x,aes(yintercept=smolt_capacity),linetype='dashed',lwd=0.8,colour='red3')+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme_bw(base_size=10)+
  #theme_bw(base_size=14)+
  labs(title=NULL,x=NULL,y="Smolts",colour=NULL)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme(legend.position='bottom')

#returns capacity
#a_returns <- spread_draws(fit, a_returns[y, fw_age, mar_age, p])|>select(-y,-fw_age,-mar_age)|>group_by(p,.chain,.iteration,.draw)|>summarize(a_returns = mean(a_returns))|>ungroup()
#c_returns <- c_fry|>left_join(a_returns)|>mutate(returns_capacity = smolt_capacity * a_returns)
a_returns_logit <- spread_draws(fit, a_returns_logit[a, p])
c_returns <- a_returns_logit|>left_join(c_fry_mu)|>mutate(returns_capacity = smolt_capacity * plogis(a_returns_logit))|>group_by(cu_name,.chain,.iteration,.draw)|>summarize(returns_capacity = sum(returns_capacity))

ggplot()+geom_violin(data=c_returns, aes(x=cu_name, y=returns_capacity),fill='dodgerblue1', position="dodge", alpha=0.5)+
  scale_y_log10()+
  labs(x=NULL, y="Capacity")+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))

c_returns_x <- c_returns|>select(cu_name,returns_capacity)|>group_by(cu_name)|>summarize(returns_capacity = mean(returns_capacity,na.rm=TRUE))

ggplot()+
  stat_pointinterval(data=predictions,aes(x = year, y = returns_total,colour=domain),shape=19,.width=0.66,linewidth=1,size=7)+
  #stat_pointinterval(data=predictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = returns_total,colour=domain),shape=19,.width=0.66,linewidth=2,size=9)+
  #stat_pointinterval(data=predictions%>%filter(cu_name=="Sproat Lake"),aes(x = year, y = returns_total,colour=domain),shape=NA,.width=0.66,linewidth=NA,size=9)+
  scale_colour_manual(values=c('South'='darkorchid3','Fraser'='goldenrod3','North'='seagreen3','observed'='grey29'))+
  geom_point(data=adult_data,aes(x=year,y=returns_total,colour=domain),shape=18,size=3)+
  #geom_point(data=adult_data%>%filter(cu_name=="Sproat Lake"),aes(x=year,y=returns_total,colour=domain),shape=18,size=5)+
  geom_hline(data=c_returns_x,aes(yintercept=returns_capacity),linetype='dashed',lwd=0.8,colour='red3')+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme_bw(base_size=10)+
  #theme_bw(base_size=14)+
  labs(title=NULL,x=NULL,y="Returns",colour=NULL)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  facet_wrap(~cu_name, scales = "free_y", ncol=3)+
  theme(legend.position='bottom')

#Simulation####
##model fit on simulated data
fit <- readRDS(file="./output/sockeye_climate_model_MVN_env_sim3_20260630.rds")

b_FTwre<-spread_draws(fit, b_FTwre[p,w])
b_FTsre<-spread_draws(fit, b_FTsre[p,s])
b_TEMPcs<-spread_draws(fit, b_TEMPcs[p])
b_MLDcs<-spread_draws(fit, b_MLDcs[p])
b_TEMPwoo<-spread_draws(fit, b_TEMPwoo[p,w])
b_TEMPsoo<-spread_draws(fit, b_TEMPsoo[p,s])
b_FTum<-spread_draws(fit, b_FTum[p])
b_FDum<-spread_draws(fit, b_FDum[p])

b_FTwre<-b_FTwre%>%mutate(variable="FTwre")%>%select(.chain,.iteration,.draw,p,value=b_FTwre,variable)
b_FTsre<-b_FTsre%>%mutate(variable="FTsre")%>%select(.chain,.iteration,.draw,p,value=b_FTsre,variable)
b_TEMPcs<-b_TEMPcs%>%mutate(variable="TEMPcs")%>%select(.chain,.iteration,.draw,p,value=b_TEMPcs,variable)
b_MLDcs<-b_MLDcs%>%mutate(variable="MLDcs")%>%select(.chain,.iteration,.draw,p,value=b_MLDcs,variable)
b_TEMPsoo<-b_TEMPsoo%>%mutate(variable="TEMPsoo")%>%select(.chain,.iteration,.draw,p,value=b_TEMPsoo,variable)
b_TEMPwoo<-b_TEMPwoo%>%mutate(variable="TEMPwoo")%>%select(.chain,.iteration,.draw,p,value=b_TEMPwoo,variable)
b_FTum<-b_FTum%>%mutate(variable="FTum")%>%select(.chain,.iteration,.draw,p,value=b_FTum,variable)
b_FDum<-b_FDum%>%mutate(variable="FDum")%>%select(.chain,.iteration,.draw,p,value=b_FDum,variable)

pdf<-data.frame(rbind(b_FTwre[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_FTsre[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_TEMPcs[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_MLDcs[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_TEMPwoo[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_TEMPsoo[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_FTum[,c(".chain",".iteration",".draw","p","value","variable")],
                      b_FDum[,c(".chain",".iteration",".draw","p","value","variable")]))

#pdf <- pdf|>group_by(p,variable)|>summarize(value = mean(value,na.rm=TRUE))|>ungroup()

pdf$cu_name <- cuid_datx$cu_name[pdf$p]
#pdf <- pdf|>mutate(cu_name=p)
pdf <- pdf%>%mutate(d = case_when(p %in% c("Great Central Lake","Sproat Lake", "Osoyoos Lake", "Wenatchee Lake") ~ "South",
                                  p %in% c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake") ~ "Fraser", 
                                  p  %in% c("Babine Lake", "Tahltan Lake", "Tatsamenie Lake") ~ "North"))
pdf$d <- factor(pdf$d, levels=c("Hyper","South","Fraser","North"))
pdf$variable <- factor(pdf$variable, levels=c("FTwre","FTsre","TEMPcs","MLDcs","TEMPwoo","TEMPsoo","FTum","FDum"))
pdf <- pdf%>%mutate(v = case_when(variable == "FTwre" ~ "Winter Freshwater Temperature",
                                  variable == "FTsre" ~ "Summer Freshwater Temperature",
                                  variable == "TEMPcs" ~ "Coastal Temperature",
                                  variable == "MLDcs" ~ "Coastal Mixed Layer Depth",
                                  variable == "TEMPwoo" ~ "Winter Open Ocean Temperature",
                                  variable == "TEMPsoo" ~ "Summer Open Ocean Temperature",
                                  variable == "FTum" ~ "Return Migration Temperature",
                                  variable == "FDum" ~ "Return Migration Discharge"))
pdf$variable = factor(pdf$v,levels=c("Winter Freshwater Temperature","Summer Freshwater Temperature","Coastal Temperature","Coastal Mixed Layer Depth","Winter Open Ocean Temperature","Summer Open Ocean Temperature","Return Migration Temperature","Return Migration Discharge"))


##simulation parameters
sim_results <- readRDS("./output/sim_results_20260630.rds")
sim.parameters <- sim_results$env_slopes

#add preset parameters from the simulation
preset <- sim.parameters|>mutate(cu_name=cuid_datx$cu_name[pop_id],FTwre=rowMeans(across(c(b_FTwre1,b_FTwre2,b_FTwre3))),FTsre=rowMeans(across(c(b_FTsre1,b_FTsre2))),TEMPwoo=rowMeans(across(c(b_TEMPwoo1,b_TEMPwoo2,b_TEMPwoo3))),TEMPsoo=rowMeans(across(c(b_TEMPsoo1,b_TEMPsoo2))))|>
  select(cu_name,FTwre,FTsre,TEMPcs=b_TEMPcs,MLDcs=b_MLDcs,TEMPwoo,TEMPsoo,FTum=b_FTum,FDum=b_FDum)|>
  gather(-cu_name,key=variable,value=preset)|>
  mutate(v = case_when(variable == "FTwre" ~ "Winter Freshwater Temperature",
                       variable == "FTsre" ~ "Summer Freshwater Temperature",
                       variable == "TEMPcs" ~ "Coastal Temperature",
                       variable == "MLDcs" ~ "Coastal Mixed Layer Depth",
                       variable == "TEMPwoo" ~ "Winter Open Ocean Temperature",
                       variable == "TEMPsoo" ~ "Summer Open Ocean Temperature",
                       variable == "FTum" ~ "Return Migration Temperature",
                       variable == "FDum" ~ "Return Migration Discharge"))
preset$variable = factor(preset$v,levels=c("Winter Freshwater Temperature","Summer Freshwater Temperature","Coastal Temperature","Coastal Mixed Layer Depth","Winter Open Ocean Temperature","Summer Open Ocean Temperature","Return Migration Temperature","Return Migration Discharge"))

#Accuracy & Precision
apx <- pdf|>select(cu_name,variable=v,value)|>group_by(cu_name,variable)|>summarize(meanx=mean(value,rm.na=TRUE), varx=var(value))|>ungroup()
#apx <- apx|>left_join(preset|>select(cu_name,variable=v,preset))|>mutate(accx = abs(meanx-preset), prex=1/varx)|>group_by(variable)|>mutate(acc_stnd = (accx-mean(accx))/sd(accx), pre_stnd = (prex-mean(prex))/sd(prex), var_stnd = (varx-mean(varx))/sd(varx))
apx <- apx|>left_join(preset|>select(cu_name,variable=v,preset))|>mutate(accx = abs(meanx-preset), prex=1/varx)|>group_by(variable)|>mutate(acc_stnd = (accx - min(accx, na.rm = TRUE)) / (max(accx, na.rm = TRUE) - min(accx, na.rm = TRUE)), pre_stnd = (prex - min(prex, na.rm = TRUE)) / (max(prex, na.rm = TRUE) - min(prex, na.rm = TRUE)), var_stnd = (varx - min(varx, na.rm = TRUE)) / (max(varx, na.rm = TRUE) - min(varx, na.rm = TRUE)))
apx$variable = factor(apx$variable,levels=c("Winter Freshwater Temperature","Summer Freshwater Temperature","Coastal Temperature","Coastal Mixed Layer Depth","Winter Open Ocean Temperature","Summer Open Ocean Temperature","Return Migration Temperature","Return Migration Discharge"))

#95% interval
credx <- pdf|>select(cu_name,variable,value)|>group_by(cu_name,variable)|>summarize(q025 = quantile(value,0.025), q975 = quantile(value,0.975))
credx <- credx|>left_join(preset|>select(cu_name,variable,preset))
credx$int <- credx$preset>=credx$q025 & credx$preset<=credx$q975
apx <- apx|>left_join(credx)

ggplot(apx)+
  geom_point(aes(x = variable, y = accx, colour = int, size = varx))+
  labs(title=NULL,x=NULL, y="Accuracy", size="Variance",colour="C95%")+
  scale_size_continuous(range=c(5,9))+
  scale_colour_manual(values=c('TRUE'="dodgerblue1",'FALSE'="indianred1"))+
  theme_bw(base_size = 13)+
  theme(legend.position='right')+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  coord_cartesian(expand = TRUE)

pdfa <- pdf|>select(cu_name,variable,.chain,.iteration,.draw,value)|>left_join(apx|>select(cu_name,variable,int))|>mutate(int=factor(int,levels=c(TRUE,FALSE)))
ggplot()+geom_violin(data=pdfa,aes(x=cu_name, y=value,fill=int), position="dodge", alpha=0.5)+
  scale_fill_manual(values=c('TRUE'="dodgerblue1",'FALSE'="indianred1"))+
  geom_point(data=preset, aes(x=cu_name,y=preset),colour='black',shape=3,size=3,stroke=1.5)+
  labs(title=NULL,x=NULL, y="Parameter",fill="C95%")+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.1))+
  facet_wrap(~variable, ncol=2, scale="free_y")

