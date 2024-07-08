####Visualization of the results###
###Make meta data including all category####
reads.close100.meta0 <- as.data.frame(t(reads.close100.df)) %>%
  dplyr::mutate(ASV = rownames(.))
reads.close100.meta1 <- left_join(reads.close100.meta0, DA.res.df, by = "ASV")
reads.close100.meta2 <- left_join(reads.close100.meta1, ASV_16S.Phylum.RE, by = "ASV")
reads.close100.meta3 <- left_join(reads.close100.meta2, ASV_TaxID, by = "ASV")
reads.close100.meta <- reads.close100.meta3 %>%
  mutate(Category = replace_na(as.character(Category), "Unclassified")) %>%
  dplyr::select(ASV, everything()) 
head(reads.close100.meta)
names(reads.close100.meta)
write.table(reads.close100.meta, "./table/reads.close100.meta.txt", quote = F, sep = "\t", col.names = T, row.names = F, append = F)



##Grouping by Taxa in all samples
Taxa.all.df <- reads.close100.meta %>% 
  group_by(Taxa) %>% 
  summarise_if(is.numeric, list(sum))
Taxa.all.df <- as.data.frame(Taxa.all.df) %>%
  dplyr::select(contains("H209_"), Taxa)
#Taxa.all.df:RA data frame grouped by Taxa

Taxa.all.ga <- Taxa.all.df %>%
  gather(SampleID, RA, -Taxa) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  dplyr::mutate(Taxa = gsub(Taxa, pattern="_",replacement = " ")) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = Fraction.level)) %>%
  transform(Taxa= factor(Taxa, levels = Taxonomy.level))
head(Taxa.all.ga) 
#Taxa.all.ga:SampleID, Taxa, Station, Layer, Fraction, and RA


##Grouping by Category in all samples
Category.all.df <- reads.close100.meta %>% 
  group_by(Category) %>% 
  summarise_if(is.numeric, list(sum))
Category.all.df <- as.data.frame(Category.all.df) %>%
  dplyr::select(contains("H209_"), Category)
#Category.all.df:RA data frame grouped by Category

Category.all.ga <- Category.all.df %>%
  gather(SampleID, RA, -Category) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  dplyr::mutate(Category = gsub(Category, pattern="_",replacement = " ")) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = Fraction.level)) %>%
  transform(Category= factor(Category, levels = Category.level))
head(Category.all.ga) 
#Category.all.ga:SampleID, Category, Station, Layer, Fraction, and RA



###Make meta data including only Category SK####
reads.close100.meta.SK <- reads.close100.meta %>% dplyr::filter(Category == "SK")
rownames(reads.close100.meta.SK) = reads.close100.meta.SK$ASV
##Grouping by Taxa
Taxa.SK.df <- reads.close100.meta.SK %>% 
  group_by(Taxa) %>% 
  summarise_if(is.numeric, list(sum))
Taxa.SK.df <- as.data.frame(Taxa.SK.df) %>%
  dplyr::select(contains("H209_"), Taxa)
#Taxa.SK.df:RA data frame grouped by Taxa

Taxa.SK.ga <- Taxa.SK.df %>%
  gather(SampleID, RA, -Taxa) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  dplyr::mutate(Taxa = gsub(Taxa, pattern="_",replacement = " ")) %>%
  unite(Fraction.Taxa, c("Fraction", "Taxa"), sep = "_", remove = F) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = Fraction.level)) %>%
  transform(Taxa= factor(Taxa, levels = Taxonomy.level))
head(Taxa.SK.ga) 
#Taxa.SK.ga:SampleID, Taxa, Station, Layer, Fraction, and RA


###Make meta data including only Category SS####
reads.close100.meta.SS <- reads.close100.meta %>% dplyr::filter(Category == "SS")
rownames(reads.close100.meta.SS) = reads.close100.meta.SS$ASV
##Grouping by Taxa
Taxa.SS.df <- reads.close100.meta.SS %>% 
  group_by(Taxa) %>% 
  summarise_if(is.numeric, list(sum))
Taxa.SS.df <- as.data.frame(Taxa.SS.df) %>%
  dplyr::select(contains("H209_"), Taxa)
#Taxa.SS.df:RA data frame grouped by Taxa

Taxa.SS.ga <- Taxa.SS.df %>%
  gather(SampleID, RA, -Taxa) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  dplyr::mutate(Taxa = gsub(Taxa, pattern="_",replacement = " ")) %>%
  unite(Fraction.Taxa, c("Fraction", "Taxa"), sep = "_", remove = F) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = Fraction.level)) %>%
  transform(Taxa= factor(Taxa, levels = Taxonomy.level))
head(Taxa.SS.ga) 
#Taxa.SS.ga:SampleID, Taxa, Station, Layer, Fraction, and RA


###Make meta data including only Category SKSS####
reads.close100.meta.SKSS <- reads.close100.meta %>% dplyr::filter(Category == "SK/SS")
rownames(reads.close100.meta.SKSS) = reads.close100.meta.SKSS$ASV
##Grouping by Taxa
Taxa.SKSS.df <- reads.close100.meta.SKSS %>% 
  group_by(Taxa) %>% 
  summarise_if(is.numeric, list(sum))
Taxa.SKSS.df <- as.data.frame(Taxa.SKSS.df) %>%
  dplyr::select(contains("H209_"), Taxa)
#Taxa.SKSS.df:RA data frame grouped by Taxa

Taxa.SKSS.ga <- Taxa.SKSS.df %>%
  gather(SampleID, RA, -Taxa) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  dplyr::mutate(Taxa = gsub(Taxa, pattern="_",replacement = " ")) %>%
  unite(Fraction.Taxa, c("Fraction", "Taxa"), sep = "_", remove = F) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = Fraction.level)) %>%
  transform(Taxa= factor(Taxa, levels = Taxonomy.level))
head(Taxa.SKSS.ga) 
#Taxa.SKSS.ga:SampleID, Taxa, Station, Layer, Fraction, and RA


###Make meta data including only Category FL####
reads.close100.meta.FL <- reads.close100.meta %>% dplyr::filter(Category == "FL")
rownames(reads.close100.meta.FL) = reads.close100.meta.FL$ASV
##Grouping by Taxa
Taxa.FL.df <- reads.close100.meta.FL %>% 
  group_by(Taxa) %>% 
  summarise_if(is.numeric, list(sum))
Taxa.FL.df <- as.data.frame(Taxa.FL.df) %>%
  dplyr::select(contains("H209_"), Taxa)
#Taxa.FL.df:RA data frame grouped by Taxa

Taxa.FL.ga <- Taxa.FL.df %>%
  gather(SampleID, RA, -Taxa) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  dplyr::mutate(Taxa = gsub(Taxa, pattern="_",replacement = " ")) %>%
  unite(Fraction.Taxa, c("Fraction", "Taxa"), sep = "_", remove = F) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = Fraction.level)) %>%
  transform(Taxa= factor(Taxa, levels = Taxonomy.level))
head(Taxa.FL.ga) 
#Taxa.FL.ga:SampleID, Taxa, Station, Layer, Fraction, and RA


###Make meta data including only Category FLPA####
reads.close100.meta.FLPA <- reads.close100.meta %>% dplyr::filter(Category == "FL/PA")
rownames(reads.close100.meta.FLPA) = reads.close100.meta.FLPA$ASV
##Grouping by Taxa
Taxa.FLPA.df <- reads.close100.meta.FLPA %>% 
  group_by(Taxa) %>% 
  summarise_if(is.numeric, list(sum))
Taxa.FLPA.df <- as.data.frame(Taxa.FLPA.df) %>%
  dplyr::select(contains("H209_"), Taxa)
#Taxa.FLPA.df:RA data frame grouped by Taxa

Taxa.FLPA.ga <- Taxa.FLPA.df %>%
  gather(SampleID, RA, -Taxa) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  dplyr::mutate(Taxa = gsub(Taxa, pattern="_",replacement = " ")) %>%
  unite(Fraction.Taxa, c("Fraction", "Taxa"), sep = "_", remove = F) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = Fraction.level)) %>%
  transform(Taxa= factor(Taxa, levels = Taxonomy.level))
head(Taxa.FLPA.ga) 
#Taxa.FLPA.ga:SampleID, Taxa, Station, Layer, Fraction, and RA



#####


###Fig.2: Taxonomic compositions of prokaryotic community####
plot_Fig.2 <- ggplot(Taxa.all.ga, aes(x=Station, y=RA, fill=Taxa)) + 
  facet_grid(Layer ~ Fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#373737", linewidth = 1), 
        strip.text.x = element_text(color = "white", face = "bold", size = 12), 
        strip.text.y = element_text(color = "white", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour ="black"))+
  geom_bar(stat = "identity") + 
  #scale_fill_nejm() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")+
  scale_fill_manual(values=c("green", "#8b0000", "red", "yellow", "orange", 
                             "blue", "purple", "cyan", "#f0f0f0", "black"))+
  theme(axis.text.x = element_text(size = 6))+
  xlab("Station")+
  ylab("Relative abundance (%)")+
  NULL

graphics.off()
plot(plot_Fig.2)
graphics.off()
png("./figure/plot_Fig.2.png", height=6, width=6, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.2)
dev.off()

graphics.off()
tiff("./figure/plot_Fig.2.tif", height=6, width=6, units="in", res=300)
par(mfrow=c(1, 1))
plot(plot_Fig.2)
dev.off()



###Fig.3 : Category####
plot_Fig.3 <- ggplot(Category.all.ga, aes(x=Station, y=RA, fill=Category)) + 
  facet_grid(Layer ~ Fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#373737", linewidth = 1), 
        strip.text.x = element_text(color = "white", face = "bold", size = 12), 
        strip.text.y = element_text(color = "white", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour ="black"))+
  geom_bar(stat = "identity") + 
  #scale_fill_nejm() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")+
  scale_fill_manual(values=c("#42ff42", "#4242ff", "#ff4242", "#f9f9f9", "#2f2f2f", "lightgrey"))+
  theme(axis.text.x = element_text(size = 6))+
  xlab("Station")+
  ylab("Relative abundance (%)")+
  NULL

graphics.off()
plot(plot_Fig.3)
graphics.off()
png("./figure/plot_Fig.3.png", height=6, width=6, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.3)
dev.off()

graphics.off()
tiff("./figure/plot_Fig.3.tif", height=6, width=6, units="in", res=300)
par(mfrow=c(1, 1))
plot(plot_Fig.3)
dev.off()


###Fig. 4 (RDA plot)####
p_RDA_map <- ggplot(scores.df, mapping = aes(x = RDA1, y = RDA2)) + 
  xlab(paste0("RDA1 (", format(axes.perc[1], nsmall = 1), "%)")) +
  ylab(paste0("RDA2 (", format(axes.perc[2], nsmall = 1), "%)")) +
  NULL

plot_Fig.4 <- p1_RDA_map +
  theme_classic(base_size = 14) +
  geom_point(data=species.score.df, aes(x=RDA1, y=RDA2, color=Category, fill = Category), 
             shape = 21, stroke = 1, alpha = 0.6, size = 2.5) +
  scale_color_manual(values = c("green", "blue", "#dcdcdc", "#dcdcdc"))+
  scale_fill_manual(values = c("green", "blue", "#dcdcdc", "#dcdcdc"))+
  geom_point(data=sites.score.df, aes(x = RDA1, y = RDA2, shape = Fraction_Layer), 
             size = 4, color = "#dc143c") + 
  geom_text_repel(data=sites.score.df, size = 3.5, colour = "#770B20", aes(label = Station))+
  scale_shape_manual(values = c(19, 15, 10, 7)) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(), 
               data = ExpVar.cor.df, linewidth =2, alpha = 0.3, colour = "black") +
  geom_text(data = ExpVar.cor.df, aes(x = text.posi.RDA1, y = text.posi.RDA2), color = "black", 
            size = 5, fontface = "bold", label = ExpVar.cor.df$Variables) + 
  theme(legend.text = element_text(size = 16, colour ="black")) +
  labs(shape = "Fraction-Layer") + 
  xlim(-1.6, 2.2) +
  ylim(-2.7, 2.0) +
  NULL

graphics.off()
plot(plot_Fig.4)
graphics.off()
png("./figure/plot_Fig.4.png", height=6, width=8, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.4)
dev.off()

graphics.off()
tiff("./figure/plot_Fig.4.tif", height=6, width=8, units="in", res=300)
par(mfrow=c(1, 1))
plot(plot_Fig.4)
dev.off()


###plot_Fig.5####
plot_Fig.5.Alpha <- p1_RDA_map +
  geom_point(data=species.score.Alpha, aes(x=RDA1, y=RDA2), , color="#dc143c", fill = "#dc143c", 
             shape = 21, stroke = 1, alpha = 0.9, size = 2.5) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(), 
               data = ExpVar.cor.df, linewidth =2, alpha = 0.3, colour = "black") +
  geom_text(data = ExpVar.cor.df, aes(x = text.posi.RDA1, y = text.posi.RDA2), color = "black", 
            fontface = "bold", label = ExpVar.cor.df$Variables, size = 5) + 
  xlim(-1.6, 2.2) +
  ylim(-2.7, 2.0) +
  annotate("text", x = -0.4, y = -2.2, label="Alphaptoreobacteria",
           family ="serif", fontface="bold", colour="black", size=8, alpha = 0.7) +
  theme_classic(base_size = 16) +
  NULL
graphics.off()
plot(plot_Fig.5.Alpha)
graphics.off()


plot_Fig.5.Gamma <- p1_RDA_map +
  geom_point(data=species.score.Gamma, aes(x=RDA1, y=RDA2), , color="#dc143c", fill = "#dc143c", 
             shape = 21, stroke = 1, alpha = 0.9, size = 2.5) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(), 
               data = ExpVar.cor.df, linewidth =2, alpha = 0.3, colour = "black") +
  geom_text(data = ExpVar.cor.df, aes(x = text.posi.RDA1, y = text.posi.RDA2), color = "black", 
            fontface = "bold", label = ExpVar.cor.df$Variables, size = 5) + 
  xlim(-1.6, 2.2) +
  ylim(-2.7, 2.0) +
  annotate("text", x = -0.35, y = -2.2, label="Gammaptoreobacteria",
           family ="serif", fontface="bold", colour="black", size = 8, alpha = 0.7) +
  theme_classic(base_size = 16) +
  NULL
graphics.off()
plot(plot_Fig.5.Gamma)
graphics.off()


plot_Fig.5.Bac <- p1_RDA_map +
  geom_point(data=species.score.Bac, aes(x=RDA1, y=RDA2), , color="#dc143c", fill = "#dc143c", 
             shape = 21, stroke = 1, alpha = 0.9, size = 2.5) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(), 
               data = ExpVar.cor.df, linewidth =2, alpha = 0.3, colour = "black") +
  geom_text(data = ExpVar.cor.df, aes(x = text.posi.RDA1, y = text.posi.RDA2), color = "black", 
            fontface = "bold", label = ExpVar.cor.df$Variables, size = 5) + 
  xlim(-1.6, 2.2) +
  ylim(-2.7, 2.0) +
  annotate("text", x = -0.5, y = -2.2, label="Bacteroidia",
           family ="serif", fontface="bold", colour="black", size = 8, alpha = 0.7) +
  theme_classic(base_size = 16) +
  NULL
graphics.off()
plot(plot_Fig.5.Bac)
graphics.off()


plot_Fig.5.Plan <- p1_RDA_map +
  geom_point(data=species.score.Plan, aes(x=RDA1, y=RDA2), , color="#dc143c", fill = "#dc143c", 
             shape = 21, stroke = 1, alpha = 0.9, size = 2.5) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(), 
               data = ExpVar.cor.df, linewidth =2, alpha = 0.3, colour = "black") +
  geom_text(data = ExpVar.cor.df, aes(x = text.posi.RDA1, y = text.posi.RDA2), color = "black", 
            fontface = "bold", label = ExpVar.cor.df$Variables, size = 5) + 
  xlim(-1.6, 2.2) +
  ylim(-2.7, 2.0) +
  annotate("text", x = -0.4, y = -2.2, label="Planctomycetota",
           family ="serif", fontface="bold", colour="black", size = 8, alpha = 0.7) +
  theme_classic(base_size = 16) +
  NULL
graphics.off()
plot(plot_Fig.5.Plan)
graphics.off()


plot_Fig.5.Bdell <- p1_RDA_map +
  geom_point(data=species.score.Bdell, aes(x=RDA1, y=RDA2), , color="#dc143c", fill = "#dc143c", 
             shape = 21, stroke = 1, alpha = 0.9, size = 2.5) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(), 
               data = ExpVar.cor.df, linewidth =2, alpha = 0.3, colour = "black") +
  geom_text(data = ExpVar.cor.df, aes(x = text.posi.RDA1, y = text.posi.RDA2), color = "black", 
            fontface = "bold", label = ExpVar.cor.df$Variables, size = 5) + 
  xlim(-1.6, 2.2) +
  ylim(-2.7, 2.0) +
  annotate("text", x = -0.4, y = -2.2, label="Bdellovibrionota",
           family ="serif", fontface="bold", colour="black", size = 8, alpha = 0.7) +
  theme_classic(base_size = 16) +
  NULL
graphics.off()
plot(plot_Fig.5.Bdell)
graphics.off()


plot_Fig.5.Desul <- p1_RDA_map +
  geom_point(data=species.score.Desul, aes(x=RDA1, y=RDA2), , color="#dc143c", fill = "#dc143c", 
             shape = 21, stroke = 1, alpha = 0.9, size = 2.5) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(), 
               data = ExpVar.cor.df, linewidth =2, alpha = 0.3, colour = "black") +
  geom_text(data = ExpVar.cor.df, aes(x = text.posi.RDA1, y = text.posi.RDA2), color = "black", 
            fontface = "bold", label = ExpVar.cor.df$Variables, size = 5) + 
  xlim(-1.6, 2.2) +
  ylim(-2.7, 2.0) +
  annotate("text", x = -0.4, y = -2.2, label="Desulfobacterota",
           family ="serif", fontface="bold", colour="black", size = 8, alpha = 0.7) +
  theme_classic(base_size = 16) +
  NULL
graphics.off()
plot(plot_Fig.5.Desul)
graphics.off()

#
plot_Fig.5 <- ggpubr::ggarrange(plot_Fig.5.Alpha, plot_Fig.5.Gamma, plot_Fig.5.Bac,
                                plot_Fig.5.Bdell, plot_Fig.5.Desul, plot_Fig.5.Plan,
                                 ncol = 3, nrow = 2)

graphics.off()
plot(plot_Fig.5)
graphics.off()
png("./figure/plot_Fig.5.png", height=10, width=16, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.5)
dev.off()

graphics.off()
tiff("./figure/plot_Fig.5.tif", height=10, width=16, units="in", res=300)
par(mfrow=c(1, 1))
plot(plot_Fig.5)
dev.off()


##plot_Fig.6####
head(species.score.df)
species.score.sel <- species.score.df %>%
  dplyr::filter(Taxa == "Alphaproteobacteria" | Taxa == "Gammaproteobacteria" | Taxa == "Bacteroidia" 
                | Taxa == "Planctomycetota" | Taxa == "Bdellovibrionota" | Taxa == "Desulfobacterota") %>%
  transform(Taxa = factor(Taxa, levels = Taxonomy.level))
#Fig.6a
plot_Fig.6a <- ggplot(species.score.sel, aes(x = Taxa, y = RDA1)) +
  geom_boxplot(width=0.5,lwd=1.5) +
  geom_hline(yintercept=0, color="#373737",linewidth=1) +
  theme_light()+
  theme(axis.text.x = element_text(colour = "black", family ="serif", face = "bold", size = 12, 
                                   angle = 60, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 10), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"))+
  labs(x = NULL, y = "RDA1")  +
  NULL


#Fig.6b
plot_Fig.6b <- ggplot(species.score.sel, aes(x = Taxa, y = RDA2)) +
  geom_boxplot(width=0.5,lwd=1.5) +
  geom_hline(yintercept=0, color="#373737",linewidth=1) +
  theme_light()+
  theme(axis.text.x = element_text(colour = "black", family ="serif", face = "bold", size = 12, 
                                   angle = 60, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 10), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"))+
  labs(x = NULL, y = "RDA2")  +
  NULL

#Fig.6
plot_Fig.6 <- ggpubr::ggarrange(plot_Fig.6a, plot_Fig.6b,
                                  ncol = 2,
                                  labels = c("(A)", "(B)"))

graphics.off()
plot(plot_Fig.6)
graphics.off()
png("./figure/plot_Fig.6.png", height=6, width=9, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.6)
dev.off()

graphics.off()
tiff("./figure/plot_Fig.6.tif", height=6, width=9, units="in", res=300)
par(mfrow=c(1, 1))
plot(plot_Fig.6)
dev.off()
#####


###supplementary Figures####
##NMDS
#plot_Fig.S3####
plot_Fig.S3 <- ggplot(data.scoresD, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 5, aes(shape = Fraction, colour = Layer))+
  theme(
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
    axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
    legend.text = element_text(size = 12, face ="bold", colour ="black"), 
    legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
    legend.title = element_text(size = 14, colour = "black", face = "bold"), 
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    legend.key=element_blank()) + 
  geom_text_repel(size = 3, colour = "#373737", aes(label = Station))+
  labs(x = "NMDS1", colour = "Layer", y = "NMDS2", shape = "Fraction")  + 
  scale_colour_manual(values = c("grey", "black"))+
  scale_shape_manual(values = c(19, 17, 8))+
  annotate("text", x = 0.9, y = -0.5, label="Stress = 0.173",
           fontface="bold", colour="black", size=5) +
  NULL
graphics.off()
plot(plot_Fig.S3)
graphics.off()
png("./figure/plot_Fig.S3.png", height=5, width=7, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.S3)
dev.off()

graphics.off()
tiff("./figure/plot_Fig.S3.tif", height=5, width=7, units="in", res=300)
par(mfrow=c(1, 1))
plot(plot_Fig.S3)
dev.off()


##plot_Fig.S4
plot_Fig.S4_D0 <- ggplot(Adiv.LSmeans_D0, aes(x=Layer, y=Least.Sq.Mean, colour = Fraction)) +
  geom_point(size = 4)+
  geom_line(aes(group = Fraction), linewidth = 1.5) +
  scale_colour_manual(values = viridis(3))+
  theme_bw(base_size = 14, base_family = "Arial Rounded MT Bold") +
  geom_errorbar(aes(ymin = low.SD, ymax = high.SD), width=0.1, linewidth = 1.5)+
  labs(x = "Layer", colour = "Fraction", 
       y = "Least square means of \n the effective number of species \n (q = 0)")  + 
  ylim(210, 490) +
  theme(legend.position = "none") +
  NULL

plot_Fig.S4_D1 <- ggplot(Adiv.LSmeans_D1, aes(x=Layer, y=Least.Sq.Mean, colour = Fraction)) + 
  geom_point(size = 4)+
  geom_line(aes(group = Fraction), linewidth = 1.5) +
  scale_colour_manual(values = viridis(3))+
  theme_bw(base_size = 14, base_family = "Arial Rounded MT Bold") +
  geom_errorbar(aes(ymin = low.SD, ymax = high.SD), width=0.1, linewidth = 1.5)+
  labs(x = "Layer", colour = "Fraction", 
       y = "Least square means of \n the effective number of species \n (q = 1)")  + 
  theme(legend.position = "none") +
  ylim(15, 90) +
  NULL

plot_Fig.S4_D2 <- ggplot(Adiv.LSmeans_D2, aes(x=Layer, y=Least.Sq.Mean, colour = Fraction)) + 
  geom_point(size = 4)+
  geom_line(aes(group = Fraction), linewidth = 1.5) +
  scale_colour_manual(values = viridis(3))+
  theme_bw(base_size = 14, base_family = "Arial Rounded MT Bold") +
  geom_errorbar(aes(ymin = low.SD, ymax = high.SD), width=0.1, linewidth = 1.5)+
  labs(x = "Layer", colour = "Fraction",
       y = "Least square means of \n the effective number of species \n (q = 2)")  + 
  theme(legend.position = "none") +
  NULL

plot_Fig.S4_legend.pre <- ggplot(Adiv.LSmeans_D2, aes(x=Layer, y=Least.Sq.Mean, colour = Fraction)) + 
  geom_point(size = 4)+
  geom_line(aes(group = Fraction), linewidth = 1.5) +
  scale_colour_manual(values = viridis(3))+
  theme_bw(base_size = 25)+
  NULL

plot_Fig.S4_legend <- cowplot::get_legend(plot_Fig.S4_legend.pre)

plot_Fig.S4 <- ggpubr::ggarrange(plot_Fig.S4_D0, plot_Fig.S4_D1, plot_Fig.S4_D2, plot_Fig.S4_legend,
                  nrow = 2, ncol = 2,
                  labels = c("(A)", "(B)", "(C)"))

graphics.off()
plot(plot_Fig.S4)
graphics.off()
png("./figure/plot_Fig.S4.png", height=8, width=8, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.S4)
dev.off()

graphics.off()
tiff("./figure/plot_Fig.S4.tif", height=8, width=8, units="in", res=300)
par(mfrow=c(1, 1))
plot(plot_Fig.S4)
dev.off()


####Additional figures####
###Fig.x1: Taxonomic compositions Category FL/PA####
plot_Fig.x1 <- ggplot(Taxa.FLPA.ga, aes(x=Station, y=RA, fill=Taxa)) + 
  facet_grid(Layer ~ Fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#373737", linewidth = 1), 
        strip.text.x = element_text(color = "white", face = "bold", size = 12), 
        strip.text.y = element_text(color = "white", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour ="black"))+
  geom_bar(stat = "identity") + 
  #scale_fill_nejm() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")+
  scale_fill_manual(values=c("green", "red", "yellow", "orange", 
                             "blue", "purple", "cyan", "#f0f0f0", "black"))+
  theme(axis.text.x = element_text(size = 6))+
  xlab("Station")+
  ylab("Relative abundance (%)")+
  NULL
graphics.off()
plot(plot_Fig.x1)
graphics.off()
png("./figure/plot_Fig.x1.png", height=6, width=6, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.x1)
dev.off()


###Fig.x2: Taxonomic compositions Category FL/PA####
plot_Fig.x2 <- ggplot(Taxa.FL.ga, aes(x=Station, y=RA, fill=Taxa)) + 
  facet_grid(Layer ~ Fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#373737", linewidth = 1), 
        strip.text.x = element_text(color = "white", face = "bold", size = 12), 
        strip.text.y = element_text(color = "white", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour ="black"))+
  geom_bar(stat = "identity") + 
  #scale_fill_nejm() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")+
  scale_fill_manual(values=c("green", "#8b0000", "red", "yellow", "orange", 
                             "blue", "#f0f0f0", "black"))+
  theme(axis.text.x = element_text(size = 6))+
  xlab("Station")+
  ylab("Relative abundance (%)")+
  NULL
graphics.off()
plot(plot_Fig.x2)
graphics.off()
png("./figure/plot_Fig.x2.png", height=6, width=6, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.x2)
dev.off()


###Fig.x3: Taxonomic compositions Category SK####
plot_Fig.x3 <- ggplot(Taxa.ga.SK, aes(x=Station, y=RA, fill=Taxa)) + 
  facet_grid(Layer ~ Fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#373737", linewidth = 1), 
        strip.text.x = element_text(color = "white", face = "bold", size = 12), 
        strip.text.y = element_text(color = "white", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        legend.text = element_text(size = 12, colour ="black"))+
  geom_bar(stat = "identity") + 
  #scale_fill_nejm() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")+
  scale_fill_manual(values=c("green", "red", "yellow", "orange", 
                             "blue", "purple", "cyan", "#f0f0f0", "black"))+
  theme(axis.text.x = element_text(size = 6))+
  xlab("Station")+
  ylab("Relative abundance (%)")+
  NULL
graphics.off()
plot(plot_Fig.x3)
graphics.off()
png("./figure/plot_Fig.x3.png", height=6, width=6, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.x3)
dev.off()


###Fig.x4: Taxonomic compositions Category SS####
plot_Fig.x4 <- ggplot(Taxa.ga.SS, aes(x=Station, y=RA, fill=Taxa)) + 
  facet_grid(Layer ~ Fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#373737", linewidth = 1), 
        strip.text.x = element_text(color = "white", face = "bold", size = 12), 
        strip.text.y = element_text(color = "white", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        legend.text = element_text(size = 12, colour ="black"))+
  geom_bar(stat = "identity") + 
  #scale_fill_nejm() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")+
  scale_fill_manual(values=c("red", "yellow", "orange", 
                             "blue", "purple", "cyan", "#f0f0f0", "black"))+
  theme(axis.text.x = element_text(size = 6))+
  xlab("Station")+
  ylab("Relative abundance (%)")+
  NULL
graphics.off()
plot(plot_Fig.x4)
graphics.off()
png("./figure/plot_Fig.x4.png", height=6, width=6, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.x4)
dev.off()


###Fig.x5: Taxonomic compositions Category SK####
plot_Fig.x5 <- ggplot(Taxa.ga.SKSS, aes(x=Station, y=RA, fill=Taxa)) + 
  facet_grid(Layer ~ Fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#373737", linewidth = 1), 
        strip.text.x = element_text(color = "white", face = "bold", size = 12), 
        strip.text.y = element_text(color = "white", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        legend.text = element_text(size = 12, colour ="black"))+
  geom_bar(stat = "identity") + 
  #scale_fill_nejm() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")+
  scale_fill_manual(values=c("green", "red", "yellow", "orange", 
                             "blue", "purple", "cyan", "#f0f0f0", "black"))+
  theme(axis.text.x = element_text(size = 6))+
  xlab("Station")+
  ylab("Relative abundance (%)")+
  NULL
graphics.off()
plot(plot_Fig.x5)
graphics.off()
png("./figure/plot_Fig.x5.png", height=6, width=6, units="in", res=200)
par(mfrow=c(1, 1))
plot(plot_Fig.x5)
dev.off()



