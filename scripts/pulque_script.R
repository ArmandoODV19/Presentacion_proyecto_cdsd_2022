change_bin_name(home/user01/05.Bins, "pulque")
add_names_to_seqs <- function(nombre_del_archivo){
  filenames <- unlist(strsplit(nombre_del_archivo, "/"))
  filenames <- filenames[[grep("fa", filenames)]]
  divide <- unlist(strsplit(filenames, "\\."))
  bin_name <- divide[1]
  termination <- divide[2]
  old_name <- get.fasta.name(nombre_del_archivo)
  new_name <- paste0( bin_name, "-scaffold-", old_name)
  ref2 <- data.frame(old_name, new_name)
  out_file <- paste0(bin_name, "_renamed", ".", termination)
  rename.fasta(infile = nombre_del_archivo, ref_table = ref2, outfile = out_file)
}
files <- list.files(".")
files
files <- paste0("/home/user01/05.Bins/pulque_Genoma/", files)
map(files, add_names_to_seqs)
files
####
change_bin_name<-function(ruta){
  ruta_original<-getwd()
  #setwd(ruta)
  filez <- list.files()
  file.rename(from=filez, to=sub(pattern="_renamed", replacement="", filez))
  setwd(ruta_original)
}
change_bin_name("/home/user01/05.Bins/pulque_Genoma/bins_named")
###
GTDBK<-read.table("/home/megadriveCURSO/pulqueGTDBTK/classify/gtdbtk.bac120.summary.tsv",
                  sep = "\t", header = T,
                  na.strings ="", stringsAsFactors= F)%>%
  as_tibble()
pulque_gtdbtk<-GTDBK %>%
  select(user_genome, classification) %>%
  separate(classification, c("Domain", "Phylum", "Class", "Order",
                             "Family", "Genus", "Species"), sep= ";") %>%
  rename(Bin_name=user_genome)  %>%
  unite(Bin_name_2, c("Bin_name", "Phylum"), remove = FALSE) %>%
  select(Bin_name, Domain, Phylum, Class, Order, Family, Genus,
         Species)
write.table(pulque_gtdbtk, file = "/home/megadriveCURSO/pulque/Metadatos.txt", sep="\t", quote = F,
            row.names = F, col.names = T)
####
GTDBtk<-pulque_gtdbtk %>%
  count(Domain, Phylum) %>%
  rename(Number_of_MAGs = n) %>%
  ggplot(aes(x = Domain,
             y = Number_of_MAGs, fill = Phylum)) +
  geom_bar(stat = "identity", position=position_dodge())+
  theme_minimal()
library(plotly)
GTDBtk_p_fig <- ggplotly(GTDBtk)
GTDBtk_p_fig
####
library(rbims)
library(tidyverse)
pulque_mapp<-read_ko("/home/megadriveCURSO/pulque/08.Kofamscan/") %>%
  mapping_ko()
Overview<-c("Central Metabolism", "Carbon Fixation",
            "Nitrogen Metabolism", "Sulfur Metabolism", "Fermentation",
            "Methane Metabolism")
Energy_metabolisms_pulque<-pulque_mapp %>%
  drop_na(Cycle) %>%
  get_subset_pathway(rbims_pathway, Overview)
plot_bubble(tibble_ko = Energy_metabolisms_pulque,
            x_axis = Bin_name,
            y_axis = Pathway_cycle,
            analysis="KEGG",
            calc="Percentage",
            range_size = c(1,10),
            y_labs=FALSE,
            x_labs=FALSE)
Metadatos<-read_delim("/home/megadriveCURSO/pulque/Metadatos.txt", delim="\t")
plot_bubble(tibble_ko = Energy_metabolisms_pulque,
            x_axis = Bin_name,
            y_axis = Pathway_cycle,
            analysis="KEGG",
            data_experiment = Metadatos,
            calc="Percentage",
            color_character = Class,
            range_size = c(1,10),
            y_labs=FALSE,
            x_labs=FALSE)
Fermentation_pulque<-pulque_mapp %>%
  drop_na(Cycle) %>%
  get_subset_pathway(Cycle, "Fermentation")
plot_heatmap(tibble_ko=Fermentation_pulque,
             y_axis=Genes,
             analysis = "KEGG",
             calc="Binary")
plot_heatmap(tibble_ko=Fermentation_pulque,
             y_axis=Genes,
             data_experiment = Metadatos,
             x_order = Phylum,
             analysis = "KEGG",
             calc="Binary")
?heatmap
tibble_ko
pulque_gtdbtk %>%
  count(Domain, Phylum) %>%
  rename(Number_of_MAGs = n) %>%
  ggplot(aes(x = Domain,
             y = Number_of_MAGs, fill = Phylum)) +
  geom_bar(stat = "identity", position=position_dodge())+
  theme_minimal()
######
interpro_Interpro_profile<-read_interpro(
  data_interpro = "/home/megadriveCURSO/pulque/05.Bins/pulque_Genoma/bins_named/pulque_interpro.tsv",
  database="INTERPRO", profile = T) %>%
  filter(!str_detect(INTERPRO, "-"))
important_INTERPRO<-get_subset_pca(tibble_rbims=interpro_Interpro_profile,
                                   cos2_val=0.97,
                                   analysis="INTERPRO")
plot_heatmap(important_INTERPRO, y_axis=INTERPRO, analysis = "INTERPRO", distance = T)
plot_heatmap(important_INTERPRO, y_axis=INTERPRO, analysis = "INTERPRO", distance = F)
write_metabolism("/home/megadriveCURSO/pulque/05.Bins/pulque_Genoma/bins_named/pulque_interpro.tsv",
                 "/home/megadriveCURSO/pulque/")
