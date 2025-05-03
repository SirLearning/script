library(readxl); library(dplyr)
Watkins_phenotype <- read_excel("Watseq_phenotype_data/Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx", 
                                sheet = "WGIN_Watkins_JIC_CFLN06")

# height data
attach(Watkins_phenotype)
plant_height <- data.frame(country = `COUNTRY 
of origin`,
                           Rep1 = `PH_M_cm-Rep1`, 
                           Rep2 = `PH_M_cm-Rep2`, 
                           Rep3 = `PH_M_cm-Rep3`, 
                           Rep4 = `PH_M_cm-Rep4`)
plant_height <- plant_height[!grepl("\\*", plant_height$Rep1) & 
                             !grepl("\\*", plant_height$Rep2) &
                             !grepl("\\*", plant_height$Rep3) &
                             !grepl("\\*", plant_height$Rep4) &
                             !grepl("\\*", plant_height$country), ]
plant_height <- plant_height[order(row.names(plant_height)), ]
plant_height$Rep1 <- sapply(strsplit(as.character(plant_height$Rep1), "-"), `[`, 1)
plant_height$Rep2 <- sapply(strsplit(as.character(plant_height$Rep2), "-"), `[`, 1)
plant_height$Rep3 <- sapply(strsplit(as.character(plant_height$Rep3), "-"), `[`, 1)
plant_height$Rep4 <- sapply(strsplit(as.character(plant_height$Rep4), "-"), `[`, 1)
plant_height <- plant_height %>% mutate(Rep1 = as.numeric(Rep1)) %>% arrange(Rep1)
plant_height <- plant_height %>% mutate(Rep2 = as.numeric(Rep2))
plant_height <- plant_height %>% mutate(Rep3 = as.numeric(Rep3))
plant_height <- plant_height %>% mutate(Rep4 = as.numeric(Rep4))
plant_height <- plant_height[!is.na(plant_height$country), ]
plant_height <- na.omit(plant_height)

# related data (heading data & height)
attach(Watkins_phenotype)
related_data <- data.frame(country = `COUNTRY 
of origin`,
                           HD = `Hd_dto_days-CFLN06`,
                           PH = `PH_M_cm-CFLN06`)
related_data <- related_data[!grepl("\\*", related_data$country) &
                             !grepl("\\*", related_data$HD) &
                             !grepl("\\*", related_data$PH), ]
related_data <- na.omit(related_data)
related_data <- related_data %>% mutate(HD = as.numeric(HD)) %>% arrange(HD)
related_data <- related_data %>% mutate(PH = as.numeric(PH))
related_data <- related_data[order(as.numeric(row.names(related_data))), ]


