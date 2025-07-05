demultiplex_tile_stats <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/Demultiplex_Tile_Stats.csv")

#filter na
demultiplex_tile_stats <- demultiplex_tile_stats[!is.na(demultiplex_tile_stats$Tile), ]

#average of values in Reads column
mean(demultiplex_tile_stats$Reads, na.rm = TRUE)

range(demultiplex_tile_stats$Reads, na.rm = TRUE)


#add all values in Reads column
sum(demultiplex_tile_stats$Reads, na.rm = TRUE)
