# Low vs. Mock
Day1_LOWvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day1_LOW", "Day1_MOCK"))
Day1_LOWvsMOCK_df <- as.data.frame(Day1_LOWvsMOCK)
Day1_LOWvsMOCK_df$Gene <- rownames(Day1_LOWvsMOCK_df)
Day1_LOWvsMOCK_df$Day <- "Day1"
Day1_LOWvsMOCK_df$Condition <- "Low"

Day2_LOWvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day2_LOW", "Day2_MOCK"))
Day2_LOWvsMOCK_df <- as.data.frame(Day2_LOWvsMOCK)
Day2_LOWvsMOCK_df$Gene <- rownames(Day2_LOWvsMOCK_df)
Day2_LOWvsMOCK_df$Day <- "Day2"
Day2_LOWvsMOCK_df$Condition <- "Low"

Day3_LOWvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day3_LOW", "Day3_MOCK"))
Day3_LOWvsMOCK_df <- as.data.frame(Day3_LOWvsMOCK)
Day3_LOWvsMOCK_df$Gene <- rownames(Day3_LOWvsMOCK_df)
Day3_LOWvsMOCK_df$Day <- "Day3"
Day3_LOWvsMOCK_df$Condition <- "Low"

Day4_LOWvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day4_LOW", "Day4_MOCK"))
Day4_LOWvsMOCK_df <- as.data.frame(Day4_LOWvsMOCK)
Day4_LOWvsMOCK_df$Gene <- rownames(Day4_LOWvsMOCK_df)
Day4_LOWvsMOCK_df$Day <- "Day4"
Day4_LOWvsMOCK_df$Condition <- "Low"

Day5_LOWvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day5_LOW", "Day5_MOCK"))
Day5_LOWvsMOCK_df <- as.data.frame(Day5_LOWvsMOCK)
Day5_LOWvsMOCK_df$Gene <- rownames(Day5_LOWvsMOCK_df)
Day5_LOWvsMOCK_df$Day <- "Day5"
Day5_LOWvsMOCK_df$Condition <- "Low"


# High vs. Mock
Day1_HIGHvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day1_HIGH", "Day1_MOCK"))
Day1_HIGHvsMOCK_df <- as.data.frame(Day1_HIGHvsMOCK)
Day1_HIGHvsMOCK_df$Gene <- rownames(Day1_HIGHvsMOCK_df)
Day1_HIGHvsMOCK_df$Day <- "Day1"
Day1_HIGHvsMOCK_df$Condition <- "High"


Day2_HIGHvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day2_HIGH", "Day2_MOCK"))
Day2_HIGHvsMOCK_df <- as.data.frame(Day2_HIGHvsMOCK)
Day2_HIGHvsMOCK_df$Gene <- rownames(Day2_HIGHvsMOCK_df)
Day2_HIGHvsMOCK_df$Day <- "Day2"
Day2_HIGHvsMOCK_df$Condition <- "High"


Day3_HIGHvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day3_HIGH", "Day3_MOCK"))
Day3_HIGHvsMOCK_df <- as.data.frame(Day3_HIGHvsMOCK)
Day3_HIGHvsMOCK_df$Gene <- rownames(Day3_HIGHvsMOCK_df)
Day3_HIGHvsMOCK_df$Day <- "Day3"
Day3_HIGHvsMOCK_df$Condition <- "High"


Day4_HIGHvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day4_HIGH", "Day4_MOCK"))
Day4_HIGHvsMOCK_df <- as.data.frame(Day4_HIGHvsMOCK)
Day4_HIGHvsMOCK_df$Gene <- rownames(Day4_HIGHvsMOCK_df)
Day4_HIGHvsMOCK_df$Day <- "Day4"
Day4_HIGHvsMOCK_df$Condition <- "High"


Day5_HIGHvsMOCK <- results(DE_1, alpha = 0.05, contrast = c("Collection_Date", "Day5_HIGH", "Day5_MOCK"))
Day5_HIGHvsMOCK_df <- as.data.frame(Day5_HIGHvsMOCK)
Day5_HIGHvsMOCK_df$Gene <- rownames(Day5_HIGHvsMOCK_df)
Day5_HIGHvsMOCK_df$Day <- "Day5"
Day5_HIGHvsMOCK_df$Condition <- "High"
