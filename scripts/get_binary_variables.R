mainDir <- "~/Documents/DSN/ALS/" # For running locally
setwd(mainDir)
# Step 1: Read your dictionary file
dict <- read.csv("metadata/Clinical_Dictionary_v3_2024.csv", stringsAsFactors = FALSE)

# Step 2: Standardize and clean the Values column
dict$Values_clean <- trimws(dict$Values)

# Step 3: Define a function to count valid categories
count_categories <- function(values) {
  if (is.na(values) || values == "") return(0)
  parts <- unlist(strsplit(values, "\\|"))
  codes <- sapply(strsplit(parts, ","), function(x) trimws(x[1]))
  return(length(unique(codes)))
}

# Step 4: Apply to dictionary
dict$Num_Categories <- sapply(dict$Values_clean, count_categories)

# Step 5: Filter to binary numeric variables
binary_vars <- subset(dict,
                      Num_Categories == 2 & grepl("numeric", Field.Type, ignore.case = TRUE))

# Step 6: View results
cat("Number of binary numeric variables:", nrow(binary_vars), "\n")
head(binary_vars[, c("Form.Name", "Field.Name", "Field.Label", "Values")])

write.csv(binary_vars, "binary_variables.csv", row.names = FALSE)
