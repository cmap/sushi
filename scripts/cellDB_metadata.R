# Load the required library
library(httr)
library(jsonlite)
library(sets)

# Define the API endpoint URL
api_url <- "https://api.clue.io/api/cell_sets"  

# Define api key 
#api_key <- Sys.getenv("API_KEY")
api_key <- "a0c2e1dab8bcaad34fbb269a3e7c791b"

# Function to get cell lines, pools, sets using API from CellDB
get_cell_api_info = function(api_url, api_key, filter = NULL) {
  
  if (!is.null(filter)){
    # Convert the filter to a JSON string
    filter_json <- toJSON(filter)
    
    # Encode the JSON string for URL query parameter
    filter_encoded <- URLencode(filter_json)
    
    # Create the full API URL with the filter as a query parameter
    api_url <- paste0(api_url, "?filter=", filter_encoded)
  }
  
  # Define headers
  headers <- c("Accept" = "application/json", "user_key" = api_key, "prism_key" = "prism_mts")
  
  # Make the API request with headers
  response <- GET(api_url, add_headers(headers))
  
  # Check if the request was successful (status code 200)
  if (response$status_code == 200) {
    
    # Parse the JSON response into a list
    response_json <- fromJSON(content(response, "text"), flatten = TRUE)
    
    # Convert the list to a dataframe (assuming it's a list of records)
    dataframe <- as.data.frame(response_json)
    return(dataframe)
  } else {
    cat("API request failed with status code:", http_status(response)$status_code, "\n")
    cat("Response content:\n", content(response, "text"), "\n")
  }
}

get_LUAs_from_sets <- function(cell_set_name) {
  
  cell_set_definition_api_url <- "https://api.clue.io/api/cell_set_definition_files"
  # Construct the filter object as a JSON string
  filter <- list(where = list(davepool_id = cell_set_name),
                 fields = c("davepool_id","barcode_id")
  )
  cell_set_lua_df <- get_cell_api_info(cell_set_definition_api_url, api_key, filter)
  LUAs <- list(cell_set_lua_df$barcode_id)
  return(LUAs)
}

get_LUAs_from_pools <- function(cell_pool_name) {
  v_cell_pool_api_url <- "https://api.clue.io/api/v_cell_pools"
  # Construct the filter object as a JSON string
  filter <- list(where = list(pool_id = cell_pool_name), fields = c("pool_id","lua"))
  cell_pool_lua_df <- get_cell_api_info(v_cell_pool_api_url, api_key, filter)
  LUAs <- list(cell_pool_lua_df$lua)
  return(LUAs)
}

get_cell_line_info <- function(all_LUAs) {
  cell_lines_api_url <- "https://api.clue.io/api/v_cell_lines"
  # Construct the filter object as a JSON string
  filter <- list(where = list(lua = list(inq = all_LUAs)))
  cell_lines_info <- get_cell_api_info(cell_lines_api_url, api_key, filter)
  return(cell_lines_info)
}


  
# Pulling/storing full collection of sets, pools, and lines
cell_sets_df <- get_cell_api_info("https://api.clue.io/api/cell_sets", "a0c2e1dab8bcaad34fbb269a3e7c791b")
# uber_pools_df <- get_cell_api_info("https://api.clue.io/api/uber_pools", "a0c2e1dab8bcaad34fbb269a3e7c791b")
cell_pools_df <- get_cell_api_info("https://api.clue.io/api/cell_pools", "a0c2e1dab8bcaad34fbb269a3e7c791b")
cell_lines_df <- get_cell_api_info("https://api.clue.io/api/cell_lines", "a0c2e1dab8bcaad34fbb269a3e7c791b")

sample_meta <- read.csv("/Users/naim/Documents/Work/Troubleshooting/SUSHI_Testing/pseq-005-miseq/sample_meta_test.csv")
raw_counts <- read.csv("/Users/naim/Documents/Work/Troubleshooting/SUSHI_Testing/pseq-005-miseq/raw_counts.csv")

# Creating cell_set_meta
unique_cell_sets <- unique(sample_meta$cell_set)
cell_set_meta <- data.frame(matrix(ncol = 2, nrow = length(unique_cell_sets)))
columns <- c("cell_set", "members")
colnames(cell_set_meta) <- columns

for (i in 1:length(unique_cell_sets)) {
  
  # Splitting cell set by LUA
  chr_unique_cell_sets <- toString(unique_cell_sets[i])
  print(paste("Expanding", unique_cell_sets[i])) 
  cs <- unlist(strsplit(chr_unique_cell_sets,";"))

  # Storing information to generate cell_set_meta
  insert_index <- i
  known_cell_sets <- list()
  all_LUAs <- list()
  
  for (j in 1:length(cs)) {
    if (cs[j] %in% cell_sets_df$name) {
    print(paste(cs[j], "is a cell set.")) 
    known_cell_sets = c(known_cell_sets, cs[j])
      
    # Collecting and storing set LUA members   
    cs_members = get_LUAs_from_sets(cs[j])
    all_LUAs = c(all_LUAs, cs_members)

  } else if (cs[j] %in% cell_pools_df$name){
    print(paste(cs[j], "is a cell pool")) 
    known_cell_sets = c(known_cell_sets, cs[j])
    
    # Collecting pool LUA members
    pool_members = get_LUAs_from_pools(cs[j])
    all_LUAs = c(all_LUAs, pool_members)

  } else if (cs[j] %in% cell_lines_df$lua){
    print(paste(cs, "is a cell line")) 
    all_LUAs = c(all_LUAs, cs[j])
    known_cell_sets = c(known_cell_sets, cs[j])
    
  } else {
    print(paste(cs[j], "is not in our collection of cell sets, pools, or lines")) 
  }}
  
  # Should duplicates be removed before adding to cell_set_meta?
  # How should situations where sets/pools/lines within a collection do not exist in the CellDB be handled?
    # Currently, known cell instances are concatenated together for unique cell sets
    # Should it instead throw warnings, not include the unique cell set entry entirely, etc?
  if (length(all_LUAs) > 0) {
    known_cell_sets <- paste(known_cell_sets, collapse = ";")
    all_LUAs <- paste(all_LUAs, collapse = ";")
    cell_set_meta$cell_set[insert_index] <- known_cell_sets
    cell_set_meta$members[insert_index] <- all_LUAs   
  }}


# cellDBexpansion(sample_meta, raw_counts)

