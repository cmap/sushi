#' cellDB_metadata
#' 
#' takes the user's API information to pull cell set information from CellDB
#' based on the the given project metadata, a curated project cell_set_meta is created
#'
#' @param sample_meta - master metadata of cell lines
#' @param api_key - personal api_key to Clue
#' @return - list with the following elements
#' #' \itemize{
#'   \item cell_set_meta: metadata of cell sets specific to provided project sample metadata
#'   \item cell_line_meta: metadata of cell lines 
#' }
#' @export 
#' 

# Function to get cell lines, pools, sets using API from CellDB
get_cell_api_info = function(api_url, api_key, filter = NULL) {
  
  if (!is.null(filter)){
    # Convert the filter to a JSON string
    filter_json <- toJSON(filter)
    filter_encoded <- URLencode(filter_json)
    api_url <- paste0(api_url, "?filter=", filter_encoded)
  }
  
  # Define headers
  headers <- c("Accept" = "application/json", "user_key" = api_key, "prism_key" = "prism_mts")
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
  
  cell_set_definition_api_url <- "https://api.clue.io/api/e_cell_set_definition_files"
  # Construct the filter object as a JSON string
  filter <- list(where = list(davepool_id = cell_set_name),
                 fields = c("davepool_id","depmap_id","barcode_id")
  )
  cell_set_DepMap_df <- get_cell_api_info(cell_set_definition_api_url, api_key, filter)
  depmap_ids <- list(cell_set_DepMap_df$depmap_id)
  luas <- list(cell_set_DepMap_df$barcode_id)
  res <- data.table(
    cell_set = cell_set_name,
    depmap_id = unlist(depmap_ids),
    lua = unlist(luas)
  )
  return(res)
}

get_LUAs_from_pools <- function(cell_pool_name) {
  v_assay_pool_api_url <- "https://api.clue.io/api/v_e_assay_pools"
  filter <- list(where = list(assay_pool = cell_pool_name), fields = c("assay_pool","barcode_id", "depmap_id"))
  cell_pool_DepMap_df <- get_cell_api_info(v_assay_pool_api_url, api_key, filter)
  luas <- list(cell_pool_DepMap_df$barcode_id)
  depmap_ids <- list(cell_pool_DepMap_df$depmap_id)
    res <- data.table(
        cell_set = cell_pool_name,
        depmap_id = unlist(depmap_ids),
        lua = unlist(luas)
    )
  return(res)
}

create_cell_set_meta = function(sample_meta, cell_sets_df, cell_pools_df, cell_line_meta) {
  unique_cell_sets <- unique(sample_meta$cell_set[sample_meta$cell_set != ""])
  cell_set_meta <- data.frame(matrix(ncol = 2, nrow = length(unique_cell_sets)))
  columns <- c("cell_set", "members")
  colnames(cell_set_meta) <- columns
  failed_sets <- list()
  
  for (i in 1:length(unique_cell_sets)) {
    
    # Splitting cell set by LUA/DepMap_ID
    chr_unique_cell_sets <- toString(unique_cell_sets[i])
    print(paste("Expanding", unique_cell_sets[i])) 
    cs <- unlist(strsplit(chr_unique_cell_sets,";"))
    
    # Storing information to generate cell_set_meta
    insert_index <- i
    known_cell_sets <- list()
    all_DepMapIDs <- list()
    all_luas <- list()
    cell_set_failed <- FALSE
    for (j in 1:length(cs)) {
      if (cs[j] %in% cell_sets_df$name) {
        print(paste(cs[j], "is a cell set.")) 
        known_cell_sets = c(known_cell_sets, cs[j])

        # Collecting and storing set DepMap ID members   
        cs_members = get_LUAs_from_sets(cs[j])
        all_DepMapIDs = append(all_DepMapIDs, list(cs_members$depmap_id))
        all_luas = append(all_luas, list(cs_members$lua))

      } else if (cs[j] %in% cell_pools_df$name){
        print(paste(cs[j], "is a cell pool")) 
        known_cell_sets = append(known_cell_sets, cs[j])
        
        # Collecting pool lua members
        pool_members = get_LUAs_from_pools(cs[j])
        all_luas = append(all_luas, pool_members$lua)
        all_DepMapIDs = append(all_DepMapIDs, pool_members$depmap_id)

        # all_DepMapIDs = append(all_DepMapIDs, pool_members)
        print(paste(cs[j], "has a lua length of:", length(pool_members[[1]])))
        
      } else if (cs[j] %in% cell_line_meta$lua){
        print(paste(cs[j], "is a cell line"))
        member = cell_line_meta %>%
          filter(lua == cs[j]) %>%
          select(depmap_id, lua)
        all_DepMapIDs = append(all_DepMapIDs, member$depmap_id)
        all_luas = append(all_luas, member$lua)
        known_cell_sets = c(known_cell_sets, cs[j])
        
      } else {
        print(paste(cs[j], "is not in our collection of cell sets, pools, or lines"))
        print(paste(chr_unique_cell_sets, "will not be added to the cell_set_meta."))
        cell_set_failed <- TRUE
        break
      }}
    
    if (!cell_set_failed) {
      # Should duplicates be removed before adding to cell_set_meta?
      if (length(all_luas) > 0) {
        known_cell_sets <- paste(known_cell_sets, collapse = ";")
        # Deduplicate all_luas and all_DepMapIDs
        all_luas <- unique(all_luas)
        all_DepMapIDs <- unique(all_DepMapIDs)
        joined_luas <- sapply(all_luas, function(row) paste(row, collapse = ";"))
        all_luas <- paste(joined_luas, collapse = ";")
        joined_DepMapIDs <- sapply(all_DepMapIDs, function(row) paste(row, collapse = ";"))
        all_DepMapIDs <- paste(joined_DepMapIDs, collapse = ";")
        cell_set_meta$cell_set[insert_index] <- known_cell_sets
        cell_set_meta$lua_members[insert_index] <- all_luas
        cell_set_meta$depmap_id_members[insert_index] <- all_DepMapIDs
      }
    } else {
      failed_sets <- c(failed_sets, chr_unique_cell_sets)
    }
  }
  
  cell_set_meta <- na.omit(cell_set_meta)
  if (nrow(cell_set_meta) == 0) {
    print("cell_set_meta is empty. One or more pieces of cell information within a cell set in the cell_set column of sample_meta was not found.")
  }
  
  return(list(cell_set_meta, failed_sets))
}
