#' Cell Metadata API Functions
#' 
#' This file contains functions for retrieving and processing cell metadata
#' from CellDB API and creating cell set metadata for PRISM analysis.
#'
#' @author Broad Institute
#' @version 1.0

#' @title Get Cell Information from API
#' @description Retrieves cell information from CellDB API with optional filtering
#' @param api_url The API endpoint URL
#' @param api_key Authentication key for API access
#' @param filter Optional filter criteria as a list
#' @return A data frame containing the API response or NULL if request fails
#' @examples
#' \dontrun{
#' cell_sets <- get_cell_api_info("https://api.clue.io/api/cell-db/cell-sets", api_key)
#' }
get_cell_api_info <- function(api_url, api_key, filter = NULL) {
  # Apply filter if provided
  if (!is.null(filter)) {
    # Convert filter to JSON and encode for URL
    filter_json <- jsonlite::toJSON(filter, auto_unbox = TRUE)
    filter_encoded <- utils::URLencode(filter_json)
    api_url <- paste0(api_url, "?filter=", filter_encoded)
  }
  
  # Set request headers
  headers <- c(
    "Accept" = "application/json", 
    "user_key" = api_key, 
    "prism_key" = "prism_mts"
  )
  
  # Make API request
  tryCatch({
    response <- httr::GET(api_url, httr::add_headers(headers))
    
    # Check if request was successful
    if (httr::status_code(response) == 200) {
      # Parse JSON response
      response_json <- jsonlite::fromJSON(httr::content(response, "text"), flatten = TRUE)
      
      # Convert to data frame
      result_df <- as.data.frame(response_json)
      
      # Return empty data frame if no results
      if (nrow(result_df) == 0) {
        message("API request returned no results")
        return(data.frame())
      }
      
      return(result_df)
    } else {
      message(sprintf("API request failed with status code: %s", httr::http_status(response)$status_code))
      message(sprintf("Response content: %s", httr::content(response, "text")))
      return(NULL)
    }
  }, error = function(e) {
    message(sprintf("Error in API request: %s", e$message))
    return(NULL)
  })
}

#' @title Get DepMap IDs from Cell Set
#' @description Retrieves DepMap IDs associated with a specific cell set
#' @param cell_set_name The name of the cell set
#' @param api_key Authentication key for API access
#' @param api_base_url Base URL for the API
#' @return A list of DepMap IDs or empty list if not found
get_depmap_ids_from_cell_set <- function(cell_set_name, api_key, 
                                         api_base_url = "https://api.clue.io/api") {
  # Construct API URL
  cell_set_definition_api_url <- paste(api_base_url, "cell-db/cell-sets", sep = "/")
  
  # Create filter for specific cell set
  filter <- list(
    where = list(davepool_id = cell_set_name),
    fields = c("davepool_id", "depmap_id")
  )
  
  # Get cell set data
  cell_set_depmap_ids <- get_cell_api_info(cell_set_definition_api_url, api_key, filter)
  
  # Return empty list if no results
  if (is.null(cell_set_depmap_ids) || nrow(cell_set_depmap_ids) == 0) {
    message(sprintf("No DepMap IDs found for cell set: %s", cell_set_name))
    return(list())
  }
  
  return(list(cell_set_depmap_ids$depmap_id))
}

#' @title Get DepMap IDs from Cell Pool
#' @description Retrieves DepMap IDs associated with a specific cell pool
#' @param cell_pool_name The name of the cell pool
#' @param api_key Authentication key for API access
#' @param api_base_url Base URL for the API
#' @return A list of DepMap IDs or empty list if not found
get_depmap_ids_from_pool <- function(cell_pool_name, api_key, 
                                     api_base_url = "https://api.clue.io/api") {
  # Construct API URL
  assay_pool_api_url <- paste(api_base_url, "cell-db/assay-pools", sep = "/")
  
  # Create filter for specific pool
  filter <- list(
    where = list(assay_pool = cell_pool_name), 
    fields = c("assay_pool", "depmap_id")
  )
  
  # Get pool data
  cell_pool_depmap_ids <- get_cell_api_info(assay_pool_api_url, api_key, filter)
  
  # Return empty list if no results
  if (is.null(cell_pool_depmap_ids) || nrow(cell_pool_depmap_ids) == 0) {
    message(sprintf("No DepMap IDs found for cell pool: %s", cell_pool_name))
    return(list())
  }
  
  return(list(cell_pool_depmap_ids$depmap_id))
}

#' @title Create Cell Set Metadata
#' @description Creates metadata for cell sets by expanding cell set definitions
#' @param sample_meta Sample metadata containing cell set information
#' @param cell_sets_df Data frame of cell sets
#' @param cell_pools_df Data frame of cell pools
#' @param cell_line_meta Data frame of cell line metadata
#' @param api_key Authentication key for API access
#' @param api_base_url Base URL for the API
#' @param verbose Whether to print detailed progress messages
#' @return A list containing cell set metadata and any failed cell sets
create_cell_set_meta <- function(sample_meta, cell_sets_df, cell_pools_df, cell_line_meta,
                                api_key, api_base_url = "https://api.clue.io/api",
                                verbose = TRUE) {
  # Extract unique cell sets from sample metadata
  unique_cell_sets <- unique(sample_meta$cell_set[sample_meta$cell_set != ""])
  
  if (length(unique_cell_sets) == 0) {
    message("No cell sets found in sample metadata")
    return(list(
      cell_set_meta = data.frame(cell_set = character(), members = character()),
      failed_sets = list()
    ))
  }
  
  if (verbose) {
    message(sprintf("Processing %d unique cell sets", length(unique_cell_sets)))
  }
  
  # Initialize result data frame
  cell_set_meta <- data.frame(
    cell_set = character(length(unique_cell_sets)),
    members = character(length(unique_cell_sets)),
    stringsAsFactors = FALSE
  )
  
  # Track failed cell sets
  failed_sets <- list()
  
  # Process each cell set
  for (i in seq_along(unique_cell_sets)) {
    current_cell_set <- as.character(unique_cell_sets[i])
    
    if (verbose) {
      message(sprintf("Processing cell set %d/%d: %s", 
                     i, length(unique_cell_sets), current_cell_set))
    }
    
    # Split cell set by semicolon to get components
    components <- unlist(strsplit(current_cell_set, ";"))
    
    # Track cell groups and their DepMap IDs
    known_cell_groups <- character()
    all_depmap_ids <- list()
    cell_group_failed <- FALSE
    
    # Process each component
    for (component in components) {
      if (component %in% cell_sets_df$davepool_id) {
        # Component is a cell set
        if (verbose) {
          message(sprintf("  - %s is a cell set", component))
        }
        
        known_cell_groups <- c(known_cell_groups, component)
        
        # Get DepMap IDs for this cell set
        component_members <- get_depmap_ids_from_cell_set(component, api_key, api_base_url)
        
        if (length(component_members) > 0 && length(component_members[[1]]) > 0) {
          all_depmap_ids <- append(all_depmap_ids, component_members)
          
          if (verbose) {
            message(sprintf("    Found %d DepMap IDs", length(component_members[[1]])))
          }
        } else {
          message(sprintf("    No DepMap IDs found for cell set: %s", component))
        }
        
      } else if (component %in% cell_pools_df$cell_pool) {
        # Component is a cell pool
        if (verbose) {
          message(sprintf("  - %s is a cell pool", component))
        }
        
        known_cell_groups <- c(known_cell_groups, component)
        
        # Get DepMap IDs for this pool
        pool_members <- get_depmap_ids_from_pool(component, api_key, api_base_url)
        
        if (length(pool_members) > 0 && length(pool_members[[1]]) > 0) {
          all_depmap_ids <- append(all_depmap_ids, pool_members)
          
          if (verbose) {
            message(sprintf("    Found %d DepMap IDs", length(pool_members[[1]])))
          }
        } else {
          message(sprintf("    No DepMap IDs found for cell pool: %s", component))
        }
        
      } else if (component %in% cell_line_meta$depmap_id) {
        # Component is a cell line
        if (verbose) {
          message(sprintf("  - %s is a cell line", component))
        }
        
        all_depmap_ids <- append(all_depmap_ids, component)
        known_cell_groups <- c(known_cell_groups, component)
        
      } else {
        # Component not found
        message(sprintf("  - %s is not in our collection of cell sets, pools, or lines", component))
        message(sprintf("  - %s will not be added to the cell_set_meta", current_cell_set))
        cell_group_failed <- TRUE
        break
      }
    }
    
    # Add to results if all components were found
    if (!cell_group_failed && length(all_depmap_ids) > 0) {
      # Join cell groups and DepMap IDs
      known_cell_groups_str <- paste(known_cell_groups, collapse = ";")
      
      # Flatten and join DepMap IDs
      joined_depmap_ids <- sapply(all_depmap_ids, function(row) {
        if (is.character(row) && length(row) == 1) {
          return(row)
        } else {
          return(paste(row, collapse = ";"))
        }
      })
      all_depmap_ids_str <- paste(joined_depmap_ids, collapse = ";")
      
      # Store in result data frame
      cell_set_meta$cell_set[i] <- known_cell_groups_str
      cell_set_meta$members[i] <- all_depmap_ids_str
      
      if (verbose) {
        message(sprintf("  Added cell set with %d components and %d total DepMap IDs", 
                       length(known_cell_groups), 
                       length(unlist(strsplit(all_depmap_ids_str, ";")))))
      }
    } else if (cell_group_failed) {
      failed_sets <- c(failed_sets, current_cell_set)
    }
  }
  
  # Remove NA rows
  cell_set_meta <- stats::na.omit(cell_set_meta)
  
  # Check if result is empty
  if (nrow(cell_set_meta) == 0) {
    message("WARNING: cell_set_meta is empty. One or more pieces of cell information within a cell set in the cell_set column of sample_meta was not found.")
  } else if (verbose) {
    message(sprintf("Created cell set metadata with %d entries", nrow(cell_set_meta)))
  }
  
  # Return results
  return(list(
    cell_set_meta = cell_set_meta,
    failed_sets = failed_sets
  ))
}
