# Luca's scraping code

library(tidyverse)
library(httr)
library(jsonlite)
library(anytime)

# --- REPLACE IF NECESSARY ---
# you can get the poll name from the URL (https://civiqs.com/results/pollname)
# example poll names: "coronavirus_concern", "coronavirus_response", coronavirus_response_local"

# loop over pollnames

polls <- c("coronavirus_concern","coronavirus_response_local","approve_president_trump","economy_family_retro")

over_polls <- lapply(polls, function(pollname) {
  # --- REPLACE IF NECESSARY ---
  
  # link to get a run ID as preparation for the query. This ensures that the results are fresh!
  url_metadata <- paste0("https://results-api.civiqs.com/results_api/results/", pollname, "/metadata")
  # link to the national data
  url_national <- paste0("https://results-api.civiqs.com/results_api/results/", pollname, "/trendline/%7B%7D?run_id=")
  
  # --- fetching national data ---
  # send a GET request to the API server for a run ID
  response <- GET(url_metadata,
                  accept_json(),
                  add_headers(accept = "application/vnd.questionator.v3"))
  run_id <- fromJSON(content(response, "text"))$run_id
  
  # send a GET request to the API server for polling data
  response <- GET(paste0(url_national, run_id),
                  accept_json(),
                  add_headers(accept = "application/vnd.questionator.v3"))
  # get date information
  date <- fromJSON(content(response, "text"))$trendline %>%
    names() %>%
    lapply(as.numeric) %>%
    lapply(function(x) {x/1000}) %>%
    unlist() %>%
    lapply(anydate) %>%
    lapply(as.character) %>%
    unlist()
  result <- fromJSON(content(response, "text")) %>%
    map_if(is.data.frame, list) %>%
    as_tibble() %>%
    apply(1, unlist) %>%
    t()
  result <- cbind(date, result)
  # standardize column names
  colnames(result) <- colnames(result) %>%
    lapply(tolower) %>%
    gsub("[. ]", "_", .)
  result <- as.data.frame(result)
  
  # edit this path if necessary
  write_csv(result, paste0(pollname, "_national_", Sys.Date(), ".csv"))
  
  # --- fetching state data ---
  # states tracked by Civiqs
  states <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
              "Connecticut", "Delaware", "Florida", "Georgia", "Hawaii", "Idaho",
              "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana",
              "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
              "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
              "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota",
              "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island",
              "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah",
              "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming")
  
  
  id_response <- GET(url_metadata,
                     accept_json(),
                     add_headers(accept = "application/vnd.questionator.v3"))
  run_id <- fromJSON(content(id_response, "text"))$run_id
  
  # build url for state-specific data
  url_state_start <- paste0("https://results-api.civiqs.com/results_api/results/", pollname, "/trendline/%7B%22home_state%22%3A%22")
  url_state_end <- "%22%7D?run_id="
  result_state <- data.frame()
  for (state in states) {
    url <- paste0(url_state_start, gsub(" ", "%20", state), url_state_end, run_id)
    response <- GET(url, accept_json(), add_headers(accept = "application/vnd.questionator.v3"))
    df <- data.frame()
    df <- fromJSON(content(response, "text")) %>%
      map_if(is.data.frame, list) %>%
      as_tibble() %>%
      apply(1, unlist) %>%
      t()
    state_col <- unlist(lapply(date, function(x){state}))
    df <- cbind(state_col, df)
    df <- cbind(date, df)
    colnames(df) <- colnames(df) %>%
      lapply(tolower) %>%
      gsub("[. ]", "_", .)
    df <- as.data.frame(df)
    result_state <- rbind(result_state, df)
  }
  
  write_csv(result_state,paste0("../data/simulation/Civiqs/",pollname,".csv"))

})

